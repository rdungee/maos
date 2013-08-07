/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
  This file is part of Multithreaded Adaptive Optics Simulator (MAOS).

  MAOS is free software: you can redistribute it and/or modify it under the
  terms of the GNU General Public License as published by the Free Software
  Foundation, either version 3 of the License, or (at your option) any later
  version.

  MAOS is distributed in the hope that it will be useful, but WITHOUT ANY
  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
  A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along with
  MAOS.  If not, see <http://www.gnu.org/licenses/>.
*/
extern "C"
{
#include <cuda.h>
#include "gpu.h"
}
#include "moao.h"
#include "utils.h"
#include "recon.h"
#include "accphi.h"
#include "pcg.h"
#include "cudata.h"
extern int *wfsgpu;
extern int *evlgpu;

namespace cuda_recon{
cumoao_l::cumoao_l(const PARMS_T *parms, MOAO_T *moao, curecon_geom *_grid)
    :cucg_t(parms?parms->fit.maxit:0, parms?parms->recon.warm_restart:0),grid(_grid),
     NW(0),dotNW(0),amap(0),actslave(0),opdfit(0),opdfit2(0),ha(0){
    amap=new cugrid_t(moao->amap);
    if(moao->NW){
	cp2gpu(&NW, moao->NW->p[0]);
	dotNW=curnew(NW->ny, 1);
    }
    if(moao->actslave){
	actslave=new cusp(moao->actslave->p[0], 1);
    }

    dir_t dir0={0,0,INFINITY,0};
    ha=new map_l2d(grid->fmap, &dir0, 1, amap, 1);
    opdfit=curcellnew(1,1,grid->fmap.nx,grid->fmap.ny);
    opdfit2=curcellnew(1,1,grid->fmap.nx,grid->fmap.ny);
}

cumoao_t::cumoao_t(const PARMS_T *parms, MOAO_T *moao, dir_t *dir, int _ndir, curecon_geom *_grid)
    :cumoao_l(parms,moao,_grid),ndir(_ndir){
    hxp=new map_ray*[ndir];
    hap=new map_ray*[ndir];
    for(int idir=0; idir<ndir; idir++){
	hxp[idir]=new map_l2d(grid->fmap, dir+idir, 1, grid->xmap, grid->npsr);
	hap[idir]=new map_l2d(grid->fmap, dir+idir, 1, grid->amap, grid->ndm);
    }
    rhs=curcellnew(1,1,amap->nx,amap->ny);
}
float cumoao_t::moao_solve(curcell **xout, const curcell *xin, const curcell *ain, stream_t &stream){
    static int count=-1; count++;
    for(int idir=0; idir<ndir; idir++){
	opdfit->m->zero(stream);
	hxp[idir]->forward(opdfit->pm, xin->pm, 1.f, NULL, stream);//tomography	
	hap[idir]->forward(opdfit->pm, ain->pm, -1.f, NULL, stream);//minus common DM.
	grid->W01->apply(opdfit2->m->p, opdfit->m->p, opdfit->nx, stream);
	rhs->m->zero(stream);
	ha->backward(opdfit2->pm, rhs->pm, 1, NULL, stream);
	solve(&xout[idir], rhs, stream);
    }
    return 0;
}
void cumoao_l::L(curcell **xout, float beta, const curcell *xin, float alpha, stream_t &stream){
    if(!*xout){
	*xout=curcellnew(1, 1, amap->nx, amap->ny);
    }else{
	curscale((*xout)->m, beta, stream);
    }
    opdfit->m->zero(stream);
    ha->forward(opdfit->pm, xin->pm, 1, NULL, stream);
    grid->W01->apply(opdfit2->m->p, opdfit->m->p, opdfit->nx, stream);
    ha->backward(opdfit2->pm, (*xout)->pm, alpha, NULL, stream);
    if(NW){
	curmv(dotNW->p, 0, NW, xin->m->p, 't', 1, stream);
	curmv((*xout)->m->p, 1, NW, dotNW->p, 'n', alpha, stream);
    }
    if(actslave){
	cuspmul((*xout)->m->p, actslave, xin->m->p, 1,'n', alpha, stream);
    }
}
}//namespace
/*
  Embed and copy DM commands to GPU.
*/
static void gpu_dm2gpu_embed(curmat *dmgpu, dmat *dmcpu, long *embed, int nx, int ny){
    assert(dmcpu->ny==1);
    double *pin=dmcpu->p;
    float *pout=(float*)calloc(nx*ny, sizeof(float));
    for(int i=0; i<dmcpu->nx; i++){
	pout[embed[i]]=pin[i];
    }
    DO(cudaMemcpy(dmgpu->p, pout, nx*ny*sizeof(float), cudaMemcpyHostToDevice));
    free(pout);
}

/**
   Copy MOAO DM commands from CPU to GPU.*/
void gpu_moao_2gpu(SIM_T *simu){
    const PARMS_T *parms=simu->parms;
    const RECON_T *recon=simu->recon;
    if(parms->gpu.moao){
	error("Invalid use\n");
    }
    const int nwfs=parms->nwfs;
    const int nevl=parms->evl.nevl;
    if(parms->gpu.wfs && simu->dm_wfs){
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfs[iwfs].powfs;
	    int imoao=parms->powfs[ipowfs].moao;
	    if(imoao<0) continue;
	    MOAO_T *moao=recon->moao+imoao;
	    gpu_set(wfsgpu[iwfs]);
	    if(parms->fit.square){
		cp2gpu(&cudata->dm_wfs[iwfs]->p, simu->dm_wfs->p[iwfs]);
	    }else{
		gpu_dm2gpu_embed(cudata->dm_wfs[iwfs]->p, simu->dm_wfs->p[iwfs],
				 moao->aembed, moao->amap->nx, moao->amap->ny);
	    }
	}
    }
    if(parms->gpu.evl && simu->dm_evl){
	int imoao=parms->evl.moao;
	MOAO_T *moao=recon->moao+imoao;
	for(int ievl=0; ievl<nevl; ievl++){
	    gpu_set(evlgpu[ievl]);
	    if(parms->fit.square){
		cp2gpu(&cudata->dm_evl[ievl]->p, simu->dm_evl->p[ievl]);
	    }else{
		gpu_dm2gpu_embed(cudata->dm_evl[ievl]->p, simu->dm_evl->p[ievl],
				 moao->aembed, moao->amap->nx, moao->amap->ny);
	    }
	}
    }
}

