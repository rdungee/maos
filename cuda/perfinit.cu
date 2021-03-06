/*
  Copyright 2009-2018 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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
#include "utils.h"
#include "accphi.h"
#include "cucmat.h"
#include "kernel.h"
#include "cudata.h"
#include "perf.h"

cuarray<int> cuperf_t::nembed;
cuarray<int> cuperf_t::psfsize;
cuarray<Real> cuperf_t::wvls;
cuarray<stream_t> cuperf_t::stream;
cuarray<cufftHandle> cuperf_t::plan;
int cuperf_t::nevl;
curcell cuperf_t::surf;
curcell cuperf_t::opd;
curcell cuperf_t::psfcl;
curcell cuperf_t::psfcl_ngsr;
curcell cuperf_t::opdcov;
curcell cuperf_t::opdcov_ngsr;
curcell cuperf_t::opdmean;
curcell cuperf_t::opdmean_ngsr;
curcell cuperf_t::cc_cl;
curcell cuperf_t::cc_ol;
curcell cuperf_t::coeff;
Real **cuperf_t::ccb_cl=0;
Real **cuperf_t::ccb_ol=0;
pthread_mutex_t cuperf_t::perfmutex=PTHREAD_MUTEX_INITIALIZER;
cuperf_t::~cuperf_t(){
    //The static members are shared across devices. Need to lock mutex before reinitializing.
    lock_t tmp(perfmutex);
    surf=curcell();
    opd=curcell();
    psfcl=curcell();
    psfcl_ngsr=curcell();
    opdcov=curcell();
    opdcov_ngsr=curcell();
    opdmean=curcell();
    opdmean_ngsr=curcell();
    cc_ol=curcell();
    cc_cl=curcell();
    coeff=curcell();
    if(cuperf_t::ccb_cl){
	for(int ievl=0; ievl<nevl; ievl++){
	    free4async(cuperf_t::ccb_cl[ievl]);
	    free4async(cuperf_t::ccb_ol[ievl]);
	}
        free(cuperf_t::ccb_cl); cuperf_t::ccb_cl=0;
	free(cuperf_t::ccb_ol); cuperf_t::ccb_ol=0;
    }
}
/**
   Initialize perfevl
*/
void gpu_perfevl_init(const PARMS_T *parms, APER_T *aper){
    if(!parms->gpu.evl){
	return;
    }
    const int nevl=parms->evl.nevl;
    const int nwvl=parms->evl.nwvl;
    /*The following lives in CPU. */
    if(parms->evl.psfmean || parms->evl.psfhist){
	cuperf_t::nembed =cuarray<int>(nwvl, 1);
	cuperf_t::psfsize=cuarray<int>(nwvl, 1);
	cuperf_t::wvls   =cuarray<Real>(nwvl, 1);
    
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    cuperf_t::nembed[iwvl]=(int)aper->embed->nembed->p[iwvl];
	    cuperf_t::psfsize[iwvl]=parms->evl.psfsize->p[iwvl];
	    cuperf_t::wvls[iwvl]=parms->evl.wvl->p[iwvl];
	}
    }
    /*The following lives in GPU. */
    for(int im=0; im<NGPU; im++){
	gpu_set(im);
	cudata->perf.locs=culoc_t(aper->locs);
	cp2gpu(cudata->perf.amp, aper->amp);
	cp2gpu(cudata->perf.imcc, aper->imcc);
	if(parms->evl.psfmean || parms->evl.psfhist){
	    cudata->perf.embed    = (int**) calloc(nwvl, sizeof(int*));
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		cp2gpu(&cudata->perf.embed[iwvl], aper->embed->embed->p[iwvl]->p, aper->locs->nloc, 1);
	    }
	}
    }/*for igpu */
    for(int ievl=0; ievl<nevl; ievl++){
	gpu_set(cudata_t::evlgpu[ievl]);
	if(!cudata->perf.locs_dm){
	    cudata->perf.locs_dm=cuarray<cuarray<culoc_t> >(nevl, 1);
	}
	cudata->perf.locs_dm[ievl]=cuarray<culoc_t>(parms->ndm,1);
	for(int idm=0; idm<parms->ndm; idm++){
	    loc_t *loc_dm;
	    if(aper->locs_dm && aper->locs_dm->p[ievl+idm*nevl]){
		loc_dm=aper->locs_dm->p[ievl+idm*nevl];
	    }else{
		loc_dm=aper->locs;
	    }
	    cudata->perf.locs_dm[ievl][idm]=culoc_t(loc_dm);
	}
    }
    cuperf_t::stream=cuarray<stream_t>(nevl, 1);
    if(parms->evl.psfmean || parms->evl.psfhist){
	cuperf_t::plan  = cuarray<cufftHandle>(nwvl*nevl,1);
    }
    for(int ievl=0; ievl<nevl; ievl++){
	gpu_set(cudata_t::evlgpu[ievl]);
	//STREAM_NEW(cuperf_t::stream[ievl]);
	//Use stream created per GPU in order to share resource within GPU between different evl dir.
	cuperf_t::stream[ievl]=cudata->perf_stream;
	if(parms->evl.psfmean || parms->evl.psfhist){
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		DO(cufftPlan2d(&cuperf_t::plan[iwvl+nwvl*ievl],cuperf_t::nembed[iwvl],
			       cuperf_t::nembed[iwvl],FFT_T_C2C));
		DO(cufftSetStream(cuperf_t::plan[iwvl+nwvl*ievl], cuperf_t::stream[ievl]));
	    }/*for iwvl */
	}
    }
    cuperf_t::nevl=nevl;
    cuperf_t::opd=curcell(nevl,1);
    cuperf_t::cc_cl=curcell(nevl, 1);
    cuperf_t::cc_ol=curcell(nevl, 1);
    cuperf_t::coeff=curcell(nevl, 1);
    cuperf_t::ccb_ol=(Real**)malloc(sizeof(Real*)*nevl);
    cuperf_t::ccb_cl=(Real**)malloc(sizeof(Real*)*nevl);
    for(int ievl=0; ievl<nevl; ievl++){
	gpu_set(cudata_t::evlgpu[ievl]);
	cuperf_t::ccb_ol[ievl]=(Real*)malloc4async(sizeof(Real)*7);
	cuperf_t::ccb_cl[ievl]=(Real*)malloc4async(sizeof(Real)*7);
	cuperf_t::cc_cl[ievl]=curmat(7,1);
	cuperf_t::cc_ol[ievl]=curmat(7,1);
	cuperf_t::coeff[ievl]=curmat(7,1);
	cuperf_t::opd[ievl]=curmat(aper->locs->nloc, 1);
    }
    if(!parms->sim.evlol){
	if(parms->evl.cov && parms->gpu.psf){
	    cuperf_t::opdcov=curcell(nevl, 1);
	    cuperf_t::opdmean=curcell(nevl, 1);
	    cuperf_t::opdcov_ngsr=curcell(nevl, 1);
	    cuperf_t::opdmean_ngsr=curcell(nevl, 1);
	}
	if(parms->evl.psfmean || parms->evl.psfhist){
	    cuperf_t::psfcl = curcell(nwvl, parms->evl.nevl);
	    cuperf_t::psfcl_ngsr = curcell(nwvl, parms->evl.nevl);
	}
    }
    if(aper->opdadd){
	cuperf_t::surf=curcell(nevl, 1);
	for(int ievl=0; ievl<nevl; ievl++){
	    gpu_set(cudata_t::evlgpu[ievl]);
	    cp2gpu(cuperf_t::surf[ievl], aper->opdadd->p[ievl]);
	}
    }
    gpu_print_mem("perf init");
}
/*
  Initialize simulation data. Seed dependent. Create for the first seed and zero for the next.
*/
void gpu_perfevl_init_sim(const PARMS_T *parms, APER_T *aper){
    const int nevl=parms->evl.nevl;
    const int nwvl=parms->evl.nwvl;
    int nloc=aper->locs->nloc;
    if(!parms->gpu.evl){
	return;
    }
    /*first open loop ones are on every GPU.*/
    if(parms->evl.psfol){
	for(int im=0; im<NGPU; im++){
	    gpu_set(im);
	    if(parms->evl.cov && parms->gpu.psf){ /*do OL opd cov*/
		initzero(cudata->perf.opdcovol, nloc, nloc);
		initzero(cudata->perf.opdmeanol, nloc, 1);
	    }
	    if(parms->evl.psfmean || parms->evl.psfhist){
		if(cudata->perf.psfol){
		    cuzero(cudata->perf.psfol);
		}else{
		    cudata->perf.psfol=curcell(nwvl,1);
		    for(int iwvl=0; iwvl<nwvl; iwvl++){
			cudata->perf.psfol[iwvl]=curmat(cudata->perf.psfsize[iwvl], 
							cudata->perf.psfsize[iwvl]);
		    }
		}
	    }
	}
    }

    if(parms->evl.cov && parms->gpu.psf && !parms->sim.evlol){
	for(int ievl=0; ievl<nevl; ievl++){
	    if(parms->evl.psf->p[ievl]==0){
		continue;
	    }
	    gpu_set(cudata_t::evlgpu[ievl]);
	    if(parms->evl.psfngsr->p[ievl]){
		initzero(cuperf_t::opdcov_ngsr[ievl], nloc,nloc);
		initzero(cuperf_t::opdmean_ngsr[ievl], nloc,1);
	    }
	    if(parms->evl.psfngsr->p[ievl]!=2){
		initzero(cuperf_t::opdcov[ievl],nloc,nloc);
		initzero(cuperf_t::opdmean[ievl],nloc,1);
	    }
	}
    }
	
    if(parms->evl.psfmean || parms->evl.psfhist){
	if(!parms->sim.evlol){
	    for(int ievl=0; ievl<nevl; ievl++){
		if(parms->evl.psf->p[ievl]==0){
		    continue;
		}
		gpu_set(cudata_t::evlgpu[ievl]);
		for(int iwvl=0; iwvl<nwvl; iwvl++){
		    if(parms->evl.psfngsr->p[ievl]){
			initzero(cuperf_t::psfcl_ngsr[iwvl+nwvl*ievl], 
				 cuperf_t::psfsize[iwvl], cuperf_t::psfsize[iwvl]);
		    }
		    if(parms->evl.psfngsr->p[ievl]!=2){
			initzero(cuperf_t::psfcl[iwvl+nwvl*ievl],
				 cuperf_t::psfsize[iwvl], cuperf_t::psfsize[iwvl]);
		    }
		}	
	    }
	}
    }
}
