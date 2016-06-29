/*
  Copyright 2009-2016 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#include "common.h"
#include "setup_recon.h"
#include "recon_utils.h"
#include "setup_powfs.h"
#include "pywfs.h"
#include "ahst.h"
/**
   \file setup_recon_dm.c 

   Setup grid and ray tracing operators regarding DM. This is independent of
   1) WFS geometry or noise parameters
   2) Tomography
*/

/**
   Setting up PLOC grid, which is a coarse sampled (usually halves the
   subaperture spacing) grid that defines the circular aperture for tomography.*/
static void
setup_recon_ploc(RECON_T *recon, const PARMS_T *parms){
    double dxr=parms->atmr.dx/parms->tomo.pos;/*sampling of ploc */
    if(parms->load.ploc){/*optionally load ploc from the file. see dbg.conf */
	warning("Loading ploc from %s\n",parms->load.ploc);
	recon->ploc=locread("%s", parms->load.ploc);
	if(fabs(recon->ploc->dx-dxr)>dxr*1e-6){
	    warning("Loaded ploc has unexpected sampling of %g, should be %g\n",
		    recon->ploc->dx, dxr);
	}
    }else{ 
	/*
	  Create a circular PLOC with telescope diameter by calling
	  create_metapupil with height of 0. We don't add any guard points. PLOC
	  does not need to follow XLOC in FDPCG.*/
	double guard=parms->tomo.guard*dxr;
	map_t *pmap=0;
	create_metapupil(&pmap, 0, 0, parms->dirs, parms->aper.d,0,dxr,dxr,0,guard,0,0,0,parms->tomo.square);
	info2("PLOC is %ldx%ld, with sampling of %.2fm\n",pmap->nx,pmap->ny,dxr);
	recon->ploc=map2loc(pmap, 0);/*convert map_t to loc_t */
	mapfree(pmap);
    }
    if(parms->save.setup){
	locwrite(recon->ploc, "ploc");
    }
    loc_create_map_npad(recon->ploc, parms->tomo.square?0:1,0,0);
    recon->pmap=recon->ploc->map;
    loc_create_stat(recon->ploc);
    if(parms->recon.misreg_tel2wfs){
	for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	    if(parms->recon.misreg_tel2wfs[iwfs]){
		if(!recon->ploc_tel){
		    recon->ploc_tel=loccellnew(parms->nwfsr, 1);
		}
		recon->ploc_tel->p[iwfs]=loctransform(recon->ploc, parms->recon.misreg_tel2wfs[iwfs]);
	    }
	}
    }
}
static void
setup_recon_gloc(RECON_T *recon, const PARMS_T *parms, const APER_T *aper){
    //Create another set of loc/amp that can be used to build GP. It has points on edge of subapertures
    recon->gloc=loccellnew(parms->npowfs, 1);
    recon->gamp=dcellnew(parms->npowfs, 1);

    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	map_t *map=0;
	double dx=parms->powfs[ipowfs].dx;
	create_metapupil(&map, 0, 0, parms->dirs, parms->aper.d, 0, dx, dx, 0, 0, 0, 0, 0, 0);
	recon->gloc->p[ipowfs]=map2loc(map, 0); 
	mapfree(map);
	//do not use misregistration since this is the model
	recon->gamp->p[ipowfs]=mkamp(recon->gloc->p[ipowfs], aper->ampground, 
				     0,0, parms->aper.d, parms->aper.din);
	loc_reduce(recon->gloc->p[ipowfs], recon->gamp->p[ipowfs], EPS, 1, NULL);
	
	if(parms->dbg.pupmask && parms->powfs[ipowfs].lo){//for NGS WFS only.
	    int iwfs=parms->powfs[ipowfs].wfs->p[0];
	    wfspupmask(parms, recon->gloc->p[ipowfs], recon->gamp->p[ipowfs], iwfs);
	}
    }
}
/**
   Like ploc, but for DM fitting
*/
static void
setup_recon_floc(RECON_T *recon, const PARMS_T *parms){
    double dxr=parms->atmr.dx/parms->fit.pos;/*sampling of floc */
    if(parms->load.floc){
	warning("Loading floc from %s\n", parms->load.floc);
	recon->floc=locread("%s", parms->load.floc);
	if(fabs(recon->floc->dx-dxr)>dxr*1e-6){
	    warning("Loaded floc has unexpected sampling of %g, should be %g\n", 
		    recon->floc->dx, dxr);
	}
    }else{
	double guard=parms->tomo.guard*dxr;
	map_t *fmap=0;
	create_metapupil(&fmap,0,0,parms->dirs,parms->aper.d,0,dxr,dxr,0,guard,0,0,0,parms->fit.square);
	info2("FLOC is %ldx%ld, with sampling of %.2fm\n",fmap->nx,fmap->ny,dxr);
	recon->floc=map2loc(fmap, 0);/*convert map_t to loc_t */
	mapfree(fmap);
	/*Do not restrict fmap to within active pupil. */
    }
    loc_create_map_npad(recon->floc, parms->fit.square?0:1, 0, 0);
    recon->fmap=recon->floc->map;
    /*create the weighting W for bilinear influence function. See [Ellerbroek 2002] */
    if(parms->load.W){
	if(!(zfexist("W0")&&zfexist("W1"))){
	    error("W0 or W1 not exist\n");
	}
	warning("Loading W0, W1");
	recon->W0=dspread("W0");
	recon->W1=dread("W1");
    }else{
	/*
	  Compute W0,W1 weighting matrix that can be used to compute piston
	  removed wavefront variance of OPD: RMS=OPD'*(W0-W1*W1')*OPD; W0 is
	  sparse. W1 is a vector. These matrices are used for weighting in DM
	  fitting.
	*/
	double rin=0;
	if(parms->dbg.annular_W && parms->aper.din>0){
	    warning("Define the W0/W1 on annular aperture instead of circular.\n");
	    rin=parms->aper.din/2;
	}
	mkw_annular(recon->floc, 0, 0, rin, parms->aper.d/2,
		    &(recon->W0), &(recon->W1));
    }
    if(parms->save.setup){
	writebin(recon->W0, "W0");
	writebin(recon->W1, "W1");
    }
    if(parms->save.setup){
	locwrite(recon->floc, "floc");
    }
    loc_create_stat(recon->floc);
}

/**
   Setup the tomography grids xloc which is used for Tomography.
*/
static void 
setup_recon_xloc(RECON_T *recon, const PARMS_T *parms){
    const int npsr=recon->npsr;
    long nin0=0;
    if(parms->load.xloc){
	char *fn=parms->load.xloc;
	warning("Loading xloc from %s\n",fn);
	recon->xloc=loccellread("%s",fn);
	int nxloc=recon->xloc->nx;
	if(nxloc!=npsr) 
	    error("Invalid saved file. npsr=%d, nxloc=%d\n",npsr,nxloc);
	for(int ips=0; ips<npsr; ips++){
	    double dxr=recon->dx->p[ips];
	    if(fabs(recon->xloc->p[ips]->dx-dxr)>0.01*dxr){
		warning("xloc[%d]: sampling is %g, expected %g\n", ips, recon->xloc->p[ips]->dx, dxr);
	    }
	}
    }else{
	recon->xloc=loccellnew(npsr, 1);
	info2("Tomography grid is %ssquare:\n", parms->tomo.square?"":"not ");
	/*FFT in FDPCG prefers power of 2 dimensions. for embeding and fast FFT*/
	if(parms->tomo.nxbase){
	    nin0=parms->tomo.nxbase;
	}else if(!parms->sim.idealfit && (parms->tomo.precond==1 || parms->tomo.square==2)){
	    /*same square grid dimension in meter on all layers.*/
	    long nxmin=LONG_MAX, nymin=LONG_MAX;
	    long nxmax=0, nymax=0;
	    for(int ips=0; ips<npsr; ips++){
		long nxi, nyi;
		double dxr=recon->dx->p[ips];
		create_metapupil(0, &nxi, &nyi, parms->dirs, parms->aper.d, recon->ht->p[ips], dxr, dxr, 0,
				 dxr*parms->tomo.guard, 0, 0, 0, 1);
		nxi/=recon->os->p[ips];
		nyi/=recon->os->p[ips];
		if(nxmax<nxi) nxmax=nxi;
		if(nymax<nyi) nymax=nyi;
		if(nxmin>nxi) nxmin=nxi;
		if(nymin>nyi) nymin=nyi;
	    }
	    /*FFT grid Must be at least 1.5 times the smallest (on pupil) to avoid
	     * severe aliasing penalty*/
	    long nx=MAX(nxmax, nxmin*1.5);
	    long ny=MAX(nymax, nymin*1.5);
	    nin0=nextfftsize(MAX(nx, ny));
	
	    while (parms->tomo.precond==1 && (nin0 & 1)){//FFT need even number
		nin0=nextfftsize(nin0+1);
	    }
	}
	for(int ips=0; ips<npsr; ips++){
	    const double ht=recon->ht->p[ips];
	    double dxr=(parms->sim.idealfit)?parms->atm.dx:recon->dx->p[ips];
	    const double guard=parms->tomo.guard*dxr;
	    long nin=nin0*recon->os->p[ips];
	    map_t *map=0;
	    create_metapupil(&map, 0, 0, parms->dirs, parms->aper.d,ht,dxr,dxr,0,guard,nin,nin,0,parms->tomo.square);
	    recon->xloc->p[ips]=map2loc(map, 0);
	    loc_create_stat(recon->xloc->p[ips]);
	    info2("layer %d: xloc grid is %3ld x %3ld, sampling is %.3f m, %5ld points\n",
		  ips, map->nx,map->ny,dxr, recon->xloc->p[ips]->nloc);
	    mapfree(map);
	}
    }
    if(parms->gpu.fit==2 && parms->fit.cachex){//to cache x on grid matching floc.
	recon->xcmap=mapcellnew(npsr, 1);
	for(int ips=0; ips<npsr; ips++){
	    const double ht=recon->ht->p[ips];
	    double dxr=parms->atmr.dx/parms->fit.pos;
	    const double guard=parms->tomo.guard*dxr;
	    create_metapupil(&recon->xcmap->p[ips], 0, 0, parms->dirs, parms->aper.d,ht,dxr,dxr,0,guard,0,0,0,parms->fit.square);
	    free(recon->xcmap->p[ips]->p);recon->xcmap->p[ips]->p=NULL;
	    free(recon->xcmap->p[ips]->nref);recon->xcmap->p[ips]->nref=NULL;
	}
    }
    recon->xmap=mapcellnew(npsr, 1);
    recon->xnx=lnew(recon->npsr, 1);
    recon->xny=lnew(recon->npsr, 1);
    recon->xnloc=lnew(recon->npsr, 1);
    for(long i=0; i<recon->npsr; i++){
	recon->xloc->p[i]->iac=parms->tomo.iac;
	loc_create_map_npad(recon->xloc->p[i], (nin0||parms->tomo.square)?0:1, 
			    nin0*recon->os->p[i], nin0*recon->os->p[i]);
	recon->xmap->p[i]=mapref(recon->xloc->p[i]->map);
	recon->xmap->p[i]->h=recon->ht->p[i];
	recon->xnx->p[i]=recon->xmap->p[i]->nx;
	recon->xny->p[i]=recon->xmap->p[i]->ny;
	recon->xnloc->p[i]=recon->xloc->p[i]->nloc;
    }
    recon->xmcc=dcellnew(npsr,1);
    for(int ipsr=0; ipsr<npsr; ipsr++){
	recon->xmcc->p[ipsr]=loc_mcc_ptt(recon->xloc->p[ipsr],NULL);
	dinvspd_inplace(recon->xmcc->p[ipsr]);
    }
    if(parms->save.setup){
	writebin(recon->xloc, "xloc");
	writebin(recon->xmap, "xmap");
    }
}

/**
   Setup the deformable mirrors grid aloc. This is used for DM fitting.
*/
static void
setup_recon_aloc(RECON_T *recon, const PARMS_T *parms){
    const int ndm=parms->ndm;
    if(ndm==0) return;
    if(parms->fit.cachedm){
	recon->acmap=mapcellnew(ndm, 1);
    }
    if(parms->load.aloc){
	char *fn=parms->load.aloc;
	warning("Loading aloc from %s\n",fn);
	recon->aloc=loccellread("%s",fn);
	if(recon->aloc->nx!=ndm || recon->aloc->ny!=1) {
	    error("Loaded aloc should have %dx1 cells\n", ndm);
	}
	for(int idm=0; idm<ndm; idm++){
	    if(fabs(parms->dm[idm].dx-recon->aloc->p[idm]->dx)>1e-7){
		error("DM[%d]: loaded aloc has dx=%g while dm.dx=%g\n", idm, 
		      recon->aloc->p[idm]->dx, parms->dm[idm].dx);
	    }
	    double max,min;
	    dmaxmin(recon->aloc->p[idm]->locx,recon->aloc->p[idm]->nloc, &max, &min);
	    if(max-min<parms->aper.d){
		warning("DM[%d]: loaded aloc is too small: diameter is %g while aper.d is %g\n", 
		      idm, max-min, parms->aper.d); 
	    }
	}
    }else{
	recon->aloc=loccellnew(ndm, 1);
	/*int nxmax=0, nymax=0; */
	for(int idm=0; idm<ndm; idm++){
	    double ht=parms->dm[idm].ht;
	    double dx=parms->dm[idm].dx;
	    double dy=parms->dm[idm].dy;
	    double offset=parms->dm[idm].offset+((int)round(parms->dm[idm].order)%2)*0.5;
	    double guard=parms->dm[idm].guard*MAX(dx,dy);
	    map_t *map;
	    if(parms->dbg.dmfullfov && !parms->fit.square){//DM covers full fov
		double D=(parms->sim.fov*ht+parms->aper.d+guard*2);
		long nx=D/dx+1;
		long ny=D/dy+1;
		map=mapnew(nx, ny, dx, dy, 0);
		map->h=ht;
		map->ox+=offset*dx;
		mapcircle_symbolic(map, D*0.5);
	    }else{
		create_metapupil(&map,0,0,parms->dirs, parms->aper.d,ht,dx,dy,offset,guard,0,0,0,parms->fit.square);
	    }
	    info2("DM %d: grid is %ld x %ld\n", idm, map->nx, map->ny);
	    recon->aloc->p[idm]=map2loc(map, 0);
	    mapfree(map);
	}
    }
    recon->amap=mapcellnew(parms->ndm, 1);
    for(int idm=0; idm<parms->ndm; idm++){
	double ht=parms->dm[idm].ht;
	double offset=parms->dm[idm].offset+((int)round(parms->dm[idm].order)%2)*0.5;
	double dx=parms->dm[idm].dx;
	recon->aloc->p[idm]->iac=parms->dm[idm].iac;
	loc_create_map_npad(recon->aloc->p[idm], parms->fit.square?0:1, 0, 0);
	recon->amap->p[idm]=recon->aloc->p[idm]->map;
	recon->amap->p[idm]->h=ht;
	if(parms->fit.cachedm){
	    const double dx2=parms->atmr.dx/parms->fit.pos;
	    const double dy2=dx2;
	    create_metapupil(&recon->acmap->p[idm],0,0, parms->dirs, parms->aper.d,
			     ht,dx2,dy2,offset*dx/dx2,dx2,0,0,0,parms->fit.square);
	    info2("amap origin is %g, %g. acmap is %g, %g\n", 
		 recon->aloc->p[idm]->map->ox, recon->aloc->p[idm]->map->oy,
		 recon->acmap->p[idm]->ox, recon->acmap->p[idm]->oy);
	}
    }

    recon->aimcc=dcellnew(ndm,1);
    for(int idm=0; idm<ndm; idm++){
	recon->aimcc->p[idm]=loc_mcc_ptt(recon->aloc->p[idm], NULL);
	dinvspd_inplace(recon->aimcc->p[idm]);
    }
    recon->anx=lnew(ndm, 1);
    recon->any=lnew(ndm, 1);
    recon->anloc=lnew(ndm, 1);
    for(int idm=0; idm<ndm; idm++){
	recon->anx->p[idm]=recon->amap->p[idm]->nx;
	recon->any->p[idm]=recon->amap->p[idm]->ny;
	recon->anloc->p[idm]=recon->aloc->p[idm]->nloc;
    }
    /*Dealing with stuck/floating actuators. */
    int anyfloat=0, anystuck=0;
    for(int idm=0; idm<ndm; idm++){
	if(parms->dm[idm].actstuck) anystuck=1;
	if(parms->dm[idm].actfloat) anyfloat=1;
    }
    if(anystuck){
	recon->actstuck=lcellnew(parms->ndm, 1);
	for(int idm=0; idm<ndm; idm++){
	    if(!parms->dm[idm].actstuck) continue;
	    recon->actstuck->p[idm]=loc_coord2ind(recon->aloc->p[idm], parms->dm[idm].actstuck);
	    double stroke=INFINITY;
	    const int nact=recon->aloc->p[idm]->nloc;
	    if(parms->dm[idm].stroke){
		if(parms->dm[idm].stroke->nx==1){
		    stroke=parms->dm[idm].stroke->p[0];
		    dfree(parms->dm[idm].stroke);
		}else if(parms->dm[idm].stroke->nx!=nact){
		    error("dm.stroke is in wrong format\n");
		}
	    }
	    if(!parms->dm[idm].stroke){
		parms->dm[idm].stroke=dnew(nact, 2);
		for(int iact=0; iact<nact; iact++){
		    parms->dm[idm].stroke->p[iact]=-stroke;
		    parms->dm[idm].stroke->p[iact+nact]=stroke;
		}
		((PARMS_T*)parms)->sim.dmclip=1;
	    }
	    for(int iact=0; iact<nact; iact++){
		double val=recon->actstuck->p[idm]->p[iact];
		if(val){
		    parms->dm[idm].stroke->p[iact]=val*1e-9;
		    parms->dm[idm].stroke->p[iact+nact]=val*1e-9;
		}
	    }
	}
    }
    if(anyfloat){
	recon->actfloat=lcellnew(parms->ndm, 1);
	for(int idm=0; idm<ndm; idm++){
	    if(!parms->dm[idm].actfloat) continue;
	    recon->actfloat->p[idm]=loc_coord2ind(recon->aloc->p[idm], parms->dm[idm].actfloat);
	}
    }
 
    if(parms->save.setup){
	writebin(recon->aloc,"aloc");
	writebin(recon->amap, "amap");
    }
}
/**
   Setup ray tracing operator from xloc to ploc for guide stars
*/

static void 
setup_recon_HXW(RECON_T *recon, const PARMS_T *parms){
    loc_t *ploc=recon->ploc;
    const int nwfs=parms->nwfsr;
    const int npsr=recon->npsr;
    if(parms->load.HXW){
	warning2("Loading saved HXW\n");
	recon->HXW=dspcellread("%s",parms->load.HXW);
	if(recon->HXW->nx!=nwfs || recon->HXW->ny!=npsr){
	    error("Wrong saved HXW\n");
	}
	dspcell* HXW=recon->HXW/*PDSPCELL*/;
	int nploc=ploc->nloc;
	for(int ips=0; ips<npsr; ips++){
	    int nloc=recon->xloc->p[ips]->nloc;
	    for(int iwfs=0; iwfs<nwfs; iwfs++){
		if(!IND(HXW,iwfs,ips) 
		   || IND(HXW,iwfs,ips)->nx!=nploc 
		   || IND(HXW,iwfs,ips)->ny!=nloc){
		    error("Wrong saved HXW\n");
		}
	    }
	}
    }else{
	info2("Generating HXW");TIC;tic;
	recon->HXW=dspcellnew(nwfs,npsr);
	dspcell* HXW=recon->HXW/*PDSPCELL*/;
    	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs = parms->wfsr[iwfs].powfs;
	    
	    if(parms->recon.split!=2 && parms->powfs[ipowfs].skip){
		//don't need HXW for low order wfs that does not participate in tomography. 
		continue;
	    }
	    const double  hs = parms->wfs[iwfs].hs;
	    loc_t *loc=recon->ploc;
	    if(recon->ploc_tel && recon->ploc_tel->p[iwfs]){
		loc=recon->ploc_tel->p[iwfs];
	    }
	    for(int ips=0; ips<npsr; ips++){
		const double  ht = recon->ht->p[ips];
		const double  scale=1. - ht/hs;
		const double dispx=parms->wfsr[iwfs].thetax*ht;
		const double dispy=parms->wfsr[iwfs].thetay*ht;
		IND(HXW,iwfs,ips)=mkh(recon->xloc->p[ips], loc, 
				   dispx,dispy,scale);
	    }
	}
	toc2(" ");
    }
    if(parms->save.setup){
	writebin(recon->HXW, "HXW");
    }
    recon->HXWtomo=dspcellnew(recon->HXW->nx, recon->HXW->ny);
    dspcell* HXWtomo=recon->HXWtomo/*PDSPCELL*/;
    dspcell* HXW=recon->HXW/*PDSPCELL*/;
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs=parms->wfsr[iwfs].powfs;
	if(!parms->powfs[ipowfs].skip){/*for tomography */
	    for(int ips=0; ips<npsr; ips++){
		IND(HXWtomo,iwfs,ips)=dspref(IND(HXW,iwfs,ips));
	    }
	} 
    }
}

/**
   Setup ray tracing operator HA from aloc to aperture ploc along DM fiting direction*/
static void
setup_recon_HA(RECON_T *recon, const PARMS_T *parms){
    if(parms->load.HA && zfexist(parms->load.HA)){
	warning("Loading saved HA\n");
	recon->HA=dspcellread("%s",parms->load.HA);
    }else{
	const int nfit=parms->fit.nfit;
	const int ndm=parms->ndm;
	recon->HA=dspcellnew(nfit, ndm);
	dspcell* HA=recon->HA/*PDSPCELL*/;
	info2("Generating HA ");TIC;tic;
	for(int ifit=0; ifit<nfit; ifit++){
	    double hs=parms->fit.hs->p[ifit];
	    for(int idm=0; idm<ndm; idm++){
		const double ht=parms->dm[idm].ht;
		const double scale=1.-ht/hs;
		double displace[2];
		displace[0]=parms->fit.thetax->p[ifit]*ht;
		displace[1]=parms->fit.thetay->p[ifit]*ht;
		loc_t *loc=recon->floc;
		if(parms->recon.misreg_dm2sci && parms->recon.misreg_dm2sci[ifit+idm*nfit]){
		    loc=loctransform(loc, parms->recon.misreg_dm2sci[ifit+idm*nfit]);
		}
		IND(HA,ifit,idm)=mkh(recon->aloc->p[idm], loc, 
				  displace[0], displace[1], scale);
		if(loc!=recon->floc){
		    locfree(loc);
		}
	    }
	}
	toc2(" ");
    }
    if(parms->save.setup){
	writebin(recon->HA,"HA");
    }
    recon->actcpl=genactcpl(recon->HA, recon->W1);
    //if(1){//new 
    //cpl accounts for floating actuators, but not stuck actuators.
    act_stuck(recon->aloc, recon->actcpl, recon->actfloat);
    //Do not modify HA by floating actuators, otherwise, HA*actinterp will not work.
    act_stuck(recon->aloc, recon->HA, recon->actstuck);
    if(parms->save.setup){
	writebin(recon->HA,"HA_float");
    }
    if(parms->fit.actinterp){
	recon->actinterp=act_extrap(recon->aloc, recon->actcpl, parms->fit.actthres);
    }else if(recon->actfloat){
	warning("There are float actuators, but fit.actinterp is off\n");
    }
    if(recon->actinterp){
	/*
	  DM fitting output a is extrapolated to edge actuators by
	  actinterp*a. The corresponding ray tracing from DM would be
	  HA*actinterp*a. We replace HA by HA*actinterp to take this into
	  account during DM fitting.
	 */
	info2("Replacing HA by HA*actinterp\n");
	
	dspcell *HA2=0;
	dcellmm(&HA2, recon->HA, recon->actinterp, "nn", 1);
	dspcellfree(recon->HA);
	recon->HA=HA2;
	if(parms->save.setup){
	    writebin(recon->HA,"HA_final");
	}
    }
    
    if(parms->save.setup){
	if(recon->actinterp){
	    writebin(recon->actinterp, "actinterp");
	}
	if(recon->actcpl){
	    writebin(recon->actcpl, "actcpl");
	}
    }
}
/**
   Setup gradient operator from ploc to wavefront sensors.
*/
static void
setup_recon_GP(RECON_T *recon, const PARMS_T *parms, const POWFS_T *powfs){
    loc_t *ploc=recon->ploc;
    const int nwfs=parms->nwfsr;
    dspcell *GP=NULL;
    if(parms->load.GP){
	warning2("Loading saved GP\n");
	GP=dspcellread("%s",parms->load.GP);
	int nploc=ploc->nloc;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    int nsa=powfs[ipowfs].saloc->nloc;
	    if(GP->p[ipowfs]->nx!=nsa*2 || GP->p[ipowfs]->ny!=nploc){
		error("Wrong saved GP: size is %ldx%ld, need %dx%d\n",
		      GP->p[ipowfs]->nx, GP->p[ipowfs]->ny, nsa*2, nploc);
	    }
	}
    }else{
	GP=dspcellnew(parms->npowfs, 1);
	info2("Generating GP with ");TIC;tic;
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    if(parms->powfs[ipowfs].nwfs==0) continue;
	    /*use ploc as an intermediate plane.  Amplitude must use assumed amp
	      (non-misregistered)*/
	    if(parms->powfs[ipowfs].type==1){
		info2(" PWFS (skip)");
		/*if(parms->recon.alg==0){
		    GP->p[ipowfs]=pywfs_mkg(powfs[ipowfs].pywfs, ploc,0,0,1);
		    }*/
	    }else if(parms->powfs[ipowfs].gtype_recon==0){
		/*Create averaging gradient operator (gtilt) from PLOC,
		  using fine sampled powfs.gloc as intermediate plane*/
		info2(" Gploc");
		GP->p[ipowfs]=mkg(ploc,recon->gloc->p[ipowfs],recon->gamp->p[ipowfs],
				  powfs[ipowfs].saloc,1,0,0,1);
	    }else if(parms->powfs[ipowfs].gtype_recon==1){
		/*Create ztilt operator from PLOC, using fine sampled
		  powfs.gloc as intermediate plane.*/
		dsp* ZS0=mkz(recon->gloc->p[ipowfs],recon->gamp->p[ipowfs]->p,
			     (loc_t*)powfs[ipowfs].pts, 1,1,0,0);
		info2(" Zploc");
		dsp *H=mkh(ploc,recon->gloc->p[ipowfs], 0,0,1);
		GP->p[ipowfs]=dspmulsp(ZS0,H,"nn");
		dspfree(H);
		dspfree(ZS0);
	    }else{
		error("Invalid gtype_recon\n");
	    }
	}
	toc2(" ");
    }
    if(parms->save.setup){
	writebin(GP,"GP");
    }
    /*assign GP for powfs to recon->GP for each wfs */
    recon->GP=GP;
    recon->GP2=dspcellnew(nwfs,1);
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs = parms->wfsr[iwfs].powfs;
	recon->GP2->p[iwfs]=dspref(recon->GP->p[ipowfs]);
    }
}
/**
   Setup gradient operator form aloc for wfs by using GP.
*/
static void
setup_recon_GA(RECON_T *recon, const PARMS_T *parms, const POWFS_T *powfs){
    loc_t *ploc=recon->ploc;
    const int nwfs=parms->nwfsr;
    const int ndm=parms->ndm;
    if(parms->load.GA){
	warning2("Loading saved GA\n");
	recon->GA=dspcellread("%s", parms->load.GA);
	if(recon->GA->nx!=nwfs || recon->GA->ny!=ndm)
	    error("Wrong saved GA\n");
	for(int idm=0; idm<ndm; idm++){
	    int nloc=recon->aloc->p[idm]->nloc;
	    for(int iwfs=0; iwfs<nwfs; iwfs++){
		int ipowfs = parms->wfsr[iwfs].powfs;
		if(parms->sim.skysim && parms->powfs[ipowfs].lo){
		    continue;
		}
		int nsa=powfs[ipowfs].saloc->nloc;
		if(IND(recon->GA,iwfs,idm)->nx!=nsa*2 || (!parms->recon.modal && IND(recon->GA,iwfs,idm)->ny!=nloc)
		   || (parms->recon.modal && IND(recon->GA,iwfs,idm)->ny!=recon->amod->p[idm]->ny)){
		    error("Wrong saved GA\n");
		}
	    }
	}
    }else{
	info2("Generating GA ");TIC;tic;
	recon->GA=dspcellnew(nwfs, ndm);
	if(parms->recon.modal) {
	    recon->GM=dcellnew(nwfs, ndm);
	    
	    recon->amod=dcellnew(ndm, 1);
	    recon->anmod=lnew(ndm, 1);
	    for(int idm=0; idm<ndm; idm++){
		int nmod=parms->recon.nmod;
		const long nloc=recon->aloc->p[idm]->nloc;
		switch(parms->recon.modal){
		case -2: {//dummy modal control, emulating zonal mode with identity modal matrix
		    if(nmod && nmod!=nloc){
			warning("recon.mod should be 0 or %ld when recon.modal=2 \n",nloc);
		    }
		    recon->amod->p[idm]=dnew(nloc, nloc);
		    double val=sqrt(nloc);
		    daddI(recon->amod->p[idm], val);
		    dadds(recon->amod->p[idm], -val/nloc);
		}
		    break;
		case -1://zernike
		    {
			if(!nmod) nmod=nloc;
			int rmax=floor((sqrt(1+8*nmod)-3)*0.5);
			recon->amod->p[idm]=zernike(recon->aloc->p[idm], 0, 0, rmax, 0);
		    }
		    break;
		case 1://Karhunen loeve. Don't limit number of modes here to make caching of G_M right.
		    recon->amod->p[idm]=KL_vonkarman(recon->aloc->p[idm], 0, parms->atmr.L0);
		    break;
		default:
		    error("Invalid recon.modal");
		}	    
		recon->anmod->p[idm]=recon->amod->p[idm]->ny;
	    }
	}
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs = parms->wfsr[iwfs].powfs;
	    if(parms->sim.skysim && parms->powfs[ipowfs].lo){
		continue;
	    }
	    if(parms->powfs[ipowfs].skip==2){//no need for TWFS
		continue;
	    }
	    double  hs = parms->wfs[iwfs].hs;
	    for(int idm=0; idm<ndm; idm++){
		double  ht = parms->dm[idm].ht;
		double  scale=1. - ht/hs;
		double  dispx=0, dispy=0;
		if(!parms->recon.glao){
		    dispx=parms->wfsr[iwfs].thetax*ht;
		    dispy=parms->wfsr[iwfs].thetay*ht;
		}
		if(parms->powfs[ipowfs].type==1){//PWFS
		    if(!parms->powfs[ipowfs].lo){
			if(!parms->recon.modal){
			    info2("\nPyWFS from aloc to saloc directly\n");
			    dmat *tmp=pywfs_mkg(powfs[ipowfs].pywfs, recon->aloc->p[idm], 0, dispx, dispy, scale);
			    IND(recon->GA, iwfs, idm)=d2sp(tmp, dmaxabs(tmp)*1e-6);
			    dfree(tmp);
			}else{
			    info2("\nPyWFS from amod to saloc directly\n");
			    //We compute the GM for full set of modes so that it is cached only once.
			    IND(recon->GM, iwfs, idm)=pywfs_mkg(powfs[ipowfs].pywfs, recon->aloc->p[idm], recon->amod->p[idm], dispx, dispy, scale);
			}
		    }
		}else{//SHWFS
		    int freeloc=0;
		    loc_t *loc=ploc;
		    if(parms->recon.misreg_dm2wfs && parms->recon.misreg_dm2wfs[iwfs+idm*nwfs]){
			loc=loctransform(loc, parms->recon.misreg_dm2wfs[iwfs+idm*nwfs]);
			freeloc=1;
		    }
		    dsp *H=mkh(recon->aloc->p[idm], loc, dispx,dispy,scale);
		    IND(recon->GA, iwfs, idm)=dspmulsp(recon->GP->p[ipowfs], H,"nn");
		    dspfree(H);
		    if(freeloc){
			locfree(loc);
		    }
		    if(parms->recon.modal){
			dspmm(PIND(recon->GM, iwfs, idm), IND(recon->GA, iwfs, idm), recon->amod->p[idm], "nn",1);
		    }
		}
	    }/*idm */
	}
    	toc2(" ");
    }
    if(parms->recon.modal && parms->recon.nmod>0){
	for(int idm=0; idm<ndm; idm++){
	    if(parms->recon.nmod<recon->amod->p[idm]->ny){
		warning("Reduce number of controlled modes from %ld to %d\n",
			recon->amod->p[idm]->ny, parms->recon.nmod);
		dresize(recon->amod->p[idm], 0, parms->recon.nmod);
		recon->anmod->p[idm]=recon->amod->p[idm]->ny;
		for(int iwfs=0; iwfs<nwfs; iwfs++){
		    if(!IND(recon->GM, iwfs, idm)) continue;
		    dresize(IND(recon->GM, iwfs, idm), 0,  parms->recon.nmod);
		}
	    }
	}
    }
    if(parms->save.setup){
	writebin(recon->GA, "GA");
	if(parms->recon.modal){
	    writebin(recon->amod, "amod");
	    writebin(recon->GM, "GM");
	}
    }
    if(parms->recon.alg==1 && !parms->recon.modal){
	recon->actcpl=genactcpl(recon->GA, 0);
	act_stuck(recon->aloc, recon->actcpl, recon->actfloat);
	
	if(recon->actstuck){
	    /*This is need for LSR reconstructor to skip stuck actuators.  GA is
	      also used to form PSOL gradients, but that one doesn't need this
	      modification because actuator extropolation was already applied.*/
	    warning2("Apply stuck actuators to GA\n");
	    act_stuck(recon->aloc, recon->GA,recon->actstuck);
	    if(parms->save.setup){
		writebin(recon->GA,"GA_stuck");
	    }
	}
	if(parms->lsr.actinterp){
	    recon->actinterp=act_extrap(recon->aloc, recon->actcpl, parms->lsr.actthres);
	}else if(recon->actfloat){
	    warning("There are float actuators, but fit.actinterp is off\n");
	}
	if(recon->actinterp){
	    dspcell *GA2=0;
	    dcellmm(&GA2, recon->GA, recon->actinterp, "nn", 1);
	    dspcellfree(recon->GA);
	    recon->GA=GA2;
	}
	if(parms->save.setup){
	    if(recon->actinterp){
		writebin(recon->actinterp, "actinterp");
	    }
	    if(recon->actcpl){
		writebin(recon->actcpl, "actcpl");
	    }
	}
    }
    /*Create GAlo that only contains GA for low order wfs */
    recon->GAlo=dspcellnew(recon->GA->nx, recon->GA->ny);
    recon->GAhi=dspcellnew(recon->GA->nx, recon->GA->ny);
    if(parms->recon.modal) recon->GMhi=dcellnew(nwfs, ndm);

    for(int idm=0; idm<ndm; idm++){
	for(int iwfs=0; iwfs<nwfs; iwfs++){
	    int ipowfs=parms->wfsr[iwfs].powfs;
	    if(parms->powfs[ipowfs].lo
	       || (parms->recon.split && parms->nlopowfs==0 && !parms->powfs[ipowfs].trs)){/*for low order wfs */
		IND(recon->GAlo, iwfs, idm)=dspref(IND(recon->GA, iwfs, idm));
	    }
	    if(!parms->powfs[ipowfs].skip){
		IND(recon->GAhi, iwfs, idm)=dspref(IND(recon->GA, iwfs, idm));
		if(parms->recon.modal){
		    IND(recon->GMhi, iwfs, idm)=dref(IND(recon->GM, iwfs, idm));
		}
	    }
	}
    }
}
/**
   Crate the xloc to wfs gradient operator.
*/
static void 
setup_recon_GX(RECON_T *recon, const PARMS_T *parms){
    const int nwfs=parms->nwfsr;
    const int npsr=recon->npsr;
    recon->GX=dspcellnew(nwfs, npsr);
    dspcell* GX=recon->GX/*PDSPCELL*/;
    dspcell* HXW=recon->HXW/*PDSPCELL*/;
    info2("Generating GX ");TIC;tic;
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	/*gradient from xloc. Also useful for lo WFS in MVST mode. */
	for(int ips=0; ips<npsr; ips++){
	    IND(GX,iwfs,ips)=dspmulsp(recon->GP2->p[iwfs], IND(HXW,iwfs,ips),"nn");
	}/*ips */
    }
    toc2(" ");
    recon->GXtomo=dspcellnew(recon->GX->nx, recon->GX->ny);
    dspcell* GXtomo=recon->GXtomo/*PDSPCELL*/;

    recon->GXlo=dspcellnew(recon->GX->nx, recon->GX->ny);
    dspcell*  GXlo=recon->GXlo/*PDSPCELL*/;

    int nlo=parms->nlopowfs;
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs=parms->wfsr[iwfs].powfs;
	for(int ips=0; ips<npsr; ips++){
	    if(!parms->powfs[ipowfs].skip){/*for tomography */
		IND(GXtomo,iwfs,ips)=dspref(IND(GX,iwfs,ips));
	    }
	    if(parms->powfs[ipowfs].lo 
	       || (parms->recon.split && nlo==0 && !parms->powfs[ipowfs].trs)){
		/*for low order wfs or extracted t/t for high order ngs wfs.*/
		IND(GXlo,iwfs,ips)=dspref(IND(GX,iwfs,ips));
	    }
	}
    }/*iwfs */
}

/**
   From focus mode to gradients
 */
static void
setup_recon_GF(RECON_T *recon, const PARMS_T *parms){
    /*Create GFall: Focus mode -> WFS grad. This is model*/
    recon->GFall=dcellnew(parms->npowfs, 1);
    recon->GFngs=dcellnew(parms->nwfs, 1);
    {
	dmat *opd=dnew(recon->ploc->nloc,1);
	loc_add_focus(opd->p, recon->ploc, 1);
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){ 
	    dspmm(&recon->GFall->p[ipowfs], recon->GP->p[ipowfs], opd, "nn", 1);
	    if(parms->powfs[ipowfs].lo && !parms->powfs[ipowfs].llt){
		for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
		    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
		    recon->GFngs->p[iwfs]=dref(recon->GFall->p[ipowfs]);
		}
	    }
	}
	dfree(opd);
    }
    if(parms->save.setup){
	writebin(recon->GFall,"GFall");
    }
}
/**
   From radial order modes to gradients.
 */
static void
setup_recon_GR(RECON_T *recon, const POWFS_T *powfs, const PARMS_T *parms){
    recon->GRall=dcellnew(parms->npowfs, 1);
    dmat *opd=zernike(recon->ploc, parms->aper.d, 3, parms->powfs[parms->itpowfs].order, 1);
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].skip==2 || parms->powfs[ipowfs].llt){
	    if(parms->powfs[ipowfs].type==1){//PWFS
		recon->GRall->p[ipowfs]=pywfs_mkg(powfs[ipowfs].pywfs, recon->ploc, opd, 0, 0, 1);
	    }else{//SHWFS
		dspmm(&recon->GRall->p[ipowfs], recon->GP->p[ipowfs], opd, "nn", 1);
	    }
	}
    }
    dfree(opd);
}
void setup_recon_dmttr(RECON_T *recon, const PARMS_T *parms){
    recon->DMTT=dcellnew(parms->ndm, 1);
    recon->DMPTT=dcellnew(parms->ndm, 1);
    if(!recon->actcpl){
	error("actcpl must not be null\n");
    }
    for(int idm=0; idm<parms->ndm; idm++){
	recon->DMTT->p[idm]=loc2mat(recon->aloc->p[idm], 0);
    }
    act_zero(recon->aloc, recon->DMTT, recon->actstuck);
    for(int idm=0; idm<parms->ndm; idm++){
	recon->DMPTT->p[idm]=dpinv(recon->DMTT->p[idm], 0);
    }
    if(parms->save.setup){
	writebin(recon->DMTT, "DMTT");
	writebin(recon->DMPTT, "DMPTT");
    }
}

/**
   Setup fitting low rank terms that are in the NULL space of DM fitting
   operator. typically include piston on each DM and tip/tilt on certain
   DMs. Becareful with tip/tilt contraint when using CBS.  */
static void 
fit_prep_lrt(RECON_T *recon, const PARMS_T *parms){
    const int ndm=parms->ndm;
    if(ndm>=3) warning("Low rank terms for 3 or more dms are not tested\n");
    recon->fitNW=dcellnew(ndm,1);
    double scl=recon->fitscl=1./recon->floc->nloc;
    if(fabs(scl)<1.e-15){
	error("recon->fitscl is too small\n");
    }
    /*computed outside. */
    int lrt_tt=parms->fit.lrt_tt;
    int nnw=0;
    if(parms->fit.lrt_piston){
	nnw+=ndm;
    }
    if(lrt_tt){
	nnw+=2*(ndm-1);
    }
    if(nnw==0) return;
    dcell* actcpl=dcelldup(recon->actcpl);
    //include stuck actuator
    act_stuck(recon->aloc, actcpl, recon->actstuck);
    for(int idm=0; idm<ndm; idm++){
	int nloc=recon->aloc->p[idm]->nloc;
	recon->fitNW->p[idm]=dnew(nloc, nnw);
    }
    int inw=0;/*current column */
    if(parms->fit.lrt_piston){
	info2("Adding piston cr to fit matrix\n");
	for(int idm=0; idm<ndm; idm++){
	    int nloc=recon->aloc->p[idm]->nloc;
	    double *p=recon->fitNW->p[idm]->p+(inw+idm)*nloc;
	    const double *cpl=actcpl->p[idm]->p;
	    for(int iloc=0; iloc<nloc; iloc++){
		if(cpl[iloc]>0.1){
		    //don't count floating or stuck actuators
		    p[iloc]=scl;
		}
	    }
	}
	inw+=ndm;
    }
    if(lrt_tt){
	double factor=0;
	if(lrt_tt==1){
	    info2("Adding TT cr on upper DMs to fit matrix.\n");
	    factor=scl*2./parms->aper.d;
	    for(int idm=1; idm<ndm; idm++){
		int nloc=recon->aloc->p[idm]->nloc;
		double *p=recon->fitNW->p[idm]->p+(inw+(idm-1)*2)*nloc;
		double *p2x=p;
		double *p2y=p+nloc;
		const double *cpl=actcpl->p[idm]->p;
		for(int iloc=0; iloc<nloc; iloc++){
		    if(cpl[iloc]>0.1){
			p2x[iloc]=recon->aloc->p[idm]->locx[iloc]*factor;/*x tilt */
			p2y[iloc]=recon->aloc->p[idm]->locy[iloc]*factor;/*y tilt */
		    }
		}
	    }
	}else if(lrt_tt==2){/*Canceling TT. only valid for 2 DMs */
	    warning("Adding Canceling TT cr to fit matrix. Deprecated\n");
	    if(ndm!=2){
		error("Only ndm=2 case is implemented\n");
	    }
	    for(int idm=0; idm<ndm; idm++){
		int nloc=recon->aloc->p[idm]->nloc;
		double *p=recon->fitNW->p[idm]->p+inw*nloc;
		if(idm==0) factor=scl*2/parms->aper.d;
		else if(idm==1) factor=-scl*2./parms->aper.d;
		double *p2x=p;
		double *p2y=p+nloc;
		const double *cpl=actcpl->p[idm]->p;
		for(int iloc=0; iloc<nloc; iloc++){
		    if(cpl[iloc]>0.1){
			p2x[iloc]=recon->aloc->p[idm]->locx[iloc]*factor;/*x tilt */
			p2y[iloc]=recon->aloc->p[idm]->locy[iloc]*factor;/*y tilt */
		    }
		}
	    }

	}
	inw+=2*(ndm-1);
    }
    if(parms->fit.actslave){
	/*
	  2011-07-19: When doing PSFR study for MVR with SCAO, NGS. Found
	  that slaving is causing mis-measurement of a few edge
	  actuators. First try to remove W1. Or lower the weight. Revert
	  back.
	  1./floc->nloc is on the same order of norm of Ha'*W*Ha. 
	*/
	TIC;tic;
	recon->actslave=slaving(recon->aloc, recon->actcpl,
				recon->fitNW, recon->actstuck,
				recon->actfloat, parms->fit.actthres, 1./recon->floc->nloc);
	toc2("slaving");
	if(parms->save.setup){
	    writebin(recon->actslave,"actslave");
	}
    }
    if(parms->save.setup){
	writebin(recon->fitNW,"fitNW");
    }
    cellfree(actcpl);
}
/**
   setting up global tip/tilt remove operator from LGS gradients.
*/

static void
setup_recon_TT(RECON_T *recon, const PARMS_T *parms, const POWFS_T *powfs){
    int nwfs=parms->nwfsr;
    recon->TT=dcellnew(nwfs,nwfs);
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].nwfs==0) continue;
	if(parms->powfs[ipowfs].trs 
	   || (parms->recon.split && !parms->powfs[ipowfs].lo)
	   || parms->powfs[ipowfs].dither==1){
	    int nsa=powfs[ipowfs].saloc->nloc;
	    dmat *TT=0;
	    if(parms->powfs[ipowfs].type==0){//SHWFS
		TT=dnew(nsa*2,2);
		double *TTx=TT->p;
		double *TTy=TT->p+nsa*2;
		for(int isa=0; isa<nsa; isa++){
		    TTx[isa]=1;
		    TTy[isa]=0;
		}
		for(int isa=nsa; isa<nsa*2; isa++){
		    TTx[isa]=0;
		    TTy[isa]=1;
		}
	    }else if(parms->powfs[ipowfs].type==1){//PYWFS
		TT=pywfs_tt(powfs[ipowfs].pywfs);
	    }else{
		error("Invalid powfs.type\n");
	    }
	    if(parms->recon.glao){
		recon->TT->p[ipowfs*(parms->npowfs+1)]=ddup(TT);
	    }else{
		for(int jwfs=0; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
		    int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
		    recon->TT->p[iwfs*(1+parms->nwfsr)]=ddup(TT);
		}
	    }
	    dfree(TT);
	}
    }
    if(parms->save.setup){
	writebin(recon->TT, "TT");
    }
}
/**
   operator to remove diffrential focus modes that might be caused by sodium layer
   horizontal structure.
*/

static void
setup_recon_DF(RECON_T *recon, const PARMS_T *parms, const POWFS_T *powfs){
    if(!recon->has_dfr) return;

    int nwfs=parms->nwfsr;
    recon->DF=dcellnew(nwfs,nwfs);
    /*Then differential focus modes. */
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].nwfs==0) continue;
	if(parms->powfs[ipowfs].dfrs){
	    if(parms->powfs[ipowfs].nwfs<2){
		error("This powfs group has only 1 wfs. Could not remove diff focus\n");
	    }
	    int nsa=powfs[ipowfs].saloc->nloc;
	    dmat* DF=dnew(nsa*2,1);
	    /*saloc is lower left corner of subaperture. don't have to be the center. */
	    memcpy(DF->p, powfs[ipowfs].saloc->locx, sizeof(double)*nsa);
	    memcpy(DF->p+nsa, powfs[ipowfs].saloc->locy, sizeof(double)*nsa);
	    /**
	       postive focus on first wfs. negative focus on diagnonal wfs.
	    */
	    for(int jwfs=1; jwfs<parms->powfs[ipowfs].nwfs; jwfs++){
		int iwfs=parms->powfs[ipowfs].wfs->p[jwfs];
		if(parms->powfs[ipowfs].skip){
		    error("This WFS %d should be included in Tomo.\n", iwfs);
		}
		dcp(&recon->DF->p[iwfs*nwfs], DF);
		dadd(&recon->DF->p[iwfs+iwfs*nwfs], 0, DF, -1);
	    }
	    dfree(DF);
	}
    }
    if(parms->save.setup){
	writebin(recon->DF, "DF");
    }
}
/**
   Create reconstruction parameters that are related to the geometry only, and
   will not be updated when estimated WFS measurement noise changes.
   
   This can be used to do NCPA calibration.
 */
RECON_T *setup_recon_prep(const PARMS_T *parms, const APER_T *aper, const POWFS_T *powfs){
    RECON_T * recon = mycalloc(1,RECON_T); 
    if(parms->recon.warm_restart){
	info2("Using warm restart\n");
    }else{
	warning2("Do not use warm restart\n");
    }
    if(parms->cn2.pair && parms->cn2.pair->nx>0 && !recon->cn2est){
	/*setup CN2 Estimator. It determines the reconstructed layer heigh can be fed to the tomography */
	recon->cn2est=cn2est_prepare(parms, powfs);
    }
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].nwfs==0) continue;
	/*we will remove tip/tilt from the high order NGS wfs in split tomo mode.*/
	if(parms->powfs[ipowfs].trs || (parms->recon.split && !parms->powfs[ipowfs].lo)){
	    recon->has_ttr=1;
	    break;
	}
    }
    if(!parms->recon.glao){
	for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	    if(parms->powfs[ipowfs].nwfs<=1) continue;
	    if(parms->powfs[ipowfs].dfrs){
		recon->has_dfr=1;
		break;
	    }
	}   
    }
    recon->ngrad=lnew(parms->nwfsr, 1);
    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	const int ipowfs=parms->wfsr[iwfs].powfs;
	if(!parms->powfs[ipowfs].skip){
	    recon->ngrad->p[iwfs]=powfs[ipowfs].saloc->nloc*2;
	}
    }
    /*to be used in tomography. */
    recon->nthread=NTHREAD;
    /*setup pupil coarse grid for gradient operator*/
    setup_recon_ploc(recon,parms);
    /*fine sampled pupil loc*/
    setup_recon_gloc(recon,parms, aper);
    /*setup DM actuator grid */
    setup_recon_aloc(recon,parms);
    /*Grid for DM fitting*/
    setup_recon_floc(recon,parms);
    /*Gradient operators*/
    setup_recon_GP(recon, parms, powfs);
    //TT Removal
    setup_recon_TT(recon,parms,powfs);
    //DF removal. //Deprecated
    if(recon->has_dfr){
	setup_recon_DF(recon,parms,powfs);
    }
    if(recon->DF){
	recon->TTF=dcellcat(recon->TT,recon->DF,2);
    }else{
	recon->TTF=dcellref(recon->TT);
    }
    if(parms->recon.alg==0 && !parms->sim.idealfit){//tomography parameters
	if(parms->cn2.tomo){
	    /*Use cn2 estimation results for tomography. Use its ht to build
	      reconstructor matrices.*/
	    cn2est_t *cn2est=recon->cn2est;
	    recon->ht=dref(cn2est->htrecon);
	    recon->os=dref(cn2est->os);
	    recon->wt=dref(cn2est->wtrecon->p[0]);
	    /*the following will be updated later in simulation. */
	    if(parms->cn2.keepht){
		for(int ips=0; ips<recon->wt->nx; ips++){
		    recon->wt->p[ips]=parms->atmr.wt->p[ips];
		}
	    }else{
		dset(recon->wt, 1./recon->wt->nx);/*evenly distributed.  */
	    }
	}else{/*use input information from atmr */
	    recon->wt=dnew(parms->atmr.nps,1);
	    recon->ht=dnew(parms->atmr.nps,1);
	    recon->os=dnew(parms->atmr.nps,1);
	    for(int ips=0; ips<recon->wt->nx; ips++){
		recon->wt->p[ips]=parms->atmr.wt->p[ips];
		recon->ht->p[ips]=parms->atmr.ht->p[ips];
		recon->os->p[ips]=parms->atmr.os->p[ips];
	    }
	}
	recon->r0=parms->atmr.r0;
	recon->L0=parms->atmr.L0;

	/*sampling of xloc */
	recon->dx=dnew(recon->ht->nx, 1);
	for(int iht=0; iht<recon->ht->nx; iht++){
	    double scale = 1.0 - recon->ht->p[iht]/parms->atmr.hs;
	    recon->dx->p[iht]=(parms->atmr.dx/recon->os->p[iht])*scale;	
	}
	/*number of reconstruction layers */
	recon->npsr= recon->ht->nx;
	/*setup atm reconstruction layer grid */
	setup_recon_xloc(recon, parms);
	/*setup xloc/aloc to WFS grad */
	setup_recon_HXW(recon,parms);
	setup_recon_GX(recon,parms);//uses HXW.
	dspcellfree(recon->HXW);/*only keep HXWtomo for tomography */

    }
    if(parms->recon.alg==0 || parms->sim.ncpa_calib){
	setup_recon_HA(recon,parms);
	fit_prep_lrt(recon, parms);
    }
    setup_recon_dmttr(recon, parms);
    return recon;
}
/**
   That may depend on GPU data.
 */
void setup_recon_prep2(RECON_T *recon, const PARMS_T *parms, const APER_T *aper, const POWFS_T *powfs){
    setup_recon_GA(recon, parms, powfs);
    setup_recon_GF(recon, parms);
    if(parms->itpowfs!=-1){
	setup_recon_GR(recon,powfs,parms);
    }
    if(parms->recon.split){
	setup_ngsmod_prep(parms,recon,aper,powfs);
    }
}
