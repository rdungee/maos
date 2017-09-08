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

/*
  The 5 NGS mode for split tomography with 2DM
  I work in not normalized zernike space. So the result is in radians,
  multiply to 2R to get zernike mode.
  2010-07-23: Added tikholnov regularization to Wa. 
*/
#include "common.h"
#include "ahst.h"
#include "pywfs.h"
/**
   \file maos/ahst.c Contains functions to setup NGS modes and reconstructor
   using AHST for one or more DMs.  Use parms->wfsr instead of parms->wfs for wfs
   information, which hands GLAO mode correctly.

   Notice that update of this file may require GPU code update accordingly
*/

static TIC;

/**
   computes the cross-coupling of NGS modes in science field.
   MCC=(M'*Ha'*W*Ha*M) where M is ngsmod on DM, Ha is propagator from DM to
   science. W is weighting in science.
*/
static dcell* ngsmod_mcc(const PARMS_T *parms, RECON_T *recon, const APER_T *aper, const double *wt){
    NGSMOD_T *ngsmod=recon->ngsmod;
    const loc_t *plocs=aper->locs;
    double *x=plocs->locx;
    double *y=plocs->locy;
    int nloc=plocs->nloc;
    double *amp=aper->amp->p;
    const int nmod=ngsmod->nmod;
    dcell *mcc=dcellnew(parms->evl.nevl,1);
    dmat *aMCC=aper->mcc;
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	dmat *MCC=mcc->p[ievl]=dnew(nmod,nmod);
	IND(MCC,0,0)=IND(aMCC,1,1);
	IND(MCC,1,1)=IND(aMCC,2,2);
	IND(MCC,0,1)=IND(MCC,1,0)=IND(aMCC,2,1);
    }
    
    if(ngsmod->nmod>2){
	tic;
	double *mod[nmod];
	mod[0]=x;
	mod[1]=y;
	for(int imod=2; imod<nmod; imod++){
	    mod[imod]=mymalloc(nloc,double);
	}
	/*dc component of the focus mod. subtract during evaluation. */
	/*this is not precisely R^2/2 due to obscuration */

	const double MCC_fcp=aper->fcp;
	const double ht=ngsmod->ht;
	const double scale=ngsmod->scale;
	const double scale1=1.-scale;
	//dmat *modvec=dnew(nmod, 1);
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    dmat *MCC=mcc->p[ievl];
	    if(fabs(wt[ievl])<1.e-12) continue;
	    double thetax=parms->evl.thetax->p[ievl];
	    double thetay=parms->evl.thetay->p[ievl];
	    /*for(int imod=2; imod<nmod; imod++){
		dmat *iopd=dnew_ref(nloc, 1, mod[imod]);
		dzero(iopd);
		dzero(modvec); 
		modvec->p[imod]=1;
		ngsmod2science(iopd, plocs, ngsmod, thetax, thetay, modvec->p, 1);
		dfree(iopd);
		}*/
	    
	    for(int iloc=0; iloc<nloc; iloc++){
		double xx=x[iloc]*x[iloc];
		double xy=x[iloc]*y[iloc];
		double yy=y[iloc]*y[iloc];
		//remove piston in focus 
		if(ngsmod->indps){
		    if(ngsmod->ahstfocus){
			mod[ngsmod->indps][iloc]=//mod[2] has no focus effect on science
			    -2.*ht*scale*(thetax*x[iloc]+thetay*y[iloc]);
		    }else{
			mod[ngsmod->indps][iloc]=scale1*(xx+yy-MCC_fcp)
			    -2.*ht*scale*(thetax*x[iloc]+thetay*y[iloc]);
		    }
		    mod[ngsmod->indps+1][iloc]=scale1*(xx-yy)
			-2.*ht*scale*(thetax*x[iloc]-thetay*y[iloc]);
		    mod[ngsmod->indps+2][iloc]=scale1*(xy)
			-ht*scale*(thetay*x[iloc]+thetax*y[iloc]);
		}
		if(ngsmod->indastig){
		    mod[ngsmod->indastig][iloc]=(xx-yy);
		    mod[ngsmod->indastig+1][iloc]=(xy);
		}
		if(ngsmod->indfocus){
		    mod[ngsmod->indfocus][iloc]=(xx+yy-MCC_fcp);//for focus tracking
		}
	    }

	    for(int jmod=0; jmod<nmod; jmod++){
		/*if(parms->save.setup){
		    writedbl(mod[jmod], nloc, 1, "ahst_mode_%d", jmod);
		    }*/
		for(int imod=jmod; imod<nmod; imod++){
		    if(imod<2&&jmod<2) continue;
		    double tmp=dotdbl(mod[imod],mod[jmod],amp,nloc);
		    IND(MCC,jmod,imod)=tmp;
		    if(imod!=jmod){
			IND(MCC,imod,jmod)=IND(MCC,jmod,imod);
		    }
		}
	    }
	}//for(ievl)
	//dfree(modvec);
	for(int imod=2; imod<nmod; imod++){
	    free(mod[imod]);
	}
	toc2("mcc");
    }
 
    return mcc;
}
/**
   Compute NGS mode aperture weighting using science field.  Wa=Ha'*W*Ha where
   Ha is from ALOC to PLOCS and W is the amplitude weighting in PLOCS when
   use_ploc==0. Otherwise, Ha is from ALOC to PLOC and W is the amplitude
   weighting in PLOC.  */
static dspcell *ngsmod_Wa(const PARMS_T *parms, RECON_T *recon, 
			 const APER_T *aper, int use_ploc){
    const double *wt=parms->evl.wt->p;
    const int ndm=parms->ndm;
    loc_t *loc;
    double *amp=NULL;
    if(use_ploc){
	loc=recon->floc;
	amp=mycalloc(loc->nloc,double);
	prop_nongrid_bin(aper->locs,aper->amp->p,loc,amp,1,0,0,1);
	normalize_sum(amp,loc->nloc,1);
    }else{
	amp=aper->amp->p;
	loc=aper->locs;
    }
    dspcell *Wa=NULL;
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	if(fabs(wt[ievl])<1.e-12) continue;
	double thetax=parms->evl.thetax->p[ievl];
	double thetay=parms->evl.thetay->p[ievl];
	dspcell *Hat=dspcellnew(ndm,1);
	for(int idm=0; idm<ndm; idm++){
	    double hc = parms->dm[idm].ht;
	    double displacex=thetax*hc;
	    double displacey=thetay*hc;
	    /*from DM to ploc (plocs) science beam */
	    Hat->p[idm]=mkhb(recon->aloc->p[idm], loc, displacex,displacey,1.);
	    dspmuldiag(Hat->p[idm], amp, wt[ievl]);
	}
	dcellmm(&Wa, Hat, Hat, "nt", 1);
	dspcellfree(Hat);
    }
    if(use_ploc){
	free(amp);
    }
    return Wa;
}
/**
   compute NGS mode removal Pngs from LGS commands using aperture
   weighting. Pngs=(MCC)^-1 (Hm'*W*Ha).

   2012-05-25: The NGS mode removal should be based on old five modes even if now focus on PS1 is merged with defocus mode
*/
static dcell* ngsmod_Pngs_Wa(const PARMS_T *parms, RECON_T *recon, 
		     const APER_T *aper, int use_ploc){

    NGSMOD_T *ngsmod=recon->ngsmod;
    const double ht=ngsmod->ht;
    const double scale=ngsmod->scale;
    const double scale1=1.-scale;
    const double *wt=parms->evl.wt->p;
    const int ndm=parms->ndm;
    const int nmod=ngsmod->nmod;
    loc_t *loc;
    double *x, *y;
    int nloc;
    double *amp=NULL;
    if(use_ploc){
	loc=recon->floc;
	amp=mycalloc(loc->nloc,double);
	prop_nongrid_bin(aper->locs,aper->amp->p,loc,amp,1,0,0,1);
	normalize_sum(amp,loc->nloc,1);
    }else{
	amp=aper->amp->p;
	loc=aper->locs;
    }
    x=loc->locx;
    y=loc->locy;
    nloc=loc->nloc;

    dcell *modc=dcellnew(1,1);/*W*Hm*M */
    modc->p[0]=dnew(nloc,nmod);
    dmat *mod=modc->p[0];
    for(int iloc=0; iloc<nloc; iloc++){
	IND(mod,iloc,0)=x[iloc]*amp[iloc];
	IND(mod,iloc,1)=y[iloc]*amp[iloc];
    }
    const double MCC_fcp=ngsmod->aper_fcp;
    /*dc component of the focus mod. subtract during evaluation. */
    /*this is not precisely R^2/2 due to obscuration */
    dcell *HatWHmt=dcellnew(ndm,1);
    dsp *HatGround=NULL;
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	if(fabs(wt[ievl])<1.e-12) continue;
	double thetax=parms->evl.thetax->p[ievl];
	double thetay=parms->evl.thetay->p[ievl];
	if(nmod>2){
	    for(int iloc=0; iloc<nloc; iloc++){
		double xx=x[iloc]*x[iloc];
		double xy=x[iloc]*y[iloc];
		double yy=y[iloc]*y[iloc];
		/*remove piston in focus */
		if(ngsmod->indps){
		    if(ngsmod->ahstfocus){
			IND(mod,iloc,ngsmod->indps)=amp[iloc]
			    *(-2.*ht*scale*(thetax*x[iloc]+thetay*y[iloc]));
		    }else{
			IND(mod,iloc,ngsmod->indps)=amp[iloc]
			    *(scale1*(xx+yy-MCC_fcp)
			      -2.*ht*scale*(thetax*x[iloc]+thetay*y[iloc]));
		    }
		    IND(mod,iloc,ngsmod->indps+1)=amp[iloc]
			*(scale1*(xx-yy)
			  -2.*ht*scale*(thetax*x[iloc]-thetay*y[iloc]));
		    IND(mod,iloc,ngsmod->indps+2)=amp[iloc]
			*(scale1*(xy)
			  -ht*scale*(thetay*x[iloc]+thetax*y[iloc]));
		}
		if(ngsmod->indastig){
		    IND(mod,iloc,ngsmod->indastig)=amp[iloc]*(xx-yy);
		    IND(mod,iloc,ngsmod->indastig+1)=amp[iloc]*(xy);
		}
		if(ngsmod->indfocus){
		    IND(mod,iloc,ngsmod->indfocus)=amp[iloc]*(xx+yy-MCC_fcp);
		}
	    }
	}
	dspcell *Hat=dspcellnew(ndm,1);
	for(int idm=0; idm<ndm; idm++){
	    double hc = parms->dm[idm].ht;
	    double displacex=thetax*hc;
	    double displacey=thetay*hc;
	    if(parms->dm[idm].isground && HatGround){
		info("Reusing HatGround\n");
		Hat->p[idm]=dspref(HatGround);
	    }else{
		/*from DM to ploc (plocs) science beam */
		Hat->p[idm]=mkhb(recon->aloc->p[idm], loc, displacex,displacey,1.);
		if(parms->dm[idm].isground){
		    HatGround=dspref(Hat->p[idm]);
		}
	    }
	}
	dcellmm(&HatWHmt,Hat,modc,"nn",wt[ievl]);
	dspcellfree(Hat);
    }
    dspfree(HatGround);
    dcell *IMCC=dcellnew(1,1);
    IMCC->p[0]=dref(ngsmod->IMCC);
    dcell *Pngs=NULL;
    dcellmm(&Pngs,IMCC,HatWHmt,"nt",1);
    dcellfree(IMCC);
    dcellfree(modc);
    if(use_ploc){
	free(amp);
    }
    dcellfree(HatWHmt);
    return Pngs;
}
/**
   compute tip/tilt mode removal from each DM commands using aperture
   weighting. Ptt=(MCC_TT)^-1 *(Hmtt * W * Ha)
*/
static dcell* ngsmod_Ptt_Wa(const PARMS_T *parms, RECON_T *recon, 
			    const APER_T *aper, int use_ploc){
    NGSMOD_T *ngsmod=recon->ngsmod;
    const double *wt=parms->evl.wt->p;
    const int ndm=parms->ndm;
    const int nmod=2;
    loc_t *loc;
    double *x, *y;
    int nloc;
    double *amp=NULL;
    if(use_ploc){
	loc=recon->floc;
	amp=mycalloc(loc->nloc,double);
	prop_nongrid_bin(aper->locs, aper->amp->p,loc,amp,1,0,0,1);
	normalize_sum(amp,loc->nloc,1);
    }else{
	amp=aper->amp->p;
	loc=aper->locs;
    }
    x=loc->locx;
    y=loc->locy;
    nloc=loc->nloc;

    dcell *modc=dcellnew(1,1);/*W*Hm*M */
    modc->p[0]=dnew(nloc,nmod);
    dmat *mod=modc->p[0];
    for(int iloc=0; iloc<nloc; iloc++){
	IND(mod,iloc,0)=x[iloc]*amp[iloc];
	IND(mod,iloc,1)=y[iloc]*amp[iloc];
    }
    dcell *HatWHmt=dcellnew(ndm,1);
    dsp *HatGround=NULL;
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	if(fabs(wt[ievl])<1.e-12) continue;
	double thetax=parms->evl.thetax->p[ievl];
	double thetay=parms->evl.thetay->p[ievl];
	dspcell *Hat=dspcellnew(ndm,1);
	for(int idm=0; idm<ndm; idm++){
	    double hc = parms->dm[idm].ht;
	    double displacex=thetax*hc;
	    double displacey=thetay*hc;
	    if(!parms->dm[idm].isground || !HatGround){
		/*from DM to ploc (plocs) science beam */
		Hat->p[idm]=mkhb(recon->aloc->p[idm], loc, displacex,displacey,1.);
		if(parms->dm[idm].isground){
		    HatGround=dspref(Hat->p[idm]);
		}
	    }else{
		Hat->p[idm]=dspref(HatGround);
	    }
	    dspmm(&HatWHmt->p[idm],Hat->p[idm],modc->p[0],"nn",wt[ievl]);
	}
	dspcellfree(Hat);
    }
    dcell *IMCC=dcellnew(1,1);
    IMCC->p[0]=dref(ngsmod->IMCC_TT);
    dcell *Ptt=dcellnew(ndm,1);
    for(int idm=0; idm<ndm; idm++){
	dmm(&(Ptt->p[idm]),0,IMCC->p[0],HatWHmt->p[idm],"nt",1);
    }
    dcellfree(IMCC);
    dcellfree(modc);
    if(use_ploc){
	free(amp);
    }
    dcellfree(HatWHmt);
    return Ptt;
}
/**
   DM modes for all the low order modes, defined on DM grid. It uses ngsmod2dm
   to define the modes*/

static dcell *ngsmod_dm(const PARMS_T *parms, RECON_T *recon){
    NGSMOD_T *ngsmod=recon->ngsmod; 
    int ndm=parms->ndm;
    int nmod=ngsmod->nmod;
    dcell *M=dcellnew(1,1);
    M->p[0]=dnew(nmod,1);
    dcell *mod=dcellnew(ndm,1);
    dcell *dmt=dcellnew(ndm,1);
    loc_t **aloc=recon->aloc->p;
    for(int idm=0; idm<ndm; idm++){
	dmt->p[idm]=dnew(aloc[idm]->nloc,1);
	mod->p[idm]=dnew(aloc[idm]->nloc,nmod);
    }
    for(int imod=0; imod<nmod; imod++){
	dcellzero(dmt);
	M->p[0]->p[imod]=1;
	ngsmod2dm(&dmt,recon,M,1.);
	M->p[0]->p[imod]=0;
	for(int idm=0; idm<ndm; idm++){
	    memcpy(PCOL(mod->p[idm], imod), dmt->p[idm]->p, 
		   sizeof(double)*aloc[idm]->nloc);
	}
    }
    dcellfree(M);
    dcellfree(dmt);
    return mod;
}

/**
   Compute NGS modes Ha*M in the science directon using ray tracing. Not used
*/
/*dcell *ngsmod_hm_accphi(const PARMS_T *parms, RECON_T *recon, const APER_T *aper){
    NGSMOD_T *ngsmod=recon->ngsmod;
    loc_t **aloc=recon->aloc->p;
    const int ndm=parms->ndm;
    dcell *dmt=dcellnew(ndm,1);
    for(int idm=0; idm<ndm; idm++){
	dmt->p[idm]=dnew(aloc[idm]->nloc,1);
    }
    dcell *M=dcellnew(1,1);
    const int nmod=ngsmod->nmod;
    M->p[0]=dnew(nmod,1);
    dcell *HMC=dcellnew(parms->evl.nevl,nmod);
    dcell* HM=HMC;
    int nloc=aper->locs->nloc;
    for(int imod=0; imod<nmod; imod++){
	dzero(M->p[0]);
	M->p[0]->p[imod]=1;
	dcellzero(dmt);
	ngsmod2dm(&dmt,recon,M,1.);
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    IND(HM,ievl,imod)=dnew(nloc,1);
	    for(int idm=0; idm<parms->ndm; idm++){
		double ht=parms->dm[idm].ht;
		double displacex=parms->evl.thetax->p[ievl]*ht;
		double displacey=parms->evl.thetay->p[ievl]*ht;
		prop_nongrid(aloc[idm], dmt->p[idm]->p,
			     aper->locs, IND(HM,ievl,imod)->p,1,
			     displacex, displacey, 1.,0,0);
	    }
	}
    }
    return HMC;
}*/
/**
   Compute NGS modes Ha*M in the science directon using analytic formula. Not used
   confirmed to agree with ngsmod_hm_accphi except DM artifacts 
*/
/*
dcell *ngsmod_hm_ana(const PARMS_T *parms, RECON_T *recon, const APER_T *aper){
    NGSMOD_T *ngsmod=recon->ngsmod;
    const double hs=ngsmod->hs;
    const double ht=ngsmod->ht;
    const double scale=pow(1.-ht/hs,-2);
    const double scale1=1.-scale;
    const int nmod=ngsmod->nmod;
    const double MCC_fcp=recon->ngsmod->aper_fcp;
    dcell *HMC=dcellnew(parms->evl.nevl,nmod);
    dcell* HM=HMC;
    double *x=aper->locs->locx;
    double *y=aper->locs->locy;
    int nloc=aper->locs->nloc;
    for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	double sx=parms->evl.thetax->p[ievl];
	double sy=parms->evl.thetay->p[ievl];
	for(int imod=0; imod<nmod; imod++){
	    IND(HM,ievl,imod)=dnew(nloc,1);
	}
	if(nmod==2){
	    for(int iloc=0; iloc<nloc; iloc++){
		IND(HM,ievl,0)->p[iloc]=x[iloc];
		IND(HM,ievl,1)->p[iloc]=y[iloc];
	    }
	}else{
	    for(int iloc=0; iloc<nloc; iloc++){
		double xx=pow(x[iloc],2);
		double yy=pow(y[iloc],2);
		double xy=x[iloc]*y[iloc];
		IND(HM,ievl,0)->p[iloc]=x[iloc];
		IND(HM,ievl,1)->p[iloc]=y[iloc];
		if(ngsmod->indps){
		    if(ngsmod->ahstfocus){
			IND(HM,ievl,ngsmod->indps)->p[iloc]=
			    -2*ht*(sx*x[iloc]+sy*y[iloc])*scale;
		    }else{
			IND(HM,ievl,ngsmod->indps)->p[iloc]=(xx+yy-MCC_fcp)*scale1
			    -2*ht*(sx*x[iloc]+sy*y[iloc])*scale;
		    }
		    IND(HM,ievl,ngsmod->indps+1)->p[iloc]=(xx-yy)*scale1-2*ht*(sx*x[iloc]-sy*y[iloc])*scale;
		    IND(HM,ievl,ngsmod->indps+2)->p[iloc]=(xy)*scale1-(sx*y[iloc]+sy*x[iloc])*ht*scale;
		}
		if(ngsmod->indastig){
		    IND(HM,ievl,ngsmod->indastig)->p[iloc]=(xx-yy);
		    IND(HM,ievl,ngsmod->indastig+1)->p[iloc]=(xy);
		}
		if(ngsmod->indfocus){
		    IND(HM,ievl,ngsmod->indfocus)->p[iloc]=(xx+yy-MCC_fcp);
		}
	    }
	}
    }
    return HMC;
}
*/
/**
   AHST parameters that are related to the geometry only, and
   will not be updated when estimated WFS measurement noise changes.
*/
void setup_ngsmod_prep(const PARMS_T *parms, RECON_T *recon, 
		       const APER_T *aper, const POWFS_T *powfs){
    if(recon->ngsmod) error("Should only be called once\n");
    NGSMOD_T *ngsmod=recon->ngsmod=mycalloc(1,NGSMOD_T);
    ngsmod->ahstfocus=parms->sim.ahstfocus;
    const int ndm=parms->ndm;	
    ngsmod->aper_fcp=aper->fcp;
    if(ndm>1 && fabs(parms->dm[0].ht)>1.e-10){
	error("Error configuration. First DM is not on ground\n");
    }
    double hs=NAN;
    for(int ipowfs=0; ipowfs<parms->npowfs; ipowfs++){
	if(parms->powfs[ipowfs].lo || parms->powfs[ipowfs].skip) continue;
	if(isnan(hs)){
	    hs=parms->powfs[ipowfs].hs;
	}else{
	    if(isfinite(hs) || isfinite(parms->powfs[ipowfs].hs)){
		if(fabs(hs-parms->powfs[ipowfs].hs)>1000){
		    //High order GS at different altitude.
		    hs=INFINITY; break;
		}
	    }
	}
    }
    ngsmod->nmod=2; //Basic tip/tilt mode.
    if(isfinite(hs)){
	//LGS WFS.
	if(ndm==2){//Plate scale mode
	    ngsmod->indps=ngsmod->nmod;
	    ngsmod->nmod+=3;
	}else if(parms->nhiwfs>1){//LGS AO
	    ngsmod->indastig=ngsmod->nmod;
	    ngsmod->nmod+=2;
	}
	if(parms->sim.mffocus || ngsmod->indastig){
	    ngsmod->indfocus=ngsmod->nmod;
	    ngsmod->nmod+=1;
	}
    }
    info2("ngsmod nmod=%d, ahstfocus=%d\n", ngsmod->nmod, parms->sim.ahstfocus);
    ngsmod->hs=hs;
    if(ndm>1){
	ngsmod->ht=parms->dm[ndm-1].ht;//last DM.
    }else{
	ngsmod->ht=0;
    }
    ngsmod->scale=pow(1.-ngsmod->ht/ngsmod->hs,-2);
    /*modal cross coupling matrix along science*/
    ngsmod->MCCP=ngsmod_mcc(parms,recon,aper, parms->evl.wt->p);
    if(ngsmod->MCCP->nx==1){
	ngsmod->MCC=dref(ngsmod->MCCP->p[0]);
    }else{
	ngsmod->MCC=NULL;
	for(int ievl=0; ievl<parms->evl.nevl; ievl++){
	    dadd(&ngsmod->MCC, 1, ngsmod->MCCP->p[ievl], parms->evl.wt->p[ievl]);
	}
    }
    dmat *MCCl=dchol(ngsmod->MCC);
    ngsmod->MCCu=dtrans(MCCl); dfree(MCCl);
    if(parms->save.setup){
	writebin(recon->ngsmod->MCC, "ahst_MCC");
	writebin(recon->ngsmod->MCCu, "ahst_MCCu");
    }
    ngsmod->IMCC=dinvspd(ngsmod->MCC);
    dmat *MCC=ngsmod->MCC;
    ngsmod->IMCC_TT=dnew(2,2);
    ngsmod->IMCC_TT->p[0]=IND(MCC,0,0);
    ngsmod->IMCC_TT->p[1]=IND(MCC,1,0);
    ngsmod->IMCC_TT->p[2]=IND(MCC,0,1);
    ngsmod->IMCC_TT->p[3]=IND(MCC,1,1);
    dinvspd_inplace(ngsmod->IMCC_TT);
    /*the ngsmodes defined on the DM.*/
    ngsmod->Modes=ngsmod_dm(parms,recon);
    if(recon->actstuck && !parms->recon.modal){
	warning2("Apply stuck actuators to ngs modes\n");
	act_zero(recon->aloc, recon->ngsmod->Modes, recon->actstuck);
    }
   /*if(recon->actfloat){
      We do extrapolation to float actuators, so no need to modify Pngs/Ptt.
      warning2("Apply float actuators to Pngs, Ptt\n");
      act_zero(recon->aloc, recon->ngsmod->Modes, recon->actfloat);
      }*/

    if(parms->recon.split==1 && !parms->tomo.ahst_idealngs && parms->ntipowfs){
	ngsmod->GM=dcellnew(parms->nwfsr, 1);
	int nttwfs=0;
	int nttfwfs=0;
	info2("Low order control includes WFS");
	for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
	    int ipowfs=parms->wfsr[iwfs].powfs;
	    if(parms->powfs[ipowfs].skip==3) continue;
	    if(parms->powfs[ipowfs].lo
	       || (parms->recon.split && parms->nlopowfs==0 && !parms->powfs[ipowfs].trs)){
		if(parms->powfs[ipowfs].order==1){
		    nttwfs++;
		}else{
		    nttfwfs++;
		}
		info2(" %d", iwfs);
		for(int idm=0; idm<parms->ndm; idm++){
		    if(parms->powfs[ipowfs].type==0){//shwfs
			dspmm(PIND(ngsmod->GM, iwfs), IND(recon->GAlo, iwfs, idm), IND(ngsmod->Modes, idm), "nn", 1);
		    }else{//pwfs.
			double  ht = parms->dm[idm].ht-parms->powfs[ipowfs].hc;
			double  scale=1. - ht/parms->wfs[iwfs].hs;
			double  dispx=0, dispy=0;
			dispx=parms->wfsr[iwfs].thetax*ht;
			dispy=parms->wfsr[iwfs].thetay*ht;
			dmat *tmp=pywfs_mkg(powfs[ipowfs].pywfs, recon->aloc->p[idm],
					    IND(ngsmod->Modes, idm), 0, dispx, dispy, scale);
			dadd(PIND(ngsmod->GM, iwfs), 1, tmp, 1);//accumulate
		    }
		}
	    }
	}
	info2("\n");
	if(ngsmod->nmod>2 && nttfwfs==0){
	    error("Only TT wfs cannot control plate scale or focus\n");
	}
	if(ngsmod->nmod==6 && nttfwfs==1 && nttwfs==0){
	    warning("There is only one wfs, remove first plate scale mode as it degenerates with focus mode");
	    for(int iwfs=0; iwfs<parms->nwfsr; iwfs++){
		if(ngsmod->GM->p[iwfs]){
		    int nx=ngsmod->GM->p[iwfs]->nx;
		    memset(ngsmod->GM->p[iwfs]->p+nx*2, 0, nx*sizeof(double));
		}
	    }
	}
    }
    if(parms->recon.modal){//convert Modes to space of amod
	for(int idm=0; idm<parms->ndm; idm++){
	    dmat *proj=dpinv(recon->amod->p[idm], NULL);
	    dmat *tmp=0;
	    dmm(&tmp, 0, proj, ngsmod->Modes->p[idm], "nn", 1);
	    dfree(ngsmod->Modes->p[idm]);
	    ngsmod->Modes->p[idm]=tmp;
	    dfree(proj);
	}
    }
    /**
       parms->tomo.ahst_wt control NGS modes removal from LGS DM commands
       if ahst_wt==1
       Rngs*GA*dmerr is zero
       if ahst_wt==2
       Doesn't perturb NGS modes in science direction.
       if ahst_wt==3
       Identity weighting.

    */
    if(parms->tomo.ahst_wt==1){
	//Do it in setup_ngsmod_recon();
    }else if(parms->tomo.ahst_wt==2){
	/*Use science based weighting. */
	if(parms->dbg.wamethod==0){
	    info("Wa using DM mode\n");

	    tic;
	    ngsmod->Wa=ngsmod_Wa(parms,recon,aper,0);
	    /*
	      Add tikhonov regularization. H is from aloc to some other loc. 
	      the eigen value of H'*amp*H is about 4/aloc->nloc.
	    */
	    int nact=0;
	    for(int idm=0; idm<parms->ndm; idm++){
		nact+=recon->aloc->p[idm]->nloc;
	    }
	    double maxeig=4./nact;
	    dcelladdI(ngsmod->Wa, 1e-9*maxeig);
	    
	    toc("Wa");
	    ngsmod->Pngs=dcellpinv(ngsmod->Modes,ngsmod->Wa);
	    toc("Pngs");
	}else{
	    info("Wa using science mode\n");
	    tic;
	    ngsmod->Pngs=ngsmod_Pngs_Wa(parms,recon,aper,0);
	    toc("Pngs_Wa");
	}
    }else if(parms->tomo.ahst_wt==3){/*Identity weighting. */
	ngsmod->Pngs=dcellpinv(ngsmod->Modes, NULL);
    }else{
	error("Invalid parms->tomo.ahst_wt=%d\n", parms->tomo.ahst_wt);
    }
    if(parms->save.setup){
	writebin(ngsmod->Modes, "ahst_Modes");
	writebin(ngsmod->GM,  "ahst_GM");
	if(ngsmod->Pngs) writebin(ngsmod->Pngs,"ahst_Pngs");
    }
}
    
/**
   setup NGS modes reconstructor in ahst mode.
 */
void setup_ngsmod_recon(const PARMS_T *parms, RECON_T *recon){
    NGSMOD_T *ngsmod=recon->ngsmod;
    if(parms->recon.split==1 && !parms->tomo.ahst_idealngs && parms->ntipowfs){
	cellfree(ngsmod->Rngs);
	/*
	  W is recon->saneai;
	  Rngs=(M'*G'*W*G*M)^-1*M'*G'*W
	  Pngs=Rngs*GA
	*/
	ngsmod->Rngs=dcellpinv(ngsmod->GM, recon->saneai);
    }
  
    if(parms->tomo.ahst_wt==1){
	/*Use gradient weighting. */
	dcellzero(ngsmod->Pngs);
	dcellmm(&ngsmod->Pngs, ngsmod->Rngs, recon->GAlo, "nn", 1);
	if(parms->save.setup){
	    writebin(ngsmod->Pngs,"ahst_Pngs");
	}
    }
 
    if(parms->save.setup){
	writebin(ngsmod->Rngs,"ahst_Rngs");
    }
}
/**
   used in performance evaluation on science opds. accumulate to out*/
void calc_ngsmod_dot(double *pttr_out, double *pttrcoeff_out,
		     double *ngsmod_out,
		     const PARMS_T *parms, const NGSMOD_T *ngsmod,
		     const APER_T *aper, const double *opd, int ievl){
    const double *restrict amp=aper->amp->p;
    const double *restrict locx=aper->locs->locx;
    const double *restrict locy=aper->locs->locy;
    double coeff[6]={0,0,0,0,0,0};
    double tot=0;
    const int nmod=ngsmod->nmod;
    if(nmod==2){
	for(int iloc=0; iloc<aper->locs->nloc; iloc++){
	    const double junk=amp[iloc]*opd[iloc];
	    tot+=junk*opd[iloc];
	    const double junkx=locx[iloc]*junk;
	    const double junky=locy[iloc]*junk;
	    coeff[0]+=junk;
	    coeff[1]+=junkx;
	    coeff[2]+=junky;
	}
    }else{
	for(int iloc=0; iloc<aper->locs->nloc; iloc++){
	    const double junk=amp[iloc]*opd[iloc];
	    tot+=junk*opd[iloc];
	    coeff[0]+=junk;//piston
	    const double junkx=locx[iloc]*junk;
	    const double junky=locy[iloc]*junk;
	    coeff[1]+=junkx;//tip
	    coeff[2]+=junky;//tilt
	    coeff[3]+=locx[iloc]*junkx;//xx
	    coeff[4]+=locy[iloc]*junky;//yy
	    coeff[5]+=locx[iloc]*junky;//xy
	}
    }
    const double thetax=parms->evl.thetax->p[ievl]; 
    const double thetay=parms->evl.thetay->p[ievl]; 
    calc_ngsmod_post(pttr_out, pttrcoeff_out, ngsmod_out, tot, coeff, ngsmod, aper, thetax, thetay);
}
/**
   Separate post processing part so that GPU code can call it.
*/
void calc_ngsmod_post(double *pttr_out, double *pttrcoeff_out, double *ngsmod_out,
		      double tot, const double *coeff,  const NGSMOD_T *ngsmod, 
		      const APER_T *aper,double thetax, double thetay){
    const double MCC_fcp=ngsmod->aper_fcp;
    const double ht=ngsmod->ht;
    const double scale=ngsmod->scale;
    const double scale1=1.-scale;

    if(pttrcoeff_out){
	memset(pttrcoeff_out, 0, sizeof(double)*3);
	dmulvec(pttrcoeff_out, aper->imcc, coeff, 1);
    }
    if(pttr_out){
	/*compute TT removed wavefront variance as a side product */
	double pis=aper->ipcc*coeff[0]*coeff[0];
	double ptt=dwdot3(coeff, aper->imcc, coeff);
	pttr_out[0]=tot-pis;/*PR */
	pttr_out[1]=ptt-pis;/*TT */
	pttr_out[2]=tot-ptt;/*PTTR */
    }
    /*don't use +=. need locking */
    ngsmod_out[0]=coeff[1];
    ngsmod_out[1]=coeff[2];

    if(ngsmod->indps){
	if(ngsmod->ahstfocus){
	    ngsmod_out[ngsmod->indps]=(-2*scale*ht*(thetax*coeff[1]+thetay*coeff[2]));
	}else{
	    ngsmod_out[ngsmod->indps]=(scale1*(coeff[3]+coeff[4]-coeff[0]*MCC_fcp)
				-2*scale*ht*(thetax*coeff[1]+thetay*coeff[2]));
	}
	ngsmod_out[ngsmod->indps+1]=(scale1*(coeff[3]-coeff[4])
			    -2*scale*ht*(thetax*coeff[1]-thetay*coeff[2]));
	ngsmod_out[ngsmod->indps+2]=(scale1*(coeff[5])
			    -scale*ht*(thetay*coeff[1]+thetax*coeff[2]));
    }
    if(ngsmod->indastig){
	ngsmod_out[ngsmod->indastig]=(coeff[3]-coeff[4]);
	ngsmod_out[ngsmod->indastig+1]=(coeff[5]);
    }
    if(ngsmod->indfocus){
	ngsmod_out[ngsmod->indfocus]=(coeff[3]+coeff[4]-coeff[0]*MCC_fcp);
    }
}
/**
   Convert NGS modes to DM actuator commands using analytical expression. For >2 DMs, we only put NGS modes on ground and top-most DM.
*/
void ngsmod2dm(dcell **dmc, const RECON_T *recon, const dcell *M, double gain){
    if(!M || !M->p[0]) return;
    const NGSMOD_T *ngsmod=recon->ngsmod;
    const int nmod=ngsmod->nmod;
    assert(M->nx==1 && M->ny==1 && M->p[0]->nx==nmod);
    double scale=ngsmod->scale;
    /*The MCC_fcp depends weakly on the aperture sampling. */
    double MCC_fcp=ngsmod->aper_fcp;
    loc_t **aloc=recon->aloc->p;
    /*convert mode vector and add to dm commands */
    const int ndm=recon->aloc->nx;
    if(!*dmc){
	*dmc=dcellnew(ndm,1);
    }
    for(int idm=0; idm<ndm; idm++){
	if(!(*dmc)->p[idm]){
	    (*dmc)->p[idm]=dnew(recon->aloc->p[idm]->nloc, 1);
	}
    }

    /*first dm */
    double *pm=M->p[0]->p;
  
    for(int idm=0; idm<ndm; idm++){
	double *p=(*dmc)->p[idm]->p;
	unsigned long nloc=aloc[idm]->nloc;
	double *xloc=aloc[idm]->locx;
	double *yloc=aloc[idm]->locy;
	if(idm==0){
	    if(nmod==2){//tip/tilt only
		for(unsigned long iloc=0; iloc<nloc; iloc++){
		    p[iloc]+=gain*(pm[0]*xloc[iloc]+pm[1]*yloc[iloc]);
		}
	    }else{//with platescale and/or focus
		double focus=0, astigx=0, astigy=0;
		if(ngsmod->indfocus){//focus is a mode
		    focus+=pm[ngsmod->indfocus];
		}
		if(ngsmod->indps){//ps mode
		    if(ngsmod->ahstfocus){
			focus+=pm[ngsmod->indps]*scale;//scaled to avoid cause focus mode in science.
		    }else{
			focus+=pm[ngsmod->indps];
		    }
		    astigx+=pm[ngsmod->indps+1];
		    astigy+=pm[ngsmod->indps+2];
		}
		if(ngsmod->indastig){
		    astigx+=pm[ngsmod->indastig];
		    astigy+=pm[ngsmod->indastig+1];
		}
		for(unsigned long iloc=0; iloc<nloc; iloc++){
		    double xx=xloc[iloc]*xloc[iloc];
		    double xy=xloc[iloc]*yloc[iloc];
		    double yy=yloc[iloc]*yloc[iloc];
		    p[iloc]+=gain*(pm[0]*xloc[iloc]
				   +pm[1]*yloc[iloc]
				   +focus*(xx+yy-MCC_fcp)
				   +astigx*(xx-yy)
				   +astigy*(xy));
		}
	    }
	}else if(idm+1==ndm && ngsmod->indps){
	    double scale2=-scale*gain;
	    for(unsigned long iloc=0; iloc<nloc; iloc++){
		double xx=xloc[iloc]*xloc[iloc];
		double xy=xloc[iloc]*yloc[iloc];
		double yy=yloc[iloc]*yloc[iloc];
		p[iloc]+=scale2*(pm[ngsmod->indps]*(xx+yy-MCC_fcp)
				 +pm[ngsmod->indps+1]*(xx-yy)
				 +pm[ngsmod->indps+2]*(xy));
	    }
	}	
    }
}
/**
   Convert NGS mode vector to aperture grid for science directions.  */
void ngsmod2science(dmat *iopd, const loc_t *loc, const NGSMOD_T *ngsmod, 
		    double thetax, double thetay,
		    const double *mod, double alpha){
    const double *locx=loc->locx;
    const double *locy=loc->locy;
    const int nmod=ngsmod->nmod;
    if(nmod==2){//tip/tilt only
	for(int iloc=0; iloc<loc->nloc; iloc++){
	    double tmp=locx[iloc]*mod[0]+locy[iloc]*mod[1];
	    iopd->p[iloc]+=tmp*alpha;
	}
    }else{
	const double ht=ngsmod->ht;
	const double scale=ngsmod->scale;
	const double scale1=1.-scale;
	double focus=0, ps1=0, ps2=0, ps3=0,astigx=0,astigy=0;
	if(ngsmod->indfocus){
	    focus+=mod[ngsmod->indfocus];
	}
	if(ngsmod->indps){
	    if(!ngsmod->ahstfocus){
		focus+=mod[ngsmod->indps]*scale1;
	    }
	    ps1=mod[ngsmod->indps];
	    ps2=mod[ngsmod->indps+1];
	    ps3=mod[ngsmod->indps+2];
	}
	if(ngsmod->indastig){
	    astigx=mod[ngsmod->indastig];
	    astigy=mod[ngsmod->indastig+1];
	}
	for(int iloc=0; iloc<loc->nloc; iloc++){
	    double x=locx[iloc];
	    double y=locy[iloc];
	    double xy=x*y;
	    double x2=x*x;
	    double y2=y*y;
	    double tmp= locx[iloc]*mod[0]
		+locy[iloc]*mod[1]
		+focus*(x2+y2)
		+ps1*(-2*scale*ht*(thetax*x+thetay*y))
		+ps2*((x2-y2)*scale1 - 2*scale*ht*(thetax*x-thetay*y))
		+ps3*(xy*scale1-scale*ht*(thetay*x+thetax*y))
		+astigx*(x2-y2)
		+astigy*(xy);
	    iopd->p[iloc]+=tmp*alpha;
	}
    }
}
void ngsmod_free(NGSMOD_T *ngsmod){
    if(!ngsmod) return;
    dcellfree(ngsmod->GM);
    dcellfree(ngsmod->Rngs);
    dcellfree(ngsmod->Pngs);
    dcellfree(ngsmod->Modes);
    dfree(ngsmod->MCC);
    dfree(ngsmod->MCCu);
    dcellfree(ngsmod->MCCP);
    dfree(ngsmod->IMCC);
    dfree(ngsmod->IMCC_TT);
    free(ngsmod);
}

/**
   remove NGS modes from LGS DM commands
   if nmod==6: make sure the global focus mode is not removed from LGS result.
*/
void remove_dm_ngsmod(SIM_T *simu, dcell *dmerr){
    if(!dmerr) return;
    const RECON_T *recon=simu->recon;
    const NGSMOD_T *ngsmod=recon->ngsmod;

    dcellzero(simu->Mngs);
    dcellmm(&simu->Mngs, ngsmod->Pngs, dmerr, "nn",1);
    double *mngs=simu->Mngs->p[0]->p;
    if(ngsmod->indastig){
	//Temporary: remove LPF'ed focus/astigmatism from LGS DM command.
	if(!simu->ngsmodlpf){
	    simu->ngsmodlpf=dnew(3,1);
	}
	const double lpfocus=simu->parms->sim.lpfocushi;
	warning_once("HPF focus/astig from DM error signal. lpfocus=%g\n", lpfocus);
	simu->ngsmodlpf->p[0]=simu->ngsmodlpf->p[0]*(1-lpfocus)+mngs[ngsmod->indfocus]*lpfocus;
	simu->ngsmodlpf->p[1]=simu->ngsmodlpf->p[1]*(1-lpfocus)+mngs[ngsmod->indastig]*lpfocus;
	simu->ngsmodlpf->p[2]=simu->ngsmodlpf->p[2]*(1-lpfocus)+mngs[ngsmod->indastig+1]*lpfocus;
	mngs[ngsmod->indfocus]=simu->ngsmodlpf->p[0];
	mngs[ngsmod->indastig]=simu->ngsmodlpf->p[1];
	mngs[ngsmod->indastig+1]=simu->ngsmodlpf->p[2];
    }
    if(ngsmod->indfocus && simu->parms->sim.mffocus){
	//zero out global focus mode if any so we don't remove it from DM commands.
	if(ngsmod->indps && ngsmod->ahstfocus){
	    const double scale=ngsmod->scale;
	    mngs[ngsmod->indfocus]=mngs[2]*(1-scale);
	}else{
	    mngs[ngsmod->indfocus]=0;
	}
    }
    dcellmm(&dmerr, ngsmod->Modes, simu->Mngs, "nn", -1);
}
