/*
  Copyright 2009, 2010, 2011 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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
#include "skyc.h"
#include "parms.h"
#include "types.h"
#include "skysim_utils.h"
/**
   \file skyc/skysim_utils.c
   Utilities for skysim.c
*/
/**
   add photon and read out noise.  pcaclib part of bkgrnd is calibrated out.
   set to 1 usually.  */
static void addnoise(dmat *A, rand_t* rstat, 
		     const double bkgrnd, const double pcalib, 
		     const double rne){
    for(int ix=0; ix<A->nx*A->ny; ix++){
	A->p[ix]=randp(rstat,A->p[ix]+bkgrnd)
	    -bkgrnd*pcalib+rne*randn(rstat);
    }
}

/**
   convert mod vector to ngs WFS cloc and add to it's opd or complex pupil function.
   Input/Output: wvf (complex pupil function).
   Input:
   wvl: wavelength. 
   modm: the mode vector
   cloc: the subaperture coarse sampled loc.
   fpc: focus piston coupling. (mean(x^2+y^2) to eliminate piston term in focus)
   thetax, thetay: direction of WFS.
   parms: parms
*/
void ngsmod2wvf(cmat *wvf,            /**<[in/out] complex pupil function*/
		double wvl,           /**<[in] the wavelength*/
		const dmat *modm,     /**<[in] the NGS mode vector*/
		const loc_t *cloc,    /**<[in] the subaperture coarse sampled loc*/
		double fpc,           /**<[in] focus piston coupling.
					 (mean(x^2+y^2) to eliminate piston term
					 in focus)*/
		double thetax,        /**<[in] direction of WFS*/
		double thetay,        /**<[in] direction of WFS*/
		const PARMS_S *parms  /**<[in] the parms*/
		){
    const double *locx=cloc->locx;
    const double *locy=cloc->locy;
    const long nloc=cloc->nloc;
    const double *mod=modm->p;
    const dcomplex ik=2*M_PI/wvl*I;
    if(modm->nx==2){
	for(int iloc=0; iloc<nloc; iloc++){
	    double tmp=locx[iloc]*mod[0]+locy[iloc]*mod[1];
	    wvf->p[iloc]*=cexp(ik*tmp);
	}
    }else{
	const double hc=parms->maos.hc;
	const double hs=parms->maos.hs;
	const double scale=pow(1.-hc/hs, -2);
	const double scale1=1.-scale;
	for(int iloc=0; iloc<nloc; iloc++){
	    double x=locx[iloc];
	    double y=locy[iloc];
	    double xy=x*y;
	    double x2=x*x;
	    double y2=y*y;
	    double tmp= locx[iloc]*mod[0]
		+locy[iloc]*mod[1]
		+mod[2]*((x2+y2-fpc)*scale1-2*scale*hc*(thetax*x+thetay*y))
		+mod[3]*((x2-y2)*scale1 - 2*scale*hc*(thetax*x-thetay*y))
		+mod[4]*(xy*scale1-scale*hc*(thetay*x+thetax*y));
	    wvf->p[iloc]*=cexp(ik*tmp);
	}
    }
}

/**
   ztilt time domain noise free simulation at LGS sampling freq.
   
   Written: 2010-06-09
   Tested OK: 2010-06-10
*/
dcell* skysim_ztilt(dmat *mideal, ASTER_S *aster, const PARMS_S *parms){
    const double gain=parms->skyc.intgain;
    const int nmod=mideal->nx;
    dmat *mint=NULL;        //integrator output.
    dmat *mreal=NULL;       //modal correction at this step.
    dmat *merr=dnew(nmod,1);//modal error
    dmat *merrm=NULL;       //measured model error
    PDMAT(mideal,pmideal);
    dcell *mres=dcellnew(1,aster->nstep);//residual mode.
    dmat *grad=NULL;
    dmat *pgm=aster->pgm->p[aster->pgm->nx-1];
    for(int istep=0; istep<aster->nstep; istep++){
	memcpy(merr->p, pmideal[istep], nmod*sizeof(double));
	dadd(&merr, 1, mreal, -1);//form error;
	dcp(&mres->p[istep],merr);//record error.
	dzero(grad);
	dmm(&grad, aster->gm, merr, "nn", 1);
	int itsa=0;
	for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
	    const int ipowfs=aster->wfs[iwfs].ipowfs;
	    const long nsa=parms->maos.nsa[ipowfs];
	    for(long isa=0; isa<nsa*2; isa++){
		grad->p[isa+itsa]+=aster->wfs[iwfs].ztiltout->p[istep]->p[isa];
	    }
	    itsa+=nsa*2;
	}
	dzero(merrm);
	dmm(&merrm, pgm, grad, "nn", 1);
	
	//Servo filter.
	dcp(&mreal, mint);//integrator output from last step.
	dadd(&mint, 1, merrm, gain);//integrator;
    }
    dfree(mint); dfree(mreal); dfree(merr); dfree(merrm);dfree(grad);
    return mres;
}

/**
   Time domain physical simulation.
   
   noisy: 
       - 0: no noise at all; 
       - 1: poisson and read out noise. 
       - 2: only poisson noise.   
*/
dmat *skysim_phy(dmat **mresout,SIM_S *simu, ASTER_S *aster, POWFS_S *powfs, 
		 const PARMS_S *parms, int idtrat, int noisy, int demotettf){
    int dtrat=parms->skyc.dtrats[idtrat];
    const int nmod=simu->mideal->nx;
    PDMAT(simu->mideal,pmideal);
    PDMAT(simu->mideal_oa, pmideal_oa);
    dmat *res=dnew(6,1);/*Results. 1-2: NGS and TT modes., 
			  3-4:On axis NGS and TT modes,
			  4-6: On axis NGS and TT wihtout considering un-orthogonality.*/
    dmat *mreal=NULL;//modal correction at this step.
    dmat *merr=dnew(nmod,1);//modal error
    dmat *merrm=NULL;//measured model error
    
    dmat *mres=dnew(nmod,aster->nstep);
    PDMAT(mres,pmres);
    PDMAT(parms->skyc.rnefs,rnefs);
    dmat *grad=dnew(aster->tsa*2,1);
    dmat *zgrad=dnew(aster->tsa*2,1);
    dmat *grads=dnew(aster->tsa*2,aster->nstep);
    dmat *zgrads=dnew(aster->tsa*2,aster->nstep);
    dcell **psf=calloc(aster->nwfs, sizeof(dcell*));
    ccell *wvf=ccellnew(aster->nwfs,1);
    ccell *wvfc=ccellnew(aster->nwfs,1);
    dcell **mtche=calloc(aster->nwfs, sizeof(dcell*));
    dcell **ints=calloc(aster->nwfs, sizeof(dcell*));
    ccell *otf=ccellnew(aster->nwfs,1);
    SERVO_T *st2t=calloc(1, sizeof(SERVO_T));
    double ngsol=simu->rmsol->p[0];
    const double dtngs=parms->maos.dt*dtrat;
    const long nwvl=parms->maos.nwvl;
    dmat *pgm;
    if(demotettf){
	pgm=aster->pgmtt->p[idtrat];
	if(!pgm){
	    warning("Demoting is forbidden\n");
	    pgm=aster->pgm->p[idtrat];
	}
    }else{
	pgm=aster->pgm->p[idtrat];
    }
    for(long iwfs=0; iwfs<aster->nwfs; iwfs++){
	const int ipowfs=aster->wfs[iwfs].ipowfs;
	const long ncomp=parms->maos.ncomp[ipowfs];
	const long nsa=parms->maos.nsa[ipowfs];
	wvf->p[iwfs]=cnew(ncomp,ncomp);
	wvfc->p[iwfs]=cnew(ncomp/2,ncomp/2);
	psf[iwfs]=dcellnew(nsa,nwvl);
	cfft2plan(wvf->p[iwfs], -1);
	mtche[iwfs]=aster->wfs[iwfs].pistat->mtche[idtrat];
	otf->p[iwfs]=cnew(ncomp,ncomp);
	cfft2plan(otf->p[iwfs],-1);
	cfft2plan(otf->p[iwfs],1);
	ints[iwfs]=dcellnew(nsa,1);
	int pixpsa=parms->skyc.pixpsa[ipowfs];
	for(long isa=0; isa<nsa; isa++){
	    ints[iwfs]->p[isa]=dnew(pixpsa,pixpsa);
	}

    }
    int phycount=0;
    for(int istep=0; istep<aster->nstep; istep++){
	memcpy(merr->p, pmideal[istep], nmod*sizeof(double));
	dadd(&merr, 1, mreal, -1);//form NGS mode error;

	if(istep>=parms->skyc.evlstart){
	    double res_ngs=dwdot(merr->p,parms->maos.mcc,merr->p);
	    if(res_ngs>ngsol*5){
		//warning2("%5.1f Hz: loop diverged. \n", parms->skyc.fss[idtrat]);
		dfree(res); res=NULL;
		break;
	    }
	    {
		res->p[0]+=res_ngs;
		res->p[1]+=dwdot2(merr->p,parms->maos.mcc_tt,merr->p);
		double dot_oa=dwdot(merr->p, parms->maos.mcc_oa, merr->p);
		double dot_res_ideal=dwdot(merr->p, parms->maos.mcc_oa, pmideal[istep]);
		double dot_res_oa=0;
		for(int imod=0; imod<nmod; imod++){
		    dot_res_oa+=merr->p[imod]*pmideal_oa[istep][imod];
		}
		res->p[2]+=dot_oa-2*dot_res_ideal+2*dot_res_oa;
		res->p[4]+=dot_oa;
	    }
	    {
		double dot_oa_tt=dwdot2(merr->p, parms->maos.mcc_oa_tt, merr->p);
		//Notice that mcc_oa_tt2 is 2x5 marix.
		double dot_res_ideal_tt=dwdot(merr->p, parms->maos.mcc_oa_tt2, pmideal[istep]);
		double dot_res_oa_tt=0;
		for(int imod=0; imod<2; imod++){
		    dot_res_oa_tt+=merr->p[imod]*pmideal_oa[istep][imod];
		}
		res->p[3]+=dot_oa_tt-2*dot_res_ideal_tt+2*dot_res_oa_tt;
		res->p[5]+=dot_oa_tt;
	    }
	}
	memcpy(pmres[istep],merr->p,sizeof(double)*nmod);
	
	if(istep<parms->skyc.phystart){
	    //Ztilt, noise free simulation for acquisition.
	    if(istep % dtrat == 0){
		dzero(zgrad);
	    }
	    dmm(&zgrad, aster->gm, merr, "nn", 1);//grad due to residual NGS mode.
	    int itsa=0;
	    for(int iwfs=0; iwfs<aster->nwfs; iwfs++){
		const int ipowfs=aster->wfs[iwfs].ipowfs;
		const long nsa=parms->maos.nsa[ipowfs];
		for(long isa=0; isa<nsa*2; isa++){//add ztilt.
		    zgrad->p[isa+itsa]+=aster->wfs[iwfs].ztiltout->p[istep]->p[isa];
		}
		itsa+=nsa*2;
	    }
	    dcp(&mreal, st2t->mint);
	    if((istep+1) % dtrat == 0){//has output
		dscale(zgrad, 1./dtrat);//averaging gradients.
		dzero(merrm);
		dmm(&merrm, pgm, zgrad, "nn", 1);
		memcpy(zgrads->p+istep*aster->tsa*2, zgrad->p, sizeof(double)*aster->tsa*2);
		switch(parms->skyc.servo){
		case 1:
		    dadd(&st2t->mint, 1, merrm, 0.3);
		    break;
		case 2:
		    servo_typeII_filter(st2t, merrm, dtngs, aster->gain->p[idtrat]);
		    break;
		default:
		    error("Invalid\n");
		}
	    }
	}else{
	    if(istep % dtrat == 0){
		for(long iwfs=0; iwfs<aster->nwfs; iwfs++){
		    dcellzero(psf[iwfs]);
		}
	    }
	    for(long iwfs=0; iwfs<aster->nwfs; iwfs++){
		const double thetax=aster->wfs[iwfs].thetax;
		const double thetay=aster->wfs[iwfs].thetay;
		const int ipowfs=aster->wfs[iwfs].ipowfs;
		const long nsa=parms->maos.nsa[ipowfs];
		PCCELL(aster->wfs[iwfs].wvfout[istep],wvfout);
		for(long iwvl=0; iwvl<nwvl; iwvl++){
		    double wvl=parms->maos.wvl[iwvl];
		    for(long isa=0; isa<nsa; isa++){
			ccp(&wvfc->p[iwfs], wvfout[iwvl][isa]);
			//Apply NGS mode error to PSF.
			ngsmod2wvf(wvfc->p[iwfs], wvl, merr, powfs[ipowfs].cloc[isa],
				   powfs[ipowfs].fpc[isa], thetax, thetay, parms);
			cembed(wvf->p[iwfs],wvfc->p[iwfs],0,C_FULL);
			cfft2(wvf->p[iwfs],-1);
			cabs22d(&psf[iwfs]->p[isa+nsa*iwvl], 1., wvf->p[iwfs], 1.);//peak in corner.
		    }//isa
		}//iwvl
	    }//iwfs
	    dcp(&mreal, st2t->mint);
	    if((istep+1) % dtrat == 0){//has output
		phycount++;
		//Form detector image
		double igrad[2];
		int itsa=0;
		for(long iwfs=0; iwfs<aster->nwfs; iwfs++){
		    dcellzero(ints[iwfs]);
		    const int ipowfs=aster->wfs[iwfs].ipowfs;
		    const long nsa=parms->maos.nsa[ipowfs];
		    for(long isa=0; isa<nsa; isa++){
			for(long iwvl=0; iwvl<nwvl; iwvl++){
			    double siglev=aster->wfs[iwfs].siglev[iwvl];
			    ccpd(&otf->p[iwfs],psf[iwfs]->p[isa+nsa*iwvl]);
			    cifft2(otf->p[iwfs], 1); //turn to OTF, peak in corner
			    ccwm(otf->p[iwfs], powfs[ipowfs].dtf[iwvl].nominal);
			    cfft2(otf->p[iwfs], -1);
			    spmulcreal(ints[iwfs]->p[isa]->p, powfs[ipowfs].dtf[iwvl].si, 
				       otf->p[iwfs]->p, siglev);
			}
		
			//Add noise and apply matched filter.
			switch(noisy){
			case 0://no noise at all.
			    break;
			case 1://both poisson and read out noise.
			    addnoise(ints[iwfs]->p[isa], &aster->rand, aster->wfs[iwfs].bkgrnd*dtrat,
				     1, rnefs[ipowfs][idtrat]);
			    break;
			case 2://there is still poisson noise.
			    addnoise(ints[iwfs]->p[isa], &aster->rand, 0, 1, 0);
			    break;
			default:
			    error("Invalid noisy\n");
			}
			for(long ix=0; ix<ints[iwfs]->p[isa]->nx*ints[iwfs]->p[isa]->ny; ix++){
			    if(ints[iwfs]->p[isa]->p[ix]<0){
				ints[iwfs]->p[isa]->p[ix]=0;
			    }
			}
			igrad[0]=0;
			igrad[1]=0;
			if(parms->skyc.mtch){
			    dmulvec(igrad, mtche[iwfs]->p[isa], ints[iwfs]->p[isa]->p, 1);
			}else{
			    double imax=dmax(ints[iwfs]->p[isa]);
			    dcog(igrad, ints[iwfs]->p[isa], 0, 0, 0.1*imax, 0.1*imax); 
			    igrad[0]*=parms->skyc.pixtheta[ipowfs];
			    igrad[1]*=parms->skyc.pixtheta[ipowfs];
			}
			grad->p[isa+itsa]=igrad[0];
			grad->p[isa+nsa+itsa]=igrad[1];
		
		    }//isa
		    itsa+=nsa*2;
		    //if(parms->skyc.dbg){
		    //		dcellwrite(ints[iwfs],"%s/ints/ints_aster%d_iwfs%d_istep%d",
		    //		   dirsetup,aster->iaster,iwfs,istep);
		    //}

		}//iwfs
		dzero(merrm);
		dmm(&merrm, pgm, grad, "nn", 1);
		memcpy(grads->p+istep*aster->tsa*2, grad->p, sizeof(double)*aster->tsa*2);
		switch(parms->skyc.servo){
		case 1:
		    dadd(&st2t->mint, 1, merrm, 0.5);
		    break;
		case 2:
		    servo_typeII_filter(st2t, merrm, dtngs, aster->gain->p[idtrat]);
		    break;
		default:
		    error("Invalid\n");
		}
	    }//if dtrat
	}//if phystart
    }//istep;
    if(parms->skyc.dbg){
	dwrite(zgrads,"%s/skysim_zgrads_aster%d_dtrat%d",dirsetup,aster->iaster,dtrat);
	dwrite(grads,"%s/skysim_grads_aster%d_dtrat%d",dirsetup, aster->iaster,dtrat);
    }
  
    dfree(mreal);
    dfree(merr);
    dfree(merrm);
    dfree(grad);
    dfree(zgrad);
    dfree(grads);//////////////remove this dmat after debugging
    dfree(zgrads);//////////////remove this dmat debugging
    for(long iwfs=0; iwfs<aster->nwfs; iwfs++){
	dcellfree(psf[iwfs]);
	dcellfree(ints[iwfs]);
    }
    ccellfree(wvf);
    ccellfree(wvfc);
    ccellfree(otf);
    dfree(st2t->mint);
    dfree(st2t->mintfirst);
    dfree(st2t->merrlast);
    dfree(st2t->mlead);
    free(st2t);
    free(psf);
    free(mtche);
    free(ints);
    if(parms->skyc.dbg){
	dwrite(mres,"%s/skysim_phy_mres_aster%d_%d",dirsetup,aster->iaster,dtrat);
    }
    //dfree(mres);
    *mresout=mres;
    dscale(res, 1./(aster->nstep-parms->skyc.evlstart));
    return res;
}

/**
   Save NGS WFS and other information for later use in MAOS simulations.*/
void skysim_save(SIM_S *simu, ASTER_S *aster, double *ipres, int selaster, int seldtrat, int isky){
    const PARMS_S* parms=simu->parms;
    const int nwvl=parms->skyc.nwvl;
    char path[PATH_MAX];
    snprintf(path,PATH_MAX,"Res%d_%d_maos/sky%d",simu->seed_maos,simu->seed_skyc,isky);
    mymkdir("%s",path);
    PDCELL(aster[selaster].nea_tot,nea_tot);
    for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
	dcell *sepsf=dcelldup(aster[selaster].wfs[iwfs].pistat->psf);
	for(int ic=0; ic<sepsf->nx*sepsf->ny; ic++){
	    dfftshift(sepsf->p[ic]);//put peak in center. required by MAOS.
	}
	dcellwrite(sepsf, "%s/pistat_wfs%d",path,iwfs+6);
	dcellfree(sepsf);
	dwrite(nea_tot[seldtrat][iwfs], "%s/nea_tot_wfs%d",path,iwfs+6);
	dwrite(aster[selaster].wfs[iwfs].pistat->sanea->p[seldtrat], 
	       "%s/nea_wfs%d",path,iwfs+6);
	dcellwrite(aster[selaster].wfs[iwfs].pistat->sanea, 
		   "%s/neafull_wfs%d",path,iwfs+6);
    }
    dwrite(aster[selaster].gain->p[seldtrat], "%s/gain",path);
    dwrite(simu->mres->p[isky], "%s/mres",path);
    dwrite(simu->psd_tt_ws,"%s/psd_tt_ws",path);
    dwrite(simu->psd_ps,"%s/psd_ps",path);
    char fnconf[PATH_MAX];
    snprintf(fnconf,PATH_MAX,"%s/base.conf",path);
    FILE *fp=fopen(fnconf,"w");

    fprintf(fp,"wfs.thetax=[0 0  -33.287 -20.5725  20.5725 33.287");
    for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
	fprintf(fp," %.4f", aster[selaster].wfs[iwfs].thetax*206265);
    }
    fprintf(fp,"]\n");
    fprintf(fp,"wfs.thetay=[0 35 10.8156 -28.3156 -28.3156 10.8156");
    for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
	fprintf(fp," %.4f", aster[selaster].wfs[iwfs].thetay*206265);
    }
    fprintf(fp,"]\n");

    fprintf(fp,"wfs.siglev=[900 900 900 900 900 900");
    for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
	fprintf(fp, " %.2f", aster[selaster].wfs[iwfs].siglevtot);
    }
    fprintf(fp,"]\n");
    fprintf(fp,"wfs.wvlwts=[1 1 1 1 1 1");
    for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
	for(int iwvl=0; iwvl<nwvl; iwvl++){
	    fprintf(fp," %.2f ", aster[selaster].wfs[iwfs].siglev[iwvl]
		    /aster[selaster].wfs[iwfs].siglevtot);
	}
    }
    fprintf(fp,"]\n");
    fprintf(fp,"wfs.powfs=[0 0 0 0 0 0");
    for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
	fprintf(fp, " %d", aster[selaster].wfs[iwfs].ipowfs+1);
    }
    fprintf(fp,"]\n");
    fprintf(fp,"sim.seeds=[%d]\n",simu->seed_maos);
    fprintf(fp,"sim.end=5000\n");
    fprintf(fp,"sim.dt=%g\n", parms->maos.dt);
    fprintf(fp,"atm.r0z=%.4f\n", parms->maos.r0z);
    fprintf(fp,"sim.zadeg=%g\n", parms->maos.zadeg);
    fprintf(fp,"atm.size=[128 128]\n");
    fprintf(fp,"tomo.split_wt=3\n");
    PDMAT(parms->skyc.rnefs,rnefs);
    double rne=rnefs[0][seldtrat];
    double bkgrnd=aster[selaster].wfs[0].bkgrnd;

    if(parms->maos.npowfs==2){
	fprintf(fp, "powfs.piinfile=[\"\" \"pistat\" \"pistat\"]\n");
	fprintf(fp, "powfs.neareconfile=[\"\" \"nea_tot\" \"nea_tot\"]\n");
	fprintf(fp, "powfs.phyusenea=[0 1 1]\n");
	fprintf(fp, "powfs.dtrat=[1 %d %d]\n", parms->skyc.dtrats[seldtrat],
		parms->skyc.dtrats[seldtrat]);
	fprintf(fp, "powfs.bkgrnd=[0 %.2f %.2f]\n", bkgrnd, bkgrnd);
	fprintf(fp, "powfs.rne=[3 %.2f %.2f]\n", rne,rne);
	fprintf(fp, "powfs.phystep=[0 %ld %ld]\n", 
		(long)50+parms->skyc.dtrats[seldtrat]*20, 
		(long)50+parms->skyc.dtrats[seldtrat]*20);
	fprintf(fp, "powfs.noisy=[1 1 1]\n");
	fprintf(fp, "powfs.pixtheta=[0.5/206265 %g/206265000 %g/206265000]\n",
		parms->skyc.pixtheta[0]*206265000,
		parms->skyc.pixtheta[1]*206265000);
	fprintf(fp, "powfs.pixpsa=[6 %d %d]\n",
		parms->skyc.pixpsa[0], 
		parms->skyc.pixpsa[1]);
	fprintf(fp, "powfs.ncomp=[64 %d %d]\n", 
		parms->maos.ncomp[0], parms->maos.ncomp[1]);
	fprintf(fp, "powfs.nwvl=[1 %d %d]\n",nwvl,nwvl);
	fprintf(fp, "powfs.wvl=[0.589e-6");
	for(int ip=0; ip<2; ip++){
	    for(int iwvl=0; iwvl<nwvl; iwvl++){
		fprintf(fp, " %.4g", parms->maos.wvl[iwvl]);
	    }
	}
	fprintf(fp,"]\n");
    }else{
	error("Fill this out please\n");
    }
    fprintf(fp,"sim.servotype_lo=2\n");//type II
    fprintf(fp, "sim.gtypeII_lo=\"gain.bin.gz\"\n");
    if(parms->maos.wddeg){
	fprintf(fp, "atm.wddeg=[");
	for(int ips=0; ips<parms->maos.nwddeg; ips++){
	    fprintf(fp, "%.2f ", parms->maos.wddeg[ips]);
	}
	fprintf(fp, "]\n");
    }
    fclose(fp);
    snprintf(fnconf,PATH_MAX,"%s/skyres.txt",path);
    fp=fopen(fnconf,"w");
    fprintf(fp, "TotAll\tNGS\tTT\n");
    fprintf(fp, "%g\t%g\t%g\n",
	    sqrt(ipres[0])*1e9, sqrt(ipres[1])*1e9, sqrt(ipres[2])*1e9);
    fclose(fp);
}
