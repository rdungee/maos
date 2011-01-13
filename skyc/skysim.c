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

/**
   \file skysim.c

   Contains the skysim() function that does the skycoverage simulation.
   
   The steps:
   *) Read in parameters
   *) Read in necessary struct saved by MAOS and/or prepared by user
   *) Read in star catalog.
   *) Read in Focus tracking error as a function of sampling frequency
   *) Read in NGS Mode time vector
   *) For each possible asterism in a star field.
   #) Compute NGS Gradient operator from NGS modes.
   #) Read in saved PSF and interpolate to the location of NGS. apply extra strehl degradation.
   #) Compute NGS signal level and measurement error from matched filter (be careful about PSF quality).
   #) For each sample frequency.
   $) Compute the reconstructor
   $) Optimize for the loop gain.
   $) Estimate residual windshake error.
   $) Do a rough open loop estimate of the total error. (atm + windshake +focus tracking ...)
   $) Guess for the best sampling frequency (compare with closed loop results to validate this model).
   *) Compare the asterisms and pick only a few that are within ~200% of the best asterism at their best sampling freq.
   *) For each remaining asterism in the star field.
   #) Do a temporal simulation to obtain the final result.
   #) save error due to focus tracking, telescope windshake, separately.
   #) Find the best sampling frequency.
   *) Find the best asterism.
   *) Repeat for another star field.


   Don't forget: 

   *) Sampling frequency dependent rne. input a precomputed matrix for all
   possible sampling frequencies.
   *) Implementation error.
   *) Windshake (handle separately), focus tracking jitter
   *) On axis performance.

   2010-06-16: Single asterism debugged.
   2010-06-21: Bug found: wfs.thetax, thetay not copied from star.
   differences:
   MATLAB starts with a guess using QC model for NEA when regenerating pistat
   maos starts with loaded pistat. This is 1 difference.
   The gains are not sensitive to which PSD is loaded, or whether PSD is scaled or not.
*/

#include "skyc.h"
#include "parms.h"
#include "skysim.h"
#include "types.h"
#include "setup_powfs.h"
#include "setup_aster.h"
#include "skysim_utils.h"
#include "mtch.h"
#include "genstars.h"
#include "genpistat.h"
#include "setup_star.h"
#include "servo.h"
#include "utils.h"
#define VERBOSE 2

/**
   Compute Open loop NGS mode wavefront error from mode vectors.  */
static dmat* calc_rmsol(dmat *mideal, const PARMS_S *parms){
    double rms=0, rmstt=0;
    PDMAT(mideal, pmideal);
    for(long istep=0; istep<mideal->ny; istep++){
	rms+=dwdot(pmideal[istep], parms->maos.mcc, pmideal[istep]);
	rmstt+=dwdot2(pmideal[istep],parms->maos.mcc_tt,pmideal[istep]);
    }
    dmat *rmsol=dnew(2,1);
    rmsol->p[0]=rms/mideal->ny;
    rmsol->p[1]=rmstt/mideal->ny;
    return rmsol;
}

/**
   The actual work horse that does the physical optics time domain simulation.
   */
static void skysim_isky(SIM_S *simu){
    int isky;
    POWFS_S *powfs=simu->powfs;
    const PARMS_S *parms=simu->parms;
    const dcell *stars=simu->stars;
    const int noisy=parms->skyc.noisy;
    const int do_demote_end=parms->skyc.demote?2:1;
    int nstar;
    int nstep;
    const int seed_skyc=simu->seed_skyc;
    const int seed_maos=simu->seed_maos;

    PDMAT(simu->res, pres);
    PDMAT(simu->res_oa, pres_oa);
    PDMAT(simu->res_ol, pres_ol);
    while(LOCK(simu->mutex_isky),isky=simu->isky++,UNLOCK(simu->mutex_isky),isky<simu->isky_end){
	double tk_1=myclockd();
	//Setup star parameters.
	STAR_S *star=setup_star(&nstar, simu,stars->p[isky],seed_maos);
	if(!star){
	    info2("Field %d, No stars available\n", isky);
	    continue;
	}
	int naster;
	ASTER_S *aster=setup_aster_comb(&naster, nstar, parms);
	if(!aster || aster==0){
	    info2("Field %d, Aster is empty. skip\n", isky);
	    continue;
	}
	TIC;tic;
	for(int iaster=0; iaster<naster; iaster++){
	    //Parallelizing over aster gives same random stream.
	    seed_rand(&aster[iaster].rand, seed_skyc+iaster+40);
	    //Compute signal level.
	    setup_aster_copystar(&aster[iaster], star, parms);
	    //setup gradient operator.
	    setup_aster_g(&aster[iaster], star, powfs, parms);
	    //Compute the reconstructor, nea, sigman
	    setup_aster_recon(&aster[iaster], star,parms);
	    //Optimize servo gains.
	    setup_aster_servo(simu, &aster[iaster], parms);
	}
	if(parms->skyc.verbose){
	    toc2("Estimation");
	}
	//Select asters that have good performance.
	setup_aster_select(pres_ol[isky],aster, naster, star, 
			   simu->rmsol->p[0]/9,parms); 
	//Read in physical optics data (wvf)
	nstep=setup_star_read_wvf(star,nstar,parms,seed_maos);
	//Physical Optics Simulations.

	double skymini=INFINITY;
	int selaster=0;
	int seldtrat=0;
	for(int iaster=0; iaster<naster; iaster++){
	    ASTER_S *asteri=&aster[iaster];
	    asteri->nstep=nstep;
	    if(!asteri->use || (parms->skyc.dbgaster>-1 && iaster!=parms->skyc.dbgaster)){
		continue;
	    }
	    if(parms->skyc.verbose>1){
		for(int iwfs=0; iwfs<aster[iaster].nwfs; iwfs++){
		    info2("wfs %d: istar=%d, ipowfs=%d\n",iwfs,aster[iaster].wfs[iwfs].istar,
			  aster[iaster].wfs[iwfs].ipowfs);
		    info2("wfs %d: at (%g,%g). siglev=%g\n",iwfs,
			  aster[iaster].wfs[iwfs].thetax*206265,
			  aster[iaster].wfs[iwfs].thetay*206265, aster[iaster].wfs[iwfs].siglevtot);
		}
	    }
	    if(parms->skyc.verbose){
		info2("Aster %d, Estimated minimum error is %.2fnm at %.1f Hz\n", iaster,
		      sqrt(asteri->mresol)*1e9, parms->skyc.fss[asteri->mdtrat]);
	    }
	    //Copy wvf from star to aster
	    setup_aster_wvf(asteri, star, parms);
	    //Regenerate PSD for this combination only.
	    setup_aster_regenpsf(simu->mideal, asteri,powfs,parms);
	    //Redo matched filter.
	    setup_aster_redomtch(asteri, powfs, parms);
	    //Compute the reconstructor, nea, sigman
	    setup_aster_recon(asteri, star, parms);
	    //Optimize servo gains.
	    setup_aster_servo(simu, asteri, parms);

	    double mini=INFINITY;
	    dmat *pmini=NULL;
	    dmat *min_imres=NULL;
	    int mdtrat=0;
	    int demote=-1;
	    for(int idtrat=aster->idtratmin; idtrat<=aster->idtratmax; idtrat++){
		//focus and windshake residual;
		double resadd=asteri->res_ws->p[idtrat] + parms->skyc.resfocus->p[idtrat];
		dmat *ires=NULL;
		dmat *imres=NULL;
		for(int do_demote=0; do_demote<do_demote_end; do_demote++){
		    ires=skysim_phy(&imres,simu, asteri, powfs, parms, idtrat, noisy, do_demote);
		    if(ires){
			if(parms->skyc.verbose){
			    if(do_demote){
				info2("  DemoteTTF: ");
			    }
			    info2("%5.1f Hz %7.2f +%7.2f =%7.2f", parms->skyc.fss[idtrat], 
				  sqrt(ires->p[0])*1e9, sqrt(resadd)*1e9,
				  sqrt(ires->p[0]+resadd)*1e9);
			}
			//Add windshake contribution.
			double tot_1=ires->p[0] + resadd;
			if(tot_1 < mini){
			    mini=tot_1;
			    mdtrat=idtrat;
			    demote=do_demote;
			    dfree(pmini); pmini=dref(ires);
			    dfree(min_imres); min_imres=dref(imres);
			}
		    }
		}
		if(ires && parms->skyc.verbose){
		    info2("\n");
		}
		dfree(ires);
		dfree(imres);
	    }
	    if(mini<skymini){
		selaster=iaster;
		skymini=mini;
		//Field Averaged Performance.
		pres[isky][1]=pmini->p[0];
		pres[isky][2]=pmini->p[1];
		pres[isky][3]=asteri->res_ws->p[mdtrat];
		pres[isky][4]=parms->skyc.resfocus->p[mdtrat];
		pres[isky][0]=pres[isky][1]+pres[isky][3]+pres[isky][4];
		//On axis performance.
		pres_oa[isky][1]=pmini->p[2];
		pres_oa[isky][2]=pmini->p[3];
		pres_oa[isky][3]=asteri->res_ws->p[mdtrat];
		pres_oa[isky][4]=parms->skyc.resfocus->p[mdtrat];
		pres_oa[isky][0]=pres_oa[isky][1]+pres_oa[isky][3]+pres_oa[isky][4];
		seldtrat = mdtrat;
		simu->fss->p[isky]=parms->skyc.fss[mdtrat];
		simu->demote->p[isky]=demote;
		if(parms->skyc.verbose){
		    info2("%5.1f Hz: Update Tot: %6.2f nm NGS: %6.2f nm TT: %6.2f nm\n", 
			  simu->fss->p[isky],
			  sqrt(pres[isky][0])*1e9, sqrt(pres[isky][1])*1e9, sqrt(pres[isky][2])*1e9);
		}
		dcp(&simu->mres->p[isky], min_imres);
	    }
	    dfree(pmini);
	    dfree(min_imres);
	}//iaster

	PDMAT(simu->sel->p[isky],psel);
	for(int iwfs=0; iwfs<aster[selaster].nwfs; iwfs++){
	    psel[iwfs][0]=aster[selaster].wfs[iwfs].thetax;
	    psel[iwfs][1]=aster[selaster].wfs[iwfs].thetay;
	    for(int iwvl=0; iwvl<parms->maos.nwvl; iwvl++){
		psel[iwfs][iwvl+2]=aster[selaster].wfs[iwfs].mags->p[iwvl];
	    }
	}
	if(parms->skyc.save){
	    skysim_save(simu, aster, pres[isky], selaster, seldtrat, isky);
	}
	free_aster(aster, naster, parms);
	free_star(star, nstar, parms);
	double tk_2=myclockd();
	LOCK(simu->mutex_save);
	simu->status->isim=isky;
	simu->status->tot=tk_2-tk_1;//per step
	simu->status->laps=tk_2-simu->tk_0;
	int nsky_seed=simu->isky_end-simu->isky_start;
	int nseed_tot=parms->maos.nseed*parms->skyc.nseed;
	int nsky_left=simu->isky_end-simu->isky-1+nsky_seed*(nseed_tot-simu->iseed-1);
	int nsky_laps=simu->isky-simu->isky_start+1+nsky_seed*simu->iseed;
	simu->status->rest=simu->status->laps*nsky_left/nsky_laps;
	simu->status->clerrlo=sqrt(pres[isky][1])*1e9;
	simu->status->clerrhi=sqrt(pres[isky][0])*1e9;
	scheduler_report(simu->status);
	UNLOCK(simu->mutex_save);
	long totm=(long)floor(simu->status->tot/60.);
	long tots=(long)simu->status->tot-totm*60;
	long laps_h=simu->status->laps/3600;
	long laps_m=simu->status->laps/60-laps_h*60;
	long rest_h=simu->status->rest/3600;
	long rest_m=simu->status->rest/60-rest_h*60;
	info2("Field %d, %2d stars, %5.1f Hz: Final Tot: %6.2f nm NGS: %6.2f nm TT: %6.2f nm"
	      "Tot %ld:%2ld. Used %ld:%2ld, Left %ld:%2ld\n",
	      isky, nstar, simu->fss->p[isky],
	      sqrt(pres[isky][0])*1e9, sqrt(pres[isky][1])*1e9, sqrt(pres[isky][2])*1e9,
	      totm, tots, laps_h, laps_m, rest_h, rest_m);
    }//while
}
/**
   Setup the stars fields and then calls skysim_isky() to handle each star field.
*/
void skysim(const PARMS_S *parms){
    if(parms->skyc.dbg){
	dwrite(parms->skyc.resfocus, "%s/resfocus",dirsetup);
    }
    const int npowfs=parms->maos.npowfs;
    SIM_S *simu=calloc(1, sizeof(SIM_S));
    simu->status=calloc(1, sizeof(STATUS_T));
    simu->status->info=S_RUNNING;
    simu->status->scale=1;
    simu->status->nseed=parms->maos.nseed*parms->skyc.nseed;
    simu->status->nthread=parms->skyc.nthread;
    simu->status->timstart=myclocki();
    simu->powfs=calloc(npowfs, sizeof(POWFS_S));
    simu->parms=parms;
    setup_powfs(simu->powfs, parms);
    genpistat(parms, simu->powfs); 
    PINIT(simu->mutex_isky);
    PINIT(simu->mutex_save);
    simu->tk_0=myclockd();
    for(int iseed_maos=0; iseed_maos<parms->maos.nseed; iseed_maos++){
	simu->seed_maos=parms->maos.seeds[iseed_maos];//loop over seed
	//Read ideal NGS mode.
	simu->mideal=dread("%s_%d.bin",parms->maos.fnmideal,simu->seed_maos);
	//Read ideal NGS mode dot product for on axis.
	dcell *midealp=dcellread("%s_%d.bin",parms->maos.fnmidealp,simu->seed_maos);
	simu->mideal_oa=dref(midealp->p[parms->maos.evlindoa]);
	dcellfree(midealp);
	simu->rmsol=calc_rmsol(simu->mideal, parms);
	simu->psd_ws=ddup(parms->skyc.psd_ws);
	simu->psd_ngs=ddup(parms->skyc.psd_ngs);
	simu->psd_ps=ddup(parms->skyc.psd_ps);
	simu->psd_tt=ddup(parms->skyc.psd_tt);
	prep_bspstrehl(simu);
	//renormalize PSD
	if(parms->skyc.psd_scale){
	    double rms_ngs=psd_intelog2(simu->psd_ngs);
	    info("PSD integrates to %.2f nm\n", sqrt(rms_ngs)*1e9);
	    double rms_ratio=simu->rmsol->p[0]/rms_ngs;
	    info("Scaling PSD by %g\n", rms_ratio);
	    long nx=simu->psd_ngs->nx;
	    //scale PSF in place.
	    double *p_ngs=simu->psd_ngs->p+nx;
	    double *p_tt=simu->psd_tt->p+nx;
	    double *p_ps=simu->psd_ps->p+nx;
	    for(long i=0; i<nx; i++){
		p_ngs[i]*=rms_ratio;
		p_tt[i]*=rms_ratio;
		p_ps[i]*=rms_ratio;
	    }
	    rms_ngs=psd_intelog2(simu->psd_ngs);
	    info("PSD integrates to %.2f nm after scaling\n", sqrt(rms_ngs)*1e9);
	}
	simu->psd_ngs_ws=add_psd(simu->psd_ngs, simu->psd_ws);
	simu->psd_tt_ws=add_psd(simu->psd_tt, simu->psd_ws);
	for(int iseed_skyc=0; iseed_skyc<parms->skyc.nseed; iseed_skyc++){
	    simu->status->iseed=iseed_skyc+iseed_maos*parms->skyc.nseed;
	    simu->seed_skyc=parms->skyc.seeds[iseed_skyc];
	    seed_rand(&simu->rand, simu->seed_skyc+parms->maos.zadeg);
	  

	    info2("Open loop error: NGS: %.2f TT: %.2f nm\n", 
		  sqrt(simu->rmsol->p[0])*1e9, sqrt(simu->rmsol->p[1])*1e9);
	   

	    //generate star fields.
	    if(parms->skyc.stars){
		info2("Loading stars from %s\n",parms->skyc.stars);
		simu->stars=dcellread("%s",parms->skyc.stars);
		sortstars(simu->stars);
	    }else{
		simu->stars=genstars(parms->skyc.nsky, 
				     parms->skyc.lat, parms->skyc.lon,parms->skyc.catscl,
				     parms->skyc.patfov,parms->maos.nwvl, 
				     parms->maos.wvl, &simu->rand);
	    }
	    dcellwrite(simu->stars, "Res%d_%d_stars",simu->seed_maos,simu->seed_skyc);
	    int nsky=simu->stars->nx*simu->stars->ny;
	    if(nsky > parms->skyc.nsky){
		warning("nsky is reduced from %d to %d\n", 
			nsky, parms->skyc.nsky);
		nsky=parms->skyc.nsky;
	    }
	    int seed_maos=simu->seed_maos;
	    int seed_skyc=simu->seed_skyc;
	    simu->res   =dnew_mmap(5,nsky,NULL, "Res%d_%d", seed_maos, seed_skyc);
	    simu->res_oa=dnew_mmap(5,nsky,NULL, "Res%d_%d_oa", seed_maos, seed_skyc);
	    simu->res_ol=dnew_mmap(3,nsky,NULL, "Res%d_%d_ol", seed_maos, seed_skyc);
	    simu->fss   =dnew_mmap(nsky,1,NULL, "Res%d_%d_fss", seed_maos, seed_skyc);
	    simu->demote=dnew_mmap(nsky,1,NULL, "Res%d_%d_demote", seed_maos, seed_skyc);
	    simu->sel   =dcellnewsame_mmap(nsky, 1, 2+parms->maos.nwvl, parms->skyc.nwfstot,
					   NULL,"Res%d_%d_sel", seed_maos, seed_skyc);
	    simu->mres  =dcellnewsame_mmap(nsky, 1, parms->maos.nmod, parms->maos.nstep,
					   NULL,"Res%d_%d_mres",  seed_maos, seed_skyc);
	    dset(simu->res, simu->rmsol->p[0]);
	    dset(simu->res_oa, simu->rmsol->p[0]);
	    simu->isky_start=parms->skyc.start;
	    simu->isky_end=nsky;
	    if(parms->skyc.dbgsky>-1){
		simu->isky_start=parms->skyc.dbgsky;
		simu->isky_end=parms->skyc.dbgsky+1;
	    }
	    simu->status->simstart=simu->isky_start;
	    simu->status->simend=simu->isky_end;
	    if(simu->isky_start >= simu->isky_end){
		continue;//nothing to do.
	    }
	    simu->isky=simu->isky_start;
	    CALL(skysim_isky, simu, parms->skyc.nthread);//isky iteration.
	    dcellfree(simu->stars);
	    dfree(simu->res);
	    dfree(simu->res_oa);
	    dfree(simu->res_ol);
	    dcellfree(simu->mres);
	    dcellfree(simu->sel);
	    dfree(simu->fss);
	    simu->iseed++;
	}//iseed_skyc
	dfree(simu->mideal);
	dfree(simu->mideal_oa);
	dfree(simu->rmsol);
	dfree(simu->psd_ws);
	dfree(simu->psd_ngs);
	dfree(simu->psd_ps);
	dfree(simu->psd_ws);
	dfree(simu->psd_tt_ws);
	dfree(simu->psd_ngs_ws);
	//Free the data used to do bicubic spline.
	for(int ipowfs=0; ipowfs<parms->maos.npowfs; ipowfs++){
	    for(int ic=0; ic<parms->maos.nsa[ipowfs]*parms->maos.nwvl; ic++){
		dcellfree(simu->bspstrehl[ipowfs][ic]);
	    }
	    free(simu->bspstrehl[ipowfs]);
	}
	free(simu->bspstrehl);
	dfree(simu->bspstrehlxy);
    }//iseed_maos
    free(simu->status);
    free_powfs(simu->powfs,parms);
    free(simu->powfs);
    free(simu);
}
