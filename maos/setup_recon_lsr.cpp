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
#include "ahst.h"

/**
   Setup the least square reconstruct by directly inverting GA matrix. 
   The reconstructor is simply the pseudo inverse of GA matrix:
   \f[\hat{x}=(G_a^TC_g^{-1}G_a)^{-1}G_a^TC_g^{-1}\f]

   This is very close to RR except replacing GX with GA.

   We use the tomograhy parameters for lsr, since lsr is simply "tomography" onto DM directly.
*/
void setup_recon_lsr(RECON_T *recon, const PARMS_T *parms){
    const int ndm=parms->ndm;
    const int nwfs=parms->nwfsr;
    cell *GAlsr;
    cell *GAM=parms->recon.modal?(cell*)recon->GM:(cell*)recon->GA;
    if(parms->recon.split){ //high order wfs only in split mode. 
	GAlsr=parms->recon.modal?(cell*)recon->GMhi:(cell*)recon->GAhi;
    }else{ //all wfs in integrated mode. 
	GAlsr=GAM;
    }
    int free_GAlsr=0;
    if(GAlsr->p[0]->id!=M_DBL){
	dsp *tmp=dsp_cast(GAlsr->p[0]);
	if(tmp->nzmax>tmp->nx*tmp->ny*0.2){//not very sparse
	    dcell *tmp2=0;
	    free_GAlsr=1;
	    dcelladd(&tmp2, 1, (dspcell*)GAlsr, 1);
	    GAlsr=(cell*)tmp2;
	}
    }
    info2("Building recon->LR\n");
    recon->LR.M=dcellmm2(GAlsr, recon->saneai, "tn");
    // Tip/tilt and diff focus removal low rand terms for LGS WFS.
    if(recon->TTF){
	dcellmm(&recon->LR.U, recon->LR.M, recon->TTF, "nn", 1);
	recon->LR.V=dcelltrans(recon->PTTF);
    }
    info2("Building recon->LL\n");
    recon->LL.M=dcellmm2(recon->LR.M, GAlsr, "nn");
    if(free_GAlsr){
	cellfree(GAlsr);
    }
    double maxeig=pow(recon->neamhi * recon->aloc->p[0]->dx, -2);
    if(parms->recon.modal){
	double strength=1;
	for(int idm=0; idm<ndm; idm++){
	    strength*=dnorm(recon->amod->p[idm]);
	}
	strength=pow(strength, 2./ndm);
	maxeig*=strength;
    }
    if(fabs(parms->lsr.tikcr)>EPS){
	info2("Adding tikhonov constraint of %g to LLM\n", parms->lsr.tikcr);
	info2("The maximum eigen value is estimated to be around %g\n", maxeig);
	dcelladdI(recon->LL.M, parms->lsr.tikcr*maxeig);
    }
    dcell *NW=NULL;
    if(!parms->recon.modal){
	if(parms->lsr.alg!=2){
	    /* Not SVD, need low rank terms for piston/waffle mode constraint. */
	    NW=dcellnew(ndm,1);
	    int nmod=2;/*two modes. */
	    for(int idm=0; idm<ndm; idm++){
		loc_create_map(recon->aloc->p[idm]);
		const long nloc=recon->aloc->p[idm]->nloc;
		NW->p[idm]=dnew(nloc, ndm*nmod);
		double *p=NW->p[idm]->p+nmod*idm*nloc;
		const double *cpl=recon->actcpl->p[idm]->p;
		for(long iloc=0; iloc<nloc; iloc++){
		    if(cpl[iloc]>0.1){
			p[iloc]=1;/*piston mode */
		    }
		}
		/*notice offset of 1 because map start count at 1 */
		p=NW->p[idm]->p+(1+nmod*idm)*nloc-1;
		map_t *map=recon->aloc->p[idm]->map;
		for(long iy=0; iy<map->ny; iy++){
		    for(long ix=0; ix<map->nx; ix++){
			if(IND(map,ix,iy)){
			    p[(long)IND(map,ix,iy)]=(double)2*((iy+ix)&1)-1;
			}
		    }
		}
	    }
	    /*scale it to match the magnitude of LL.M */
	    dcellscale(NW, sqrt(maxeig));
	    if(parms->save.setup){
		writebin(NW, "lsrNW");
	    }
	}
	if(parms->lsr.actslave){
	    /*actuator slaving. important. change from 0.5 to 0.1 on 2011-07-14. */
	    dspcell *actslave=slaving(recon->aloc, recon->actcpl, NW,
				      recon->actstuck, recon->actfloat, parms->lsr.actthres, maxeig);
	    if(parms->save.setup){
		if(NW){
		    writebin(NW, "lsrNW2");
		}
		writebin(actslave,"actslave");
	    }
	    dcelladd(&recon->LL.M, 1, actslave, 1);
	    cellfree(actslave);
	}
    }
    /*Low rank terms for low order wfs. Only in Integrated tomography. */
    dcell *ULo=dcellnew(ndm,nwfs);
    dcell *VLo=dcellnew(ndm,nwfs);
    dcell*  pULo=ULo/*PDELL*/;
    dcell*  pVLo=VLo/*PDELL*/;
    for(int iwfs=0; iwfs<nwfs; iwfs++){
	int ipowfs=parms->wfsr[iwfs].powfs;
	if(parms->powfs[ipowfs].skip || !parms->powfs[ipowfs].lo){
	    continue;
	}
	for(int idm=0; idm<ndm; idm++){
	    dspfull(PIND(pULo,idm,iwfs), (dsp*)IND(recon->LR.M, idm, iwfs),'n',-1);
	    dspfull(PIND(pVLo,idm,iwfs), (dsp*)IND(GAM, iwfs, idm),'t',1);
	}
    }
    recon->LL.U=dcellcat(recon->LR.U, ULo, 2);
    dcell *GPTTDF=NULL;
    dcellmm(&GPTTDF, GAM, recon->LR.V, "tn", 1);
    recon->LL.V=dcellcat(GPTTDF, VLo, 2);
    dcellfree(GPTTDF);
    dcellfree(ULo);
    dcellfree(VLo);
    if(!parms->recon.modal && NW){
	info2("Create piston and check board modes that are in NULL space of GA.\n");
	/*add to low rank terms. */
	dcell *tmp=recon->LL.U;
	recon->LL.U=dcellcat(tmp, NW, 2);
	dcellfree(tmp);
	dcellscale(NW, -1);
	tmp=recon->LL.V;
	recon->LL.V=dcellcat(tmp, NW, 2);
	dcellfree(tmp);
	dcellfree(NW);
    }
    if(parms->lsr.fnreg){
	warning("Loading LSR regularization from file %s.\n", parms->lsr.fnreg);
	dspcell *tmp=dspcellread("%s", parms->lsr.fnreg);
	dcelladd(&recon->LL.M, 1, tmp, 1);
	dspcellfree(tmp);
    }
    recon->LL.alg = parms->lsr.alg;
    recon->LL.bgs = parms->lsr.bgs;
    recon->LL.warm = parms->recon.warm_restart;
    recon->LL.maxit = parms->lsr.maxit;
    /*Remove empty cells. */
    dcelldropempty(&recon->LR.U,2);
    dcelldropempty(&recon->LR.V,2);
    dcelldropempty(&recon->LL.U,2);
    dcelldropempty(&recon->LL.V,2);
    if(parms->save.recon){
	writebin(recon->LR.M,"LRM");
	writebin(recon->LR.U,"LRU");
	writebin(recon->LR.V,"LRV");
	writebin(recon->LL.M,"LLM.bin");/*disable compression */
	writebin(recon->LL.U,"LLU");
	writebin(recon->LL.V,"LLV"); 
    }
    if(parms->lsr.alg==0 || parms->lsr.alg==2){
	if(!parms->lsr.bgs){
	    muv_direct_prep(&recon->LL, (parms->lsr.alg==2)*parms->lsr.svdthres);
	    if(parms->save.recon){
		if(recon->LL.C)
		    chol_save(recon->LL.C, "LLC.bin");
		else
		    writebin(recon->LL.MI, "LLMI.bin");
	    }
	    cellfree(recon->LL.M);
	    dcellfree(recon->LL.U);
	    dcellfree(recon->LL.V);	
	}else{
	    muv_direct_diag_prep(&(recon->LL), (parms->lsr.alg==2)*parms->lsr.svdthres);
	    if(parms->save.recon){
		for(int ib=0; ib<recon->LL.nb; ib++){
		    if(recon->LL.CB)
			chol_save(recon->LL.CB[ib],"LLCB_%d.bin", ib);
		    else
			writebin(recon->LL.MI,"LLMIB_%d.bin", ib);
		}
	    }
	    /*Don't free M, U, V */
	}
    }
}
