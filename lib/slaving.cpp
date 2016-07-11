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
#include "../math/mathdef.h"
#include "slaving.h"

/**
   Compute the actuator coupling coefficient to be used to identify non-coupled
   actuators. W1 is optional weighting function. The max value is 1.
*/
dcell *genactcpl(const dspcell *HA, const dmat *W1){
    int ndm=HA->ny;
    dcell *actcplc=dcellnew(ndm, 1);
    for(int idm=0; idm<ndm; idm++){
	for(int ifit=0; ifit<HA->nx; ifit++){
	    dsp *ha=HA->p[ifit+idm*HA->nx];
	    if(!ha) continue;
	    if(W1){
		dspmm(&actcplc->p[idm], ha, W1, "tn", 1);
	    }else{
		dmat *tmp=dspsumabs(ha, 1);
		tmp->nx=tmp->ny; tmp->ny=1;
		dadd(&actcplc->p[idm], 1, tmp, 1);
		dfree(tmp);
	    }
	}
	/* If a point is fully coupled to any direction, it is treated as fully
	   coupled. */
	double thres=dmax(actcplc->p[idm])/HA->nx;
	double scale=1./thres;
	double *p=actcplc->p[idm]->p;
	for(long i=0; i<actcplc->p[idm]->nx; i++){
	    if(p[i]>thres){
		p[i]=1;
	    }else{
		p[i]*=scale;
	    }
	}
    }
    return actcplc;
}
/**
   Compute slaving actuator regularization. HA (or use GA) is used to compute
   active actuators. If NW is non NULL, orthogonalize it with the slaving
   regularization.  When the actuators are in the NULL space of HA, we want to
   contraint their values to be close to the ones that are active. We put an
   additional term in the fitting matrix to force this. Be careful with it when
   using tip/tilt constraint and cholesky back solve.  */
dspcell *slaving(loccell *aloc,  /**<[in]The actuator grid*/
		 const dcell *actcplc,/**<[in]Actuator coupling coefficiency*/
		 dcell *NW,     /**<[in]The low rank terms that need to be orthogonal to the output (optional)*/
		 const lcell *actstuck,/**<[in]mask for stuck actuators that will not be slaved, but have value constrained.*/
		 const lcell *actfloat,/**<[in]mask for float actuators that will be slaved, but not have value constrained.*/
		 const double thres,  /**<[in]The threshold that an actuator is deemed slave*/
		 const double sclsq   /**<[in]The square of scaling of the overall strength*/
    ){
    if(!actcplc && !actfloat) {
	error("Both actcplc and actfloat are not supplied\n");
    }
    double scl=sqrt(sclsq);
    if(scl<EPS){
	error("scl=%g is too small\n", scl);
    }
    int ndm=aloc->nx;
    dspcell *actslavec=(dspcell*)cellnew(ndm, ndm);/*block diagonal. */
    dspcell*  actslave=actslavec;
    int nslavetot=0;
    /*Next process stuck and floating actuators. Adjust actcplc and compute slaving matrix.*/
    for(int idm=0; idm<ndm; idm++){
	int nact=aloc->p[idm]->nloc;
	int nslave=0;
  	double *actcpl= actcplc->p[idm]->p;
	double *actcpl0 = actcpl-1;
	for(int iact=0; iact<nact; iact++){
	    if(actstuck && actstuck->p[idm]->p[iact]){
		actcpl[iact] = 1;/*always Skip the stuck actuators. */
		nslave++;
	    }
	    if(actfloat && actfloat->p[idm]->p[iact]){
		actcpl[iact] = 0;/*always include the float actuators */
	    }
	    if(actcpl[iact]<thres){
		nslave++;
	    }
	}
	const long *stuck=actstuck?(actstuck->p[idm]?actstuck->p[idm]->p():0):0;
	const long *floated=actfloat?(actfloat->p[idm]?actfloat->p[idm]->p():0):0;

	nslavetot+=nslave;
	info2("dm %d: there are %d slave actuators\n", idm, nslave);
	if(nslave==0) {
	    continue;
	}
	loc_create_map(aloc->p[idm]);
	map_t *map=aloc->p[idm]->map;
	const double ox=map->ox;
	const double oy=map->oy;
	const double dx1=1./aloc->p[idm]->dx;
	const double dy1=1./aloc->p[idm]->dy;
	//dsp *slavet=dspnew(nact,nslave,nslave*5);
	dsp *slavet=dspnew(nact,nact,nslave*5);
	spint *pp=slavet->p;
	spint *pi=slavet->i;
	double *px=slavet->x;
	const double *locx=aloc->p[idm]->locx;
	const double *locy=aloc->p[idm]->locy;
	long count=0;
	for(int iact=0; iact<nact; iact++){
	    pp[iact]=count;
	    if(stuck && stuck[iact]){/*limit the strength of stuck actuators. */
		pi[count]=iact;
		px[count]=scl;
		count++;
	    }else if(actcpl[iact]<thres){/*slave actuators */
		long mapx=(long)round((locx[iact]-ox)*dx1);
		long mapy=(long)round((locy[iact]-oy)*dy1);
		int near_active=0;//neighbor is more coupled
		int near_exist=0;//neighbor exist
		double thres2=MAX(0.1, actcpl[iact]);//was 0.1
		for(int iy=-1; iy<2; iy++){
		    for(int ix=-1; ix<2; ix++){
			if(abs(ix+iy)!=1){
			    continue;/*skip self and corner */
			}
			int kact1=loc_map_get(map, mapx+ix, mapy+iy);
			if(kact1>0 && (!stuck || !stuck[kact1-1])){
			    near_exist++;
			    if(actcpl0[kact1]>thres2){//better than this one.
				near_active++;
			    }
			}
		
		    }
		}
		if(!near_exist){
		    error("This is an isolated actuator\n");
		}
		/*
		  neighbors are defined as the four pixels to the left, right,
		  top and bottom.  If some of the neighbors are active, use the
		  average of them for my value, otherwise, use the average of
		  all neighbors.
		*/
		double value=0;
		if(!near_active) value=-scl*0.1/near_exist;
		double valsum=0;
		/*part 1*/
		for(int iy=-1; iy<2; iy++){
		    for(int ix=-1; ix<2; ix++){
			if(abs(ix+iy)!=1){
			    continue;//skip self and corner
			}
			int kact1=loc_map_get(map, mapx+ix, mapy+iy);
			if(kact1>0 && (!stuck || !stuck[kact1-1])){
			    if(!near_active){
				valsum+=(px[count]=value);
				pi[count]=kact1-1;
				count++;
			    }else if(actcpl0[kact1]>actcpl[iact]){
				valsum+=(px[count]=-scl*MAX(0.1,actcpl0[kact1]));
				pi[count]=kact1-1;
				count++;
			    }
			}
		    }
		}
		
		/*part 2, matches negative sum of part 1*/
		pi[count]=iact;
		px[count]=-valsum;
		count++;
	    }/*if */
	}/*for iact */
	pp[nact]=count;
	dspsetnzmax(slavet, count);
	IND(actslave,idm,idm)=dspmulsp(slavet, slavet,"nt");
	if(NW && NW->p[idm] && 0){
	    /*Now we need to make sure NW is in the NULL
	      space of the slaving regularization, especially
	      the tip/tilt constraints. NW=NW-slavet*inv(slavet)*NW 

	      2014-08-26: This step has no performance difference and is slow
	      for large DMs. Disable.
	    */
	    dmat *H=NULL;
	    dspfull(&H, slavet, 'n',1);
	    dmat *Hinv=dpinv(H,NULL);
	    dmat *mod=NULL;
	    dmm(&mod, 0, Hinv, NW->p[idm], "nn", 1);
	    dmm(&NW->p[idm], 1, H, mod,"nn", -1);
	    dfree(H);
	    dfree(Hinv);
	    dfree(mod);
	    if(stuck || floated){
		dmat*  pNW=NW->p[idm];
		for(int iy=0; iy<NW->p[idm]->ny; iy++){
		    for(int iact=0; iact<nact; iact++){
			if((stuck && stuck[iact])||(floated && floated[iact])){
			    IND(pNW,iact,iy)=0;
			}
		    }
		}
	    }
	}
	dspfree(slavet);
    }/*idm */
    if(nslavetot==0){
	dspcellfree(actslavec);
	actslavec=NULL;
    }
    return actslavec;
}
/**
   When some actuators are stuck, zero the corresponding column in HA
*/
void act_stuck(loccell *aloc, void *HA_, const lcell *stuck){
    if(!stuck || !HA_) return; 
    cell *HA=(cell*)HA_;
    int ndm=aloc->nx;
    for(int idm=0; idm<ndm; idm++){
	if(!stuck->p[idm]){
	    continue;
	}
	const int nact=aloc->p[idm]->nloc;
	int nfit=0;
	if(HA->ny==ndm && HA->nx>1){
	    nfit=HA->nx;
	}else if(HA->ny==1 && HA->nx==ndm){
	    nfit=1;
	}else{
	    error("HA: Invalid format %ldx%ld\n", HA->nx, HA->ny);
	}
	for(int ifit=0; ifit<nfit; ifit++){
	    cell *HAi=HA->p[idm*nfit+ifit];
	    if(HAi->id==M_DBL){//dense
		dmat *hb=(dmat*)HAi;
		if(hb->nx>1 && hb->ny==aloc->p[idm]->nloc){
		    //modifying interaction matrix
		    for(int iact=0; iact<nact; iact++){
			if(stuck->p[idm]->p[iact]){
			    memset(hb->p+hb->nx*iact, 0, sizeof(double)*hb->nx);
			}
		    }
		}else if(hb->ny==1 && hb->nx==aloc->p[idm]->nloc){
		    //modifying coupling vector
		    for(int iact=0; iact<nact; iact++){
			if(stuck->p[idm]->p[iact]){
			    hb->p[iact]=0;
			}
		    }
		}else{
		    error("Invalid input: hb is %ldx%ld, nloc=%ld\n", hb->nx, hb->ny, aloc->p[idm]->nloc);
		}
	    }else if(HAi->id==M_DSP){//sparse
		dsp *ha=(dsp*)HAi;
		spint *pp=ha->p;
		double *px=ha->x;
		assert(ha->ny==aloc->p[idm]->nloc);
		for(int iact=0; iact<nact; iact++){
		    if(stuck->p[idm]->p[iact]){
			for(int ic=pp[iact]; ic<pp[iact+1]; ic++){
			    px[ic]=0;
			}
		    }
		}
	    }else{
		error("Invalid parameter: HA.id=%u\n", HAi->id);
	    }
	}
    }
}
/**
   Zero out rows of dead actuators in mode vector.
*/
void act_zero(loccell *aloc, const dcell *HB, const lcell *dead){
    if(!dead || !HB) return;
    for(int idm=0; idm<dead->nx; idm++){
	if(!dead->p[idm]){
	    continue;
	}
	const int nact=aloc->p[idm]->nloc;
	if(HB->nx!=aloc->nx){
	    error("HB: Invalid format\n");
	}
	for(int imod=0; imod<HB->ny; imod++){
	    dmat *hb=HB->p[idm+imod*HB->nx];
	    if(hb->nx!=aloc->p[idm]->nloc){
		error("hb: Invalid format\n");
	    }
	    dmat*  phb=hb;
	    for(int iact=0; iact<nact; iact++){
		if(dead->p[idm]->p[iact]){
		    for(int iy=0; iy<hb->ny; iy++){
			IND(phb,iact,iy)=0;
		    }
		}
	    }
	}
    }
}

/**
   When some actuators are float, remove the corresponding column in HA and/or HB,
   and add to neigh boring actuators. This is implemented using a second matrix and
   then add to the original matrix.*/
void act_float(loccell *aloc, dspcell **HA, const dcell *HB, const lcell *actfloat){
    if(!actfloat || ((!HA || !*HA) && !HB)) return;
    int ndm=actfloat->nx;
    dspcell *dHA=NULL;
    if(HA && *HA){
	int nfit=(*HA)->nx;
	dHA=dspcellnew(nfit,ndm);
    }
    for(int idm=0; idm<ndm; idm++){
	if(!actfloat->p[idm]) continue;
	loc_create_map(aloc->p[idm]);
	map_t *map=aloc->p[idm]->map;
	double ox=map->ox;
	double oy=map->oy;
	double dx1=1./aloc->p[idm]->dx;
	double dy1=1./aloc->p[idm]->dy;
	const double *locx=aloc->p[idm]->locx;
	const double *locy=aloc->p[idm]->locy;
	long nact=aloc->p[idm]->nloc;
	/*which floating actuator to assign for this one. */
	lmat *indfloat=lnew(4, nact);
	/*number of floating actuators that is assigned for this one. */
	long *nindfloat=mycalloc(nact,long);
	/*active neighbors of each dead act. */
	long *neighbor=mycalloc(nact,long);
	long nzmax=0;
	long *isfloat=actfloat->p[idm]->p;
	/*loop through the floating actuators */
	for(int iact=0; iact<nact; iact++){
	    if(!actfloat->p[idm]->p[iact]) continue;
	    long mapx=(long)round((locx[iact]-ox)*dx1);
	    long mapy=(long)round((locy[iact]-oy)*dy1);
	    /*find all its neighbors. */
	    for(int iy=-1; iy<2; iy++){
		for(int ix=-1; ix<2; ix++){
		    if((ix!=0 && iy!=0) || (ix==0 && iy==0)){
			continue;/*skip center and corner */
		    }
		    int kact=loc_map_get(map, mapx+ix, mapy+iy)-1;
		    if(kact>-1){
			if(!isfloat[kact]){
			    IND(indfloat,nindfloat[kact], kact)=iact;
			    nindfloat[kact]++;
			    neighbor[iact]++;
			    nzmax+=4;/*assume 1 point couples to 4 points max. */
			}
		    }
		}
	    }
	}
	if(HA && *HA){
	    dspcell*  pHA=*HA;
	    dspcell*  pdHA=dHA;
	    /*Create dHA to assign weights of floating actuators to neighbors. */
	    for(int ifit=0; ifit<(*HA)->nx; ifit++){
		spint *pp=IND(pHA,ifit,idm)->p;
		spint *pi=IND(pHA,ifit,idm)->i;
		double *px=IND(pHA,ifit,idm)->x;

		IND(pdHA,ifit,idm)=dspnew(IND(pHA,ifit,idm)->nx, IND(pHA,ifit,idm)->ny, nzmax);
		spint *pp2=IND(pdHA,ifit,idm)->p;
		spint *pi2=IND(pdHA,ifit,idm)->i;
		double *px2=IND(pdHA,ifit,idm)->x;
		long count=0;
		for(long iact=0; iact<nact; iact++){
		    pp2[iact]=count;
		    if(nindfloat[iact]){
			for(long in=0; in<nindfloat[iact]; in++){
			    long jact=IND(indfloat,in,iact);/*the floating act. */
			    double scale=1./neighbor[jact];
			    for(long ie=pp[jact]; ie<pp[jact+1]; ie++){
				if(count>=nzmax){
				    nzmax*=2;
				    dspsetnzmax(IND(pdHA,ifit,idm), nzmax);
				    pp2=IND(pdHA,ifit,idm)->p;
				    pi2=IND(pdHA,ifit,idm)->i;
				    px2=IND(pdHA,ifit,idm)->x;
				}	
				pi2[count]=pi[ie];
				px2[count]=px[ie]*scale;
				count++;
			    }
			}
		    }
		}
		pp2[nact]=count;
		dspsetnzmax(IND(pdHA,ifit,idm), count);
		/*Remove weights of floating actuatorsf from HA. */
		for(int jact=0; jact<actfloat->p[idm]->nx; jact++){
		    int iact=actfloat->p[idm]->p[jact];
		    for(int ic=pp[iact]; ic<pp[iact+1]; ic++){
			px[ic]=0;
		    }
		}
	    }/*ifit */
	}else{/*Do dense matrix */
	    for(long iact=0; iact<nact; iact++){
		if(nindfloat[iact]){
		    for(long in=0; in<nindfloat[iact]; in++){
			long jact=IND(indfloat,in,iact);/*the floating act. */
			double scale=1./neighbor[jact];
			for(int ifit=0; ifit<HB->nx; ifit++){
			    dmat *hbi=IND(HB,ifit,idm);
			    dmat*  phbi=hbi;
			    for(long ix=0; ix<hbi->nx; ix++){
				IND(phbi,ix,iact)+=IND(phbi,ix,jact)*scale;
			    }
			}
		    }
		}
	    }
	    for(int jact=0; jact<actfloat->p[idm]->nx; jact++){
		int iact=actfloat->p[idm]->p[jact];
		for(int ifit=0; ifit<HB->nx; ifit++){
		    dmat *hbi=IND(HB,ifit,idm);
		    memset(PCOL(hbi, iact), 0, hbi->nx*sizeof(double));
		}
	    }
	}
	free(nindfloat);
	lfree(indfloat);
	free(neighbor);	
    }/*idm */
    if(HA){
	dcelladd(HA, 1, dHA, 1);
	dspcellfree(dHA);
    }
}
/**
   Make DM actuator commands zero at stuck actuator locations.
*/
void act_stuck_cmd(loccell *aloc, /**<[in] Actuator grid array*/
		   const dcell *adm,   /**<[in,out] Actuator command to process*/
		   const lcell *stuck  /**<[in] List of stuck actuators*/
    ){
    if(!adm || !stuck) return;
    const int ndm=aloc->nx;
    for(int idm=0; idm<ndm; idm++){
	if(!stuck->p[idm]) continue;
	const int nact=aloc->p[idm]->nloc;
	for(int iy=0; iy<adm->ny; iy++){
	    dmat *pai=IND(adm,idm,iy);
	    dmat*  px=pai;
	    assert(pai->nx==aloc->p[idm]->nloc);
	    for(int icol=0; icol<pai->ny; icol++){
		for(int iact=0; iact<nact; iact++){
		    if(stuck->p[idm]->p[iact]){
			IND(px,iact,icol)=0;
		    }
		}
	    }
	}
    }
}
/**
   Create an interpreter that make floating actuators equal to their neighbors.
*/
dspcell* act_float_interp(loccell *aloc,  /**<[in] Actuator grid array*/
			  const lcell *actfloat/**<[in] List of floating actuators*/
    ){
    int ndm=aloc->nx;
    dspcell *out=dspcellnew(ndm, ndm);
    for(int idm=0; idm<ndm; idm++){
	loc_create_map(aloc->p[idm]);
	map_t *map=aloc->p[idm]->map;
	double ox=map->ox;
	double oy=map->oy;
	double dx1=1./aloc->p[idm]->dx;	
	double dy1=1./aloc->p[idm]->dy;
	const double *locx=aloc->p[idm]->locx;
	const double *locy=aloc->p[idm]->locy;
	long nact=aloc->p[idm]->nloc;
	/*actuator that is floating */
	long *isfloat=actfloat?(actfloat->p[idm]?actfloat->p[idm]->p():0):0;
	dsp *outit=dspnew(nact, nact, nact*4);
	double *px=outit->x;
	spint *pp=outit->p;
	spint *pi=outit->i;
	long count=0;
	for(long iact=0; iact<nact; iact++){
	    pp[iact]=count;
	    if(isfloat && isfloat[iact]){
		long mapx=(long)round((locx[iact]-ox)*dx1);
		long mapy=(long)round((locy[iact]-oy)*dy1);
		int count2=count;
		for(int idy=-1; idy<2; idy++){
		    for(int idx=-1; idx<2; idx++){
			if(abs(idx+idy)!=1){
			    continue;/*skip center and corner */
			}
			int kact=loc_map_get(map, mapx+idx, mapy+idy)-1;
			if(kact>-1){
			    pi[count]=kact;
			    px[count]=1.;
			    count++;
			}
		    }
		}
		if(count>count2){
		    double scl=1./(count-count2);
		    for(; count2<count; count2++){
			px[count2]=scl;
		    }
		}
	    }else{/*just copy over data */
		pi[count]=iact;
		px[count]=1;
		count++;
	    }
	}
	pp[nact]=count;
	out->p[idm+ndm*idm]=dsptrans(outit);
	dspfree(outit);
    }
    return out;
}

/**
   Create an interpreter that make inactive actuators equal avearage of active
   neighbors if exist or all other neighbors.
*/
dsp* act_extrap_do(loc_t *aloc,        /**<[in] Actuator grid array*/
		   const dmat *actcplc,/**<[in] Actuator coupling coefficiency*/
		   const double thres  /**<[in] Threshold of coupling to turn on interpolation*/
    ){
    dsp *out=0;
    const double *cpl=actcplc->p;
    const double *cpl0 = cpl-1;
    loc_create_map(aloc);
    map_t *map=aloc->map;
    double ox=map->ox;
    double oy=map->oy;
    double dx1=1./aloc->dx;
    double dy1=1./aloc->dy;
    const double *locx=aloc->locx;
    const double *locy=aloc->locy;
    long nact=aloc->nloc;
    dsp *outit=dspnew(nact, nact, nact*4);
    double *px=outit->x;
    spint *pp=outit->p;
    spint *pi=outit->i;
    long count=0;
    for(long iact=0; iact<nact; iact++){
	pp[iact]=count;
	if(cpl[iact]<thres){
	    long mapx=(long)round((locx[iact]-ox)*dx1);
	    long mapy=(long)round((locy[iact]-oy)*dy1);
	    int count2=count;
	    double sum=0;
	    /*first, interpolate from neighbors of higher cpl*/
	    int near_active=0;
	    //double thres2=0.1;
	    double thres2=MAX(0.1,cpl[iact]);
	    //Find active neighbors
	    for(int iy=-1; iy<2; iy++){
		for(int ix=-1; ix<2; ix++){
		    if(abs(ix)+abs(iy)==2){
			continue;//skip corner
		    }
		    //include self and neighbor
		    int kact1=loc_map_get(map, mapx+ix, mapy+iy);
		    if(kact1>0 && cpl0[kact1]>=thres2){
			near_active++;
		    }
		}
	    }
	    for(int iy=-1; iy<2; iy++){
		for(int ix=-1; ix<2; ix++){
		    if(abs(ix)+abs(iy)==2){
			continue;
		    }
		    int kact1=loc_map_get(map, mapx+ix, mapy+iy);
		    if(kact1>0){
			if(!near_active){
			    //there are no active neighbors, use the average
			    pi[count]=kact1-1;
			    sum+=(px[count]=1);
			    count++;
			}else if(cpl0[kact1] >=thres2 || (ix==0 && iy==0)){
			    //there are a few active neighbors
			    pi[count]=kact1-1;
			    sum+=(px[count]=cpl0[kact1]);
			    count++;
			}
		    }
		}
	    }
	    if(count>count2){
		double scl=1./sum;
		for(;count2<count; count2++){
		    px[count2]*=scl;
		}
	    }
	}else{
	    /*just copy over data */
	    pi[count]=iact;
	    px[count]=1;
	    count++;
	}
    }
    pp[nact]=count;
    out=dsptrans(outit);
    dspfree(outit);
	
    /*The above interpolation only propagate the value one step. Multiple the
      interpolator a few times to propagate longer.*/
    for(int i=0; i<5; i++){
	dsp *tmp=dspmulsp(out,out,"nn");
	dspfree(out);
	out=tmp;
    }
    return out;
}
dspcell* act_extrap(loccell *aloc,     /**<[in] Actuator grid array*/
		    const dcell *actcplc,/**<[in] Actuator coupling coefficiency*/
		    const double thres /**<[in] Threshold of coupling to turn on interpolation*/
    ){
    int ndm=actcplc->nx;
    dspcell *out=dspcellnew(ndm, ndm);
    for(int idm=0; idm<ndm; idm++){
	out->p[idm+ndm*idm]=act_extrap_do(aloc->p[idm], actcplc->p[idm], thres);
    }
    return out;
}
