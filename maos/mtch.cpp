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
#include "mtch.h"
/**
   \file maos/mtch.c
   Setting up matched filter
*/
/**
   shift i0 without wraping into i0x1 (+1) and i0x2 (-1)
*/
static void mki0shx(double *i0x1, double *i0x2, dmat *i0, double scale){
    int nx=i0->nx;
    typedef double pcol[nx];
    pcol *i0x1p=(pcol*)i0x1;
    pcol *i0x2p=(pcol*)i0x2;
    for(int iy=0; iy<i0->ny; iy++){
	for(int ix=0; ix<i0->nx-1; ix++){
	    i0x1p[iy][ix+1]=IND(i0,ix,iy)*scale;
	    i0x2p[iy][ix]=IND(i0,ix+1,iy)*scale;
	}
    }
}

/**
  shift i0 without wraping into i0y1 (+1) and i0y2 (-1)
*/
static void mki0shy(double *i0y1, double *i0y2, dmat *i0, double scale){
    int nx=i0->nx;
    typedef double pcol[nx];
    pcol *i0y1p=(pcol*)i0y1;
    pcol *i0y2p=(pcol*)i0y2;
    for(int iy=0; iy<i0->ny-1; iy++){
	for(int ix=0; ix<i0->nx; ix++){
	    i0y1p[iy+1][ix]=IND(i0,ix,iy)*scale;
	    i0y2p[iy][ix]=IND(i0,ix,iy+1)*scale;
	}
    }
}
/**
   The routine used to generate matched filter from WFS mean short exposure
   pixel intensities.
 */
void genmtch(const PARMS_T *parms, POWFS_T *powfs, const int ipowfs){
    INTSTAT_T *intstat=powfs[ipowfs].intstat;
    const double pixthetax=parms->powfs[ipowfs].radpixtheta;
    const double pixthetay=parms->powfs[ipowfs].pixtheta;
    const double kpx=1./pixthetax;
    const double kpy=1./pixthetay;
    const double rne=parms->powfs[ipowfs].rne;
    const double bkgrnd=parms->powfs[ipowfs].bkgrnd*parms->powfs[ipowfs].dtrat;
    const double bkgrnd_res=bkgrnd*(1.-parms->powfs[ipowfs].bkgrndc);
    const int sub_i0=1;/*doesn't make any difference. */
    int ni0=intstat->i0->ny;
    if(ni0!=1 && ni0!=parms->powfs[ipowfs].nwfs){
	error("ni0 should be either 1 or %d\n", parms->powfs[ipowfs].nwfs);
    }
    const int nsa=powfs[ipowfs].saloc->nloc;
    //Prevent printing of NEA during recomputing of matched filter
    const int print_nea=intstat->mtche?0:1;
    dcellfree(intstat->mtche);
    dfree(intstat->i0sum);

    intstat->mtche=dcellnew(nsa,ni0);
    dcell *sanea=dcellnew(ni0,1);
    
    intstat->i0sum=dnew(nsa,ni0);
    dcell *i0s=intstat->i0;
    dcell* gxs=intstat->gx/*PDELL*/;
    dcell* gys=intstat->gy/*PDELL*/;
    dmat *i0sum=intstat->i0sum;
    dcell *mtche=intstat->mtche;
    if(parms->powfs[ipowfs].phytype==1){//use MF nea for recon
	dcellfree(powfs[ipowfs].saneaxy);
	powfs[ipowfs].saneaxy=dcellnew(nsa,ni0);
    }
    dcell *saneaxy=powfs[ipowfs].saneaxy;
    int nllt;
    if(parms->powfs[ipowfs].llt){
	nllt=parms->powfs[ipowfs].llt->n;
    }else{
	nllt=0;
    }
    int irot_multiplier=nllt>1?1:0;
    int nmod=3;
    int mtchcrx=0;
    int mtchcry=0;
    double shiftx=0, shifty=0;
    /*always use saved i0. cyclic shift is not good */
    /*because of the wrapped ring. */
    if(fabs(parms->powfs[ipowfs].mtchcr)>1.e-10){
	shiftx=parms->powfs[ipowfs].mtchcr;
	mtchcrx=nmod;
	nmod+=2;
	if(fabs(shiftx-1)>1.e-10){
	    error("Only constraint of 1 pixel is implemented\n");
	}
    }
    if(fabs(parms->powfs[ipowfs].mtchcra)>1.e-10){
	shifty=parms->powfs[ipowfs].mtchcra;
	mtchcry=nmod;
	nmod+=2;
	if(fabs(shifty-1)>1.e-10){
	    error("Only constraint of 1 pixel is implemented\n");
	}
    }
    
    const int i0n=powfs[ipowfs].pixpsax*powfs[ipowfs].pixpsay;
    const int mtchadp=parms->powfs[ipowfs].mtchadp;
    dmat *i0m=dnew(2,nmod);
    dmat *i0g=dnew(i0n,nmod);
    dmat *i0x1=NULL, *i0x2=NULL, *i0y1=NULL, *i0y2=NULL;
    dmat *wt=dnew(i0n,1);
    double neaspeckle=parms->powfs[ipowfs].neaspeckle/206265000.;
    if(neaspeckle>pixthetax){
	error("parms->powfs[%d].neaspeckle=%g is bigger than pixel size\n",
	      ipowfs, neaspeckle);
    }
    if(neaspeckle>0){
	warning2("powfs%d: Adding speckle noise of %.2f mas\n", ipowfs, neaspeckle*206265000);
    }
    double neaspeckle2=pow(neaspeckle,2);
    for(int ii0=0; ii0<ni0; ii0++){
	int iwfs=parms->powfs[ipowfs].wfs->p[ii0];
	double *srot=NULL;
	if(powfs[ipowfs].srot){
	    int irot=ii0*irot_multiplier;
	    srot=powfs[ipowfs].srot->p[irot]->p;
	}
	sanea->p[ii0]=dnew(nsa,2);
	dmat*  psanea=sanea->p[ii0]/*PDMAT*/;
	/*Derivative is along r/a or x/y*/
	IND(i0m,0,0)=1;
	IND(i0m,1,1)=1;
	if(!parms->powfs[ipowfs].radpix || parms->powfs[ipowfs].radgx){
	    if(mtchcrx){/*constrained x(radial) */
		double shift=pixthetax*shiftx;
		/*kp is here to ensure good conditioning */
		IND(i0m,0,mtchcrx)=shift*kpx;
		IND(i0m,0,mtchcrx+1)=-shift*kpx;
	    }
	    if(mtchcry){/*constrained y(azimuthal). */
		double shift=pixthetay*shifty;
		IND(i0m,1,mtchcry)=shift*kpy;
		IND(i0m,1,mtchcry+1)=-shift*kpy;
	    }
	}
	double i0summax=0;
	int crdisable=0;/*adaptively disable mtched filter based in FWHM. */
	int ncrdisable=0;
	for(int isa=0; isa<nsa; isa++){
	    if(parms->powfs[ipowfs].radpix && !parms->powfs[ipowfs].radgx){
		/*The derivative is along x/y, but constraint is along r/a*/
		double theta=srot[isa]; 
		if(mtchcrx){/*constrained x(radial) */
		    double shift=pixthetax*shiftx;
		    IND(i0m,0,mtchcrx)=shift*kpx*cos(theta);
		    IND(i0m,1,mtchcrx)=shift*kpx*sin(theta);
		    IND(i0m,0,mtchcrx+1)=-IND(i0m,0,mtchcrx);
		    IND(i0m,1,mtchcrx+1)=-IND(i0m,1,mtchcrx);
		}
		if(mtchcry){/*constrained y(azimuthal). */
		    double shift=pixthetay*shifty;
		    IND(i0m,0,mtchcry)=-shift*kpy*sin(theta);
		    IND(i0m,1,mtchcry)= shift*kpy*cos(theta);
		    IND(i0m,0,mtchcry+1)=-IND(i0m,0,mtchcry);
		    IND(i0m,1,mtchcry+1)=-IND(i0m,1,mtchcry);
		}
	    }

	    IND(i0sum,isa,ii0)=dsum(IND(i0s,isa,ii0));
	    if(IND(i0sum,isa,ii0)>i0summax){
		i0summax=IND(i0sum,isa,ii0);
	    }
	    if(mtchadp){
		long fwhm=dfwhm(IND(i0s,isa,ii0));
		if(fwhm>4){
		    crdisable=0;
		}else{
		    crdisable=1;
		    ncrdisable++;
		}
	    }
	    
	    double* bkgrnd2=NULL;
	    double* bkgrnd2c=NULL;
	    if(powfs[ipowfs].bkgrnd && powfs[ipowfs].bkgrnd->p[ii0*nsa+isa]){
		bkgrnd2= powfs[ipowfs].bkgrnd->p[ii0*nsa+isa]->p; 
	    }
	    if(powfs[ipowfs].bkgrndc && powfs[ipowfs].bkgrndc->p[ii0*nsa+isa]){
		bkgrnd2c= powfs[ipowfs].bkgrndc->p[ii0*nsa+isa]->p; 
	    }
	    dzero(i0g);/*don't forget to zero out */
	    adddbl(PCOL(i0g,0), 1, IND(gxs,isa,ii0)->p, i0n, 1, 0);
	    adddbl(PCOL(i0g,1), 1, IND(gys,isa,ii0)->p, i0n, 1, 0);
	    adddbl(PCOL(i0g,2), 1, IND(i0s,isa,ii0)->p, i0n, kpx, bkgrnd_res);
	    adddbl(PCOL(i0g,2), 1, bkgrnd2, i0n, 1, bkgrnd_res);
	    adddbl(PCOL(i0g,2), 1, bkgrnd2c, i0n, -1, 0);/*subtract calibration */
	    if(mtchcrx && !crdisable){
		mki0shx(PCOL(i0g,mtchcrx),PCOL(i0g,mtchcrx+1),IND(i0s,isa,ii0),kpx);
		if(sub_i0){
		    adddbl(PCOL(i0g,mtchcrx),1,IND(i0s,isa,ii0)->p, i0n, -kpx, 0);
		    adddbl(PCOL(i0g,mtchcrx+1),1,IND(i0s,isa,ii0)->p, i0n, -kpx,0);
		}
		adddbl(PCOL(i0g,mtchcrx),   1, bkgrnd2,  i0n,  1, bkgrnd_res);
		adddbl(PCOL(i0g,mtchcrx),   1, bkgrnd2c, i0n, -1, 0);
		adddbl(PCOL(i0g,mtchcrx+1), 1, bkgrnd2,  i0n,  1,bkgrnd_res);
		adddbl(PCOL(i0g,mtchcrx+1), 1, bkgrnd2c, i0n, -1, 0);
	    }
	    if(mtchcry && !crdisable){
		mki0shy(PCOL(i0g,mtchcry),PCOL(i0g,mtchcry+1),IND(i0s,isa,ii0),kpy);
		if(sub_i0){
		    adddbl(PCOL(i0g,mtchcry),1,IND(i0s,isa,ii0)->p, i0n, -kpy,0);
		    adddbl(PCOL(i0g,mtchcry+1),1,IND(i0s,isa,ii0)->p,i0n, -kpy,0);
		}
		adddbl(PCOL(i0g,mtchcry),  1, bkgrnd2,  i0n,  1, bkgrnd_res);
		adddbl(PCOL(i0g,mtchcry),  1, bkgrnd2c, i0n, -1, 0);
		adddbl(PCOL(i0g,mtchcry+1),1, bkgrnd2,  i0n,  1, bkgrnd_res);
		adddbl(PCOL(i0g,mtchcry+1),1, bkgrnd2c ,i0n, -1, 0);
	    }
	  
	    if(bkgrnd2){
		/*adding rayleigh backscatter poisson noise. */
		for(int i=0; i<i0n; i++){/*noise weighting. */
		    wt->p[i]=1./(rne*rne+bkgrnd+IND(i0s,isa,ii0)->p[i]+bkgrnd2[i]);
		}	
	    }else{
		for(int i=0; i<i0n; i++){/*noise weighting. */
		    wt->p[i]=1./(rne*rne+bkgrnd+IND(i0s,isa,ii0)->p[i]);
		}
	    }
	    if(crdisable){
		/*temporarily mark matrix is only 3 col, which effectively
		  disables constraint*/
		i0g->ny=3;
		i0m->ny=3;
	    }
	    dmat *tmp=dpinv(i0g, wt);
	    dmm(PIND(mtche,isa,ii0),0,i0m, tmp, "nn", 1);
	    dfree(tmp);
	    if(crdisable){
		/*Put old values back. */
		i0g->ny=nmod;
		i0m->ny=nmod;
	    }
	    for(int i=0; i<i0n; i++){/*noise weighting. */
		wt->p[i]=1./wt->p[i];
	    }
	    dmat *nea2=dtmcc(IND(mtche,isa,ii0), wt);
	    nea2->p[0]+=neaspeckle2;
	    nea2->p[3]+=neaspeckle2;
	    if(IND(i0sum,isa,ii0)<EPS){//zero flux
		nea2->p[0]=nea2->p[3]=pixthetax*10;
	    }
	    if(parms->powfs[ipowfs].mtchcpl==0 
	       && (!parms->powfs[ipowfs].radpix || parms->powfs[ipowfs].radgx)){
		/*remove coupling between r/a (x/y) measurements. */
		nea2->p[1]=nea2->p[2]=0;
	    }
	    IND(psanea,isa,0)=nea2->p[0];
	    IND(psanea,isa,1)=nea2->p[3];
		
	    if(parms->powfs[ipowfs].radpix && parms->powfs[ipowfs].radgx){
		//Rotate matched filter to produce grads in (x/y).
		double theta=srot[isa]; 
		drotvect(IND(mtche,isa,ii0), theta);
	    }
	    if(parms->powfs[ipowfs].phytype==1){
		if(parms->powfs[ipowfs].radpix && parms->powfs[ipowfs].radgx){
		    //Rotate NEA to (x/y)
		    double theta=srot[isa]; 
		    drotvecnn(PIND(saneaxy, isa, ii0), nea2, theta);
		}else{
		    IND(saneaxy, isa, ii0)=dref(nea2);
		}
	    }
	    dfree(nea2);
	}/*isa  */
	double siglev=parms->powfs[ipowfs].dtrat*parms->wfs[iwfs].siglev;
	if(i0summax<siglev*0.1 || i0summax>siglev){
	    warning("i0 sum to maximum of %g, wfs %d has siglev of %g\n",
		    i0summax, iwfs, siglev);
	}
	if(mtchadp){
	    info2("Mtched filter contraint are disabled for %d subaps out of %d.\n",
		  ncrdisable, nsa);
	}
    }/*ii0 */
    if(print_nea){
	info2("Matched filter sanea:\n");
	if(powfs[ipowfs].sprint){/*print nea for select subapertures.*/
	    for(int ii0=0; ii0<ni0; ii0++){
		int illt=0;
		if(ni0==parms->powfs[ipowfs].llt->n){
		    illt=ii0;
		}else if(ni0==parms->powfs[ipowfs].nwfs && parms->powfs[ipowfs].llt->n==1){
		    illt=0;
		}else{
		    error("Invalid combination\n");
		}
		info2("ii0 %d, llt %d.\n", ii0, illt);
		info2("sa index   dist   noise equivalent angle\n");
		dmat*  psanea=sanea->p[ii0]/*PDMAT*/;
		for(int ksa=0; ksa<powfs[ipowfs].sprint->p[illt]->nx; ksa++){
		    int isa=(int)powfs[ipowfs].sprint->p[illt]->p[ksa];
		    if(isa>0){
			info2("sa %4d: %5.1f m, (%6.2f, %6.2f) mas\n", 
			      isa, powfs[ipowfs].srsa->p[illt]->p[isa], 
			      sqrt(IND(psanea,isa,0))*206265000,
			      sqrt(IND(psanea,isa,1))*206265000);
		    }
		}
	    }
	}else{
	    for(int ii0=0; ii0<ni0; ii0++){
		info2("ii0=%d:\n",ii0);
		dmat*  psanea=sanea->p[ii0]/*PDMAT*/;
		double dsa=powfs[ipowfs].saloc->dx;
		double llimit=-dsa/2;
		double ulimit=dsa/2;
		info2("sa index: radius   noise equivalent angle\n");
		for(int isa=0; isa<nsa; isa++){
		    double locx=powfs[ipowfs].saloc->locx[isa];
		    double locy=powfs[ipowfs].saloc->locy[isa];
		    if(nsa<10 || (locx>0&&locy>llimit&&locy<ulimit)){
			info2("sa %5d: %5.1f m, (%6.2f, %6.2f) mas\n", 
			      isa, locx, sqrt(IND(psanea,isa,0))*206265000,
			      sqrt(IND(psanea,isa,1))*206265000);
		    }
		}/*isa  */
	    }/*ii0 */
	}
    }
    if(parms->powfs[ipowfs].phytype==1 && parms->save.setup){
	writebin(sanea, "powfs%d_sanea", ipowfs);
    }
    if(parms->powfs[ipowfs].phytype==1 && parms->recon.glao && ni0>0){
	info2("Averaging saneaxy of different WFS for GLAO mode\n");
	dcell *saneaxy2=dcellnew(nsa, 1);
	double scale=1./ni0;
	for(int isa=0; isa<nsa; isa++){
	    for(int ii0=0; ii0<ni0; ii0++){
		dadd(&saneaxy2->p[isa], 1, IND(saneaxy, isa, ii0), scale);
	    }
	}
	dcellfree(powfs[ipowfs].saneaxy);
	powfs[ipowfs].saneaxy=saneaxy2;
    }
    dcellfree(sanea);
    dfree(i0m);
    dfree(i0g);
    dfree(i0x1); dfree(i0x2); dfree(i0y1); dfree(i0y2);
    dfree(wt);
}