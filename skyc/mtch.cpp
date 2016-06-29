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

#include "skyc.h"
#include "mtch.h"
/**
   \file skyc/mtch.c
   Routings regarding matched filter computation.
*/
/**
   Compute pixel intensity i0, and gradient of pixel intensity gx, gy from the
   PSF.
   \todo merge this file with maos/mtch.c and put in lib folder.
 */
void psf2i0gxgy(dmat *i0, dmat *gx, dmat *gy, dmat *psf, DTF_S *dtf){
    cmat *otf=cnew(psf->nx, psf->ny);
    cmat *otfsave=cnew(psf->nx, psf->ny);
    //cfft2plan(otf, 1);
    //cfft2plan(otf, -1);
    //cfft2plan(otfsave,1);
    ccpd(&otf, psf);/*loaded psf has peak in corner */
    cfft2i(otf, -1);/*turn to OTF, peak in corner. was 1, turn to -1 on 1/30/2013 */
    ccwm(otf, dtf->nominal);
    ccp(&otfsave, otf);
    cfft2(otf, 1);/*turn back. */
    dspmulcreal(i0->p, dtf->si, otf->p, 1);
    ccp(&otf, otfsave);
    cmat*  potf=otf;
    cmat*  potfsave=otfsave;
    /*Now derivative */
    for(int iy=0; iy<otf->ny; iy++){
	for(int ix=0; ix<otf->nx; ix++){
	    IND(potf,ix,iy)*=dtf->U->p[ix];
	    IND(potfsave,ix,iy)*=dtf->U->p[iy];
	}
    }
    cfft2(otf, 1);//was 1, changed to -1 on 1/29/2013
    cfft2(otfsave, 1);//was 1, changed to -1 on 1/29/2013
    dspmulcreal(gx->p,dtf->si,otf->p,1);
    dspmulcreal(gy->p,dtf->si,otfsave->p,1);
    cfree(otf); cfree(otfsave);
}

/**
   shift without wraping i0 into i0x1 (+1) and i0x2 (-1)
*/
static void mki0shx(double *i0x1, double *i0x2, dmat *i0, double scale){
    int nx=i0->nx;
    for(int iy=0; iy<i0->ny; iy++){
	for(int ix=0; ix<i0->nx-1; ix++){
	    i0x1[iy*nx+ix+1]=IND(i0,ix,iy)*scale;
	    i0x2[iy*nx+ix]=IND(i0,ix+1,iy)*scale;
	}
    }
}
/**
  shift without wraping i0 into i0y1 (+1) and i0y2 (-1)
*/
static void mki0shy(double *i0y1, double *i0y2, dmat *i0, double scale){
    int nx=i0->nx;
    for(int iy=0; iy<i0->ny-1; iy++){
	for(int ix=0; ix<i0->nx; ix++){
	    i0y1[(iy+1)*nx+ix]=IND(i0,ix,iy)*scale;
	    i0y2[iy*nx+ix]=IND(i0,ix,iy+1)*scale;
	}
    }
}
/**
   Compute matched filter.
 */
void mtch(dcell **mtche, dmat **sanea,
	  dcell *i0, dcell *gx, dcell *gy, double pixtheta, 
	  double rne, double bkgrnd, int cr){
    double kp=1./pixtheta;
    const double psfvar=0;
    int nmod=3;
    if(cr){
	nmod=7;
    }
    const long nsa=i0->nx;
    const long nx=i0->p[0]->nx;
    const long ny=i0->p[0]->ny;
    const long npixtot=nx*ny;
    dmat *i0m=dnew(2,nmod);
    dmat *i0g=dnew(npixtot, nmod);
    dmat *wt=dnew(npixtot, 1);
    if(!*mtche){
	*mtche=dcellnew(nsa,1);
    }
    if(!*sanea){
	*sanea=dnew(nsa*2,1);
    }
    IND(i0m,0,0)=1;
    IND(i0m,1,1)=1;

    if(cr){
	IND(i0m,0,3)=1;
	IND(i0m,0,4)=-1;
	IND(i0m,1,5)=1;
	IND(i0m,1,6)=-1;
    }
    for(long isa=0; isa<nsa; isa++){    
	dzero(i0g);
	/*kp is here to ensure good conditioning */
	adddbl(PCOL(i0g,0), 1, gx->p[isa]->p, npixtot, 1, 0);
	adddbl(PCOL(i0g,1), 1, gy->p[isa]->p, npixtot, 1, 0);
	adddbl(PCOL(i0g,2), 1, i0->p[isa]->p, npixtot, kp, 0);
	if(cr){
	    mki0shx(PCOL(i0g,3), PCOL(i0g,4), i0->p[isa], kp);
	    mki0shy(PCOL(i0g,5), PCOL(i0g,6), i0->p[isa], kp);
	}
	for(long ipix=0; ipix<npixtot; ipix++){
	    wt->p[ipix]=1./(rne*rne+bkgrnd+i0->p[isa]->p[ipix]+psfvar*i0->p[isa]->p[ipix]);
	}
	
	dmat *tmp=dpinv(i0g, wt);
	dmm(&(*mtche)->p[isa],0,i0m,tmp,"nn",1);
	dfree(tmp);
	dcwpow(wt,-1);
	dmat *nea2=dtmcc((*mtche)->p[isa], wt);
	/*Drop coupling in x/y gradients. */
	(*sanea)->p[isa]=sqrt(nea2->p[0]);
	(*sanea)->p[isa+nsa]=sqrt(nea2->p[3]);
	dfree(nea2);
    }
    dfree(wt);
    dfree(i0g);
    dfree(i0m);
}
