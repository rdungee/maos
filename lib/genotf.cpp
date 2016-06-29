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
#include "genotf.h"
/**
private data struct to mark valid pairs of points.  */
typedef struct T_VALID{
    long n;
    long (*loc)[2];
}T_VALID;
/**
 Wrap the data to genotf to have multi-thread capability.*/
typedef struct GENOTF_T{
    cmat **otf;
    loc_t *loc;     /**<the common aperture grid*/
    const dmat *amp;    /**<The amplitude map of all the (sub)apertures*/
    const dmat *opdbias;/**<The static OPD bias. */
    const dmat *area;   /**<area of a (sub)aperture*/
    double thres;/**<threshold to consider an (sub)aperture as full*/
    double wvl;  /**<The wavelength. only needef if opdbias is not null*/
    long ncompx; /**<Size of OTF*/
    long ncompy; /**<Size of OTF*/
    long nsa;    /**<Number of (sub)apertures*/
    long pttr;   /**<Remove piston/tip/tilt*/
    const dmat *B;
    const T_VALID *pval;
    long isafull;
    const cmat *otffull;
}GENOTF_T;
/**
   Remove tip/tilt from the covariance matrix.
*/
static dmat* pttr_B(const dmat *B,   /**<The B matrix. */
		    loc_t *loc,       /**<The aperture grid*/
		    const double *amp /**<The amplitude map*/
		   ){
    if(!amp) error("amplitude map has to be not empty to remove pistion/tip/tilt\n");
    double *locx=loc->locx;
    double *locy=loc->locy;
    int nloc=loc->nloc;

    dmat *B2=dnew(nloc, nloc);
    dmat*  BP=B2;
  
    double *mod[3];
    dmat *mcc=dnew(3,3);/*modal cross coupling matrix. */
    dmat*  cc=mcc;
 
    mod[0]=NULL;
    mod[1]=locx;
    mod[2]=locy;
    for(int im=0; im<3;im++){
	for(int jm=im;jm<3;jm++){
	    IND(cc,im,jm)=dotdbl(mod[im], mod[jm], amp, nloc);
	    if(im!=jm)
		IND(cc,jm,im)=IND(cc,im,jm);
	}
    }
    dinvspd_inplace(mcc);
    dmat *M   =dnew(nloc, 3);/*The tip/tilt modal matrix */
    dmat *MW  =dnew(nloc, 3);/*M*W */
    dmat *MCC =dnew(3, nloc);/*M*inv(M'*W*M) */
    dmat *Mtmp=dnew(3, nloc);/*B'*MW; */
 
    for(long iloc=0; iloc<nloc; iloc++){
	M->p[iloc]=1;
    }
    memcpy(M->p+nloc, locx, nloc*sizeof(double));
    memcpy(M->p+nloc*2, locy, nloc*sizeof(double));
    for(long iloc=0; iloc<nloc; iloc++){
	MW->p[iloc]=amp[iloc];
	MW->p[iloc+nloc]=amp[iloc]*locx[iloc];
	MW->p[iloc+nloc*2]=amp[iloc]*locy[iloc];
    }
    /* MCC = - cci' *M' */
    dmm(&MCC, 0, mcc, M,  "tt", -1);
    dmat*  pMCC=MCC;
    /* Mtmp =  MW' * B  */
    dmm(&Mtmp, 0, MW, B, "tn", 1);
    /*Remove tip/tilt from left side*/
    dmat*  pMtmp=Mtmp;
    for(long iloc=0; iloc<nloc; iloc++){
	double tmp1=IND(pMtmp,0,iloc);
	double tmp2=IND(pMtmp,1,iloc);
	double tmp3=IND(pMtmp,2,iloc);
	for(long jloc=0; jloc<nloc; jloc++){
	    IND(BP,jloc,iloc)=IND(B,jloc,iloc)+
		(IND(pMCC,0,jloc)*tmp1
		 +IND(pMCC,1,jloc)*tmp2
		 +IND(pMCC,2,jloc)*tmp3);
	}
    }
    /* Mtmp = MW' * BP' */
    dmm(&Mtmp, 0, MW, B2, "tt", 1);
    /*Remove tip/tilt from right side*/
    for(long iloc=0; iloc<nloc; iloc++){
	double tmp1=IND(pMCC,0,iloc);
	double tmp2=IND(pMCC,1,iloc);
	double tmp3=IND(pMCC,2,iloc);
	for(long jloc=0; jloc<nloc; jloc++){
	    IND(BP,jloc,iloc)+=
		tmp1*IND(pMtmp,0,jloc)
		+tmp2*IND(pMtmp,1,jloc)
		+tmp3*IND(pMtmp,2,jloc);
	}
    }
    dfree(mcc);
    dfree(M);
    dfree(MW);
    dfree(MCC);
    dfree(Mtmp);
    return B2;
}
/**
   Generate OTF from the B or tip/tilted removed B matrix.
*/
static void genotf_do(cmat **otf, long pttr, long notfx, long notfy, 
		      loc_t *loc, const double *amp, const double *opdbias, double wvl,
		      const dmat* B,  const T_VALID *pval){
    long nloc=loc->nloc;
    dmat *B2;
    if(pttr){/*remove p/t/t from the B matrix */
	B2=pttr_B(B,loc,amp);
    }else{
	B2=ddup(B);/*duplicate since we need to modify it. */
    }
    dmat*  BP=B2;
    if(!*otf){
	*otf=cnew(notfx,notfy);
    }
    cmat* OTF=*otf;
    /*Do the exponential.*/
    double k2=pow(2*M_PI/wvl,2);
    double *restrict BPD=mymalloc(nloc,double);
    for(long iloc=0; iloc<nloc; iloc++){
	for(long jloc=0; jloc<nloc; jloc++){
	    IND(BP,jloc,iloc)=exp(k2*IND(BP,jloc,iloc));
	}
	BPD[iloc]=pow(IND(BP,iloc,iloc), -0.5);
    }
    double otfnorm=0;
    if(amp){
	for(long iloc=0; iloc<nloc; iloc++){
	    otfnorm+=amp[iloc]*amp[iloc];
	}
    }else{
	otfnorm=nloc;
    }
    otfnorm=1./otfnorm;

    dcomplex wvk=COMPLEX(0, 2.*M_PI/wvl);
    for(long jm=0; jm<notfy; jm++){
	for(long im=0; im<notfx; im++){
	    long (*jloc)[2]=pval[jm*notfx+im].loc;
	    double tmp1,tmp2; dcomplex tmp3;
	    dcomplex tmp=0.;
	    for(long iloc=0; iloc<pval[jm*notfx+im].n; iloc++){
		long iloc1=jloc[iloc][0];/*iloc1 is continuous. */
		long iloc2=jloc[iloc][1];/*iloc2 is not continuous. */
		tmp1=BPD[iloc1]*IND(BP, iloc2, iloc1);
		tmp2=BPD[iloc2];
		if(amp){
		    tmp1*=amp[iloc1];
		    tmp2*=amp[iloc2];
		}
		if(opdbias){
		    tmp3=cexp(wvk*(opdbias[iloc1]-opdbias[iloc2]));
		    tmp+=tmp1*tmp2*tmp3;
		}else{
		    tmp+=tmp1*tmp2;
		}
	    }
	    IND(OTF,im,jm)=tmp*otfnorm;
	}
    }
    free(BPD);
    dfree(B2);
}
/**
   A wrapper to execute pttr parallel in pthreads
 */
static void genotf_wrap(thread_t *info){
    GENOTF_T *data=(GENOTF_T*)info->data;
    const int nsa=data->nsa;
    cmat**otf=(cmat**)data->otf;
    loc_t *loc=data->loc;
    const long nxsa=loc->nloc;
    const double wvl=data->wvl;
    const long ncompx=data->ncompx;
    const long ncompy=data->ncompy;
    const dmat *area=data->area;
    const double thres=data->thres;
    const cmat *otffull=data->otffull;
    const dmat *amp=data->amp;
    const dmat *opdbias=data->opdbias;
    const long pttr=data->pttr;
    const dmat *B=data->B;
    const T_VALID *pval=data->pval;
    assert(!area || area->nx*area->ny==nsa);
    assert(!amp || amp->nx*amp->ny==nxsa*nsa);
    assert(!opdbias || opdbias->nx*opdbias->ny==nxsa*nsa);
    for(int isa=info->start; isa<info->end; isa++){
	if(!detached && nsa>10 && info->ithread==0){
	    //info2("%3ld%%\b\b\b\b\b\b\b", (isa-info->start)*100/(info->end-info->start));
  	    info2("%6ld of %6d\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b", isa*(nsa/info->end), nsa);
	}
	const double *opdbiasi=NULL;
	if(opdbias){
	    opdbiasi=opdbias->p+isa*nxsa;
	}else{
	    opdbiasi=NULL;
	}
	if(otffull && (!area || area->p[isa]>thres)){
	    ccp(&otf[isa],otffull);/*just copy the full array */
	}else if(!area || area->p[isa]>0){ 
	    genotf_do(&otf[isa],pttr,ncompx,ncompy,loc,amp?amp->p+isa*nxsa:NULL,opdbiasi,wvl,B,pval);
	}
    }
    if(!detached && nsa>10) info2("Thread %ld done\n", info->ithread);
}
/**
   Generate pairs of overlapping points for structure function.  

   2010-11-08: removed amp. It caused wrong otf because it uses the amp of the
   first subaperture to build pval, but this one is not fully illuminated. 
 */
static T_VALID *gen_pval(long notfx, long notfy, loc_t *loc, double xsep, double ysep){
    long nloc=loc->nloc;
    double *locx=loc->locx;
    double *locy=loc->locy;
    const long pvaltot=notfx*notfy*nloc*2;
    typedef long long2[2];
    long2 *pval0=mymalloc(pvaltot,long2);
    if(!pval0){
	error("malloc for %ld failed\n", pvaltot);
    }
    T_VALID *pval=mymalloc(notfx*notfy,T_VALID);
    long count=0,count2;
    loc_create_map(loc);
    map_t *map=loc->map;
    long notfx2=notfx/2;
    long notfy2=notfy/2;
    double dx1=1./loc->dx;
    double dy1=1./loc->dy;
    for(long jm=0; jm<notfy; jm++){
	long jm2=(jm-notfy2);/*peak in the center */
	/*long jm2=jm<notfy2?jm:jm-notfy;//peak in the corner */
	for(long im=0; im<notfx; im++){
	    long im2=(im-notfx2);
	    /*long im2=im<notfx2?im:im-notfx; */
	    count2=count;
	    for(long iloc=0; iloc<loc->nloc; iloc++){
		long iy=(long)round((locy[iloc]+jm2*ysep-map->oy)*dy1);
		long ix=(long)round((locx[iloc]+im2*xsep-map->ox)*dx1);
		long iloc2=(long)loc_map_get(map, ix, iy);
		if(iloc2>0){
		    pval0[count][0]=iloc;
		    pval0[count][1]=iloc2-1;
		    count++;
		}
	    }
	    pval[jm*notfx+im].loc=pval0+count2;
	    pval[jm*notfx+im].n=count-count2;
	}
    }
    if(count>pvaltot){
	error("count=%ld > pvaltot=%ld\n", count, pvaltot);
    }
    /*loc_free_map(loc);*//*do not free map. dangerous in multi-threaded envorionment. where other threads may be visiting loc->map.*/
    /*pval0=myrealloc(pval0,count*2,int); //do not realloc. will change position. */
    return pval;
}
/**
   Generate the turbulence covariance matrix B with B(i,j)=-0.5*D(x_i, x_j). We
   are neglecting the DC part since the OTF only depends on the structure
   function.  */
static dmat* genotfB(loc_t *loc, double r0, double L0){
    (void)L0;
    long nloc=loc->nloc;
    double *locx=loc->locx;
    double *locy=loc->locy;
    dmat *B=dnew(nloc, nloc);
    const double coeff=6.88*pow(2*M_PI/0.5e-6,-2)*pow(r0,-5./3.)*(-0.5);
    for(long i=0; i<nloc; i++){
	for(long j=i; j<nloc; j++){
	    double rdiff2=pow(locx[i]-locx[j],2)+pow(locy[i]-locy[j],2);
	    IND(B,i,j)=IND(B,j,i)=coeff*pow(rdiff2,5./6.);
	}
    }
    return B;
}
/**
   Generate OTFs for an aperture or multiple subapertures. ALl these apertures
   must share the same geometry, but may come with different amplitude map and/or
   OPD biasas. if pttr is 1, the OTF will have tip/tilt removed. make r0 to
   infinity to build diffraction limited OTF. make r0 to infinity and opdbias to
   none null to build OTF for a static map.*/

void genotf(cmat **otf,    /**<The otf array for output*/
	    loc_t *loc,    /**<the aperture grid (same for all apertures)*/
	    const dmat *amp,    /**<The amplitude map of all the (sub)apertures*/
	    const dmat *opdbias,/**<The static OPD bias (complex part of amp). */
	    const dmat *area,   /**<normalized area of the (sub)apertures*/
	    double thres,  /**<The threshold to consider a (sub)aperture as full*/
	    double wvl,    /**<The wavelength. only needef if opdbias is not null*/
	    double dtheta, /**<Sampling of PSF.*/
	    const dmat *cov,/**<The covariance. If not supplied use r0 for kolmogorov spectrum.*/
	    double r0,     /**<Fried parameter*/
	    double l0,     /**<Outer scale*/
	    long ncompx,   /**<Size of OTF*/
	    long ncompy,   /**<Size of OTF*/
	    long nsa,      /**<Number of (sub)apertures*/
	    long pttr      /**<Remove piston/tip/tilt*/
	     ){
    if(amp && loc->nloc*nsa!=amp->nx*amp->ny){
	error("loc and amp mismatch. loc->nloc=%ld, amp is %ldx%ld, nsa=%ld\n", loc->nloc, amp->nx, amp->ny, nsa);
    }else if(cov && (amp->nx!=cov->nx || cov->nx!=cov->ny)){
    	error("loc and cov mismatch\n");
    }else if(nsa<1 || ncompx<1 || ncompy<1){
	error("nsa, ncompx, ncompy has to be at least 1\n");
    }
    /*creating pairs of points that both exist with given separation*/
    double duxwvl=wvl/(dtheta*ncompx);
    double duywvl=wvl/(dtheta*ncompy);
    T_VALID *pval=gen_pval(ncompx, ncompy, loc, duxwvl, duywvl);/*returns T_VALID array. */
    /* Generate the B matrix. */
    dmat *B=cov?(dmat*)cov:genotfB(loc, r0, l0);
    cmat *otffull=NULL;
    const long nloc=loc->nloc;
    long isafull=-1;
    if(!opdbias && nsa>1){
	double maxarea=0;
	for(long isa=0; isa<nsa; isa++){
	    if(area->p[isa]>maxarea){
		maxarea=area->p[isa];
		isafull=isa;
	    }
	}
	if(isafull>0){
	    genotf_do(&otffull,pttr,ncompx,ncompy,loc,amp?amp->p+isafull*nloc:NULL,NULL,wvl,B,pval);
	}
    }
    
    GENOTF_T data={otf, loc, amp, opdbias, area, thres, wvl, ncompx, ncompy, nsa, pttr, B, pval, isafull, otffull};
    /*
    data.otf=otf;
    data.loc=loc;
    data.amp=amp;
    data.opdbias=opdbias;
    data.area=area;
    data.thres=thres;
    data.wvl=wvl;
    data.ncompx=ncompx;
    data.ncompy=ncompy;
    data.nsa=nsa;
    data.pttr=pttr;
    data.B=B;
    data.pval=pval;
    data.isafull=isafull;
    data.otffull=otffull;
    */
    thread_t info[NCPU];
    thread_prep(info, 0, nsa, NCPU, genotf_wrap, &data);
    CALL_THREAD(info, 1);
    cfree(otffull);
    if(!cov) dfree(B);
    free(pval[0].loc);
    free(pval);
}

/**
   Average spatially the 4-d covariance function to create a 2-d covariance
   function. For OPD f defined on points x (2-d coordinate), the 4-d covariance
   is simply <f'f> where f is vector form of the OPD and the average is over
   time. The 2-d covariance is additionally averaged over all the points so that
   B(r)=<f(x)'f(x+r)>_x,t To compute B, we first figure out the number of
   overlapping pairs of points for each r and then compute the averaging. When
   the amplitude is less than the threshold, the point does not count.*/

dmat *mk2dcov(loc_t *loc, const dmat *amp, double ampthres, const dmat *cov, int norm){
    if(loc->nloc!=cov->nx || loc->nloc!=cov->ny){
	error("loc and cov does not match. loc->nloc=%ld, cov is %ldx%ld\n", loc->nloc, cov->nx, cov->ny);
    }
    double xmin,xmax,ymin,ymax;
    long nloc=loc->nloc;
    double *locx=loc->locx;
    double *locy=loc->locy;
    dmaxmin(locx, nloc, &xmax, &xmin);
    dmaxmin(locy, nloc, &ymax, &ymin);  
    double dx1=1./loc->dx;
    double dy1=1./loc->dy;
    long ncovx=(long) round((xmax-xmin)*dx1)*2;
    long ncovy=(long) round((ymax-ymin)*dy1)*2;
    dmat *cov2d=dnew(ncovx, ncovy);
    /*the following is adapted from gen_pval*/
    loc_create_map(loc);
    map_t *map=loc->map;
    long ncovx2=ncovx/2;
    long ncovy2=ncovy/2;
    long *map_x=mymalloc(loc->nloc,long);
    long *map_y=mymalloc(loc->nloc,long);
    for(long iloc=0; iloc<loc->nloc; iloc++){
	map_x[iloc]=(long)round((locx[iloc]-map->ox)*dx1);
	map_y[iloc]=(long)round((locy[iloc]-map->oy)*dy1);
    }
    for(long jm=0; jm<ncovy; jm++){
	long jm2=(jm-ncovy2);//peak in the center 
	/*long jm2=jm<ncovy2?jm:jm-ncovy;//peak in the corner */
	for(long im=0; im<ncovx; im++){
	    long im2=(im-ncovx2);//peak in the center 
	    /*long im2=im<ncovx2?im:im-ncovx; //peak in the corner */
	    long count=0;
	    double acc=0;
	    for(long iloc=0; iloc<loc->nloc; iloc++){
		if(amp && amp->p[iloc]<ampthres) continue;
		long ix=map_x[iloc]+im2;
		long iy=map_y[iloc]+jm2;
		long iloc2=(long)loc_map_get(map, ix, iy);
		if(iloc2>0 && (!amp || amp->p[iloc2]>=ampthres)){
		    acc+=IND(cov,iloc2-1,iloc);
		    count++;
		}
	    }
	    if(count>0){
		if(norm){/*compute the covariance*/
		    IND(cov2d,im,jm)=acc/count;
		}else{/*compute approximate PSD.*/
		    IND(cov2d,im,jm)=acc;
		}
	    }
	}
    }
    free(map_x);
    free(map_y);
    return cov2d;
}
