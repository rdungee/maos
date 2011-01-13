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

#include <sys/file.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <dirent.h>
#include "common.h"
#include "thread.h"
#include "shm.h"
#include "process.h"
#include "turbulence.h"
#include "path.h"
#include "dmat.h"
#include "cmat.h"
#include "loc.h"
#include "fft.h"
#include "hashlittle.h"
#include "fractal.h"
#include "mathmisc.h"
#include "nr/nr.h"
#include "locbin.h"
#include "misc.h"
#include "cellarr.h"
#include "sys/daemonize.h"
int disable_atm_shm=0;
/**
 *  \file turbulence.c
 *  Contains routines to generate atmospheric turbulence screens
 */
enum{
    T_VONKARMAN=0,
    T_FRACTAL,
    T_BIHARMONIC
};

/**
 * hash the data to get a unique file name
 */
static char *fnatm(GENSCREEN_T *data){
    uint32_t key;
    key=hashlittle(data->rstat, sizeof(rand_t), 0);//contains seed
    key=hashlittle(data->wt, sizeof(double)*data->nlayer, key);
    key=hashlittle(&data->dx, sizeof(double), key);
    key=hashlittle(&data->r0, sizeof(double), key);
    key=hashlittle(&data->l0, sizeof(double), key);
    key=hashlittle(&data->nx, sizeof(long), key);
    key=hashlittle(&data->ny, sizeof(long), key);
    key=hashlittle(&data->nlayer, sizeof(long), key);
    key=hashlittle(&data->ninit, sizeof(long), key);

    char fnshm[NAME_MAX];
    snprintf(fnshm,PATH_MAX,"%s/.aos/atm", HOME);
    if(!exist(fnshm)) mymkdir("%s", fnshm);
    remove_file_older(fnshm, 30*24*3600);
    char *types[]={"vonkarman","fractal","biharmonic"};
    snprintf(fnshm,PATH_MAX,"%s/.aos/atm/maos_%s_%ld_%ldx%ld_%g_%ud.bin",
	     HOME,types[data->method],data->nlayer,data->nx,data->ny,data->dx,key);
    return strdup(fnshm);
}
/**
 *   Generate turbulence screens all in memory
 */
static void spect_screen_do(GENSCREEN_T *data){
    const dmat *spect=data->spect;
    rand_t *rstat=data->rstat;
    
    const long m=data->nx;
    const long n=data->ny;
    const long nlayer=data->nlayer;
    map_t** screen=data->screen;
    const double *wt=data->wt;
    cmat* cspect=cnew(m,n);
    cfft2plan(cspect,-1);
 repeat:
    LOCK(data->mutex_ilayer);
    int ilayer=data->ilayer;
    data->ilayer+=2;//generate 2 layers at a time.
    if(ilayer>=nlayer){
	UNLOCK(data->mutex_ilayer);
	cfree(cspect);
	return ;
    }
    /*We generate random numbers inside mutex lock to make
      sure the random sequence is repeatable.*/
    
    for(long i=0; i<m*n; i++){
	cspect->p[i]=(randn(rstat)+I*randn(rstat))*spect->p[i];
    }
    UNLOCK(data->mutex_ilayer);
    cfft2(cspect,-1); 	
    double wt1=sqrt(wt[ilayer]);
    double *p=screen[ilayer]->p;

    for(long i=0; i<m*n; i++){
	p[i]=creal(cspect->p[i])*wt1;
    }
    
    if(ilayer+1<nlayer){
	double *q=screen[ilayer+1]->p;
	double wt2=sqrt(wt[ilayer+1]);
	for(long i=0; i<m*n; i++){
	    q[i]=cimag(cspect->p[i])*wt2;
	}
    }
    goto repeat;
}
/**
 * Geneate the screens sequentially and appends to file. Handles large screens well
 * without using the full storage.
 */
static void spect_screen_save(cellarr *fc, GENSCREEN_T *data){
    rand_t *rstat = data->rstat;
    dmat *spect   = data->spect;
    double* wt    = data->wt;
    int nlayer    = data->nlayer;
    dcell *dc     = dcellnew(2,1);
    long nx = data->nx;
    long ny = data->ny;
    double dx=data->dx;
    dc->p[0] = dnew(nx, ny);
    dc->p[1] = dnew(nx, ny);
    dcell_fft2plan(dc, -1, data->nthread);
    double *restrict p1=dc->p[0]->p;
    double *restrict p2=dc->p[1]->p;
    char header[1024];
    double ox=-nx/2*dx;
    double oy=-ny/2*dx;
    snprintf(header, 1024, "ox=%.15g\noy=%.15g\ndx=%.15g\nh=%.15g\nvx=%.15g\nvy=%.15g\n",
	     ox, oy, dx, 0., 0., 0.);
    dc->p[0]->header=strdup(header);
    dc->p[1]->header=strdup(header);
    for(int ilayer=0; ilayer<nlayer; ilayer+=2){
	double tk1=myclockd();
	for(long i=0; i<nx*ny; i++){
	    p1[i]=randn(rstat)*spect->p[i];
	    p2[i]=randn(rstat)*spect->p[i];
	}
	double tk2=myclockd();
	dcell_fft2(dc, -1);
	dscale(dc->p[0], sqrt(wt[ilayer]));
	double tk3=myclockd();
	cellarr_dmat(fc, dc->p[0]);
	if(ilayer+1<nlayer){
	    dscale(dc->p[1], sqrt(wt[ilayer+1]));
	    cellarr_dmat(fc, dc->p[1]);
	}
	double tk4=myclockd();
	info("%d: Randn: %.2f FFT: %.2f Save: %.2f\n", ilayer, tk2-tk1, tk3-tk2, tk4-tk3);
    }
    dcellfree(dc);
}
/**
 * Generates multiple screens from spectrum.
 */
static map_t** create_screen(GENSCREEN_T *data, 
			     void (*funsave)(cellarr *fc, GENSCREEN_T *data),
			     void (*funmem)(GENSCREEN_T *data)){
    dmat *spect=data->spect;
    map_t **screen;
    long nlayer=data->nlayer;
    spect->p[0]=0;//zero out piston.
    if(data->share){//shared with file
	char *fnshm=fnatm(data);
	char fnlock[PATH_MAX]; 
	snprintf(fnlock, PATH_MAX, "%s.lock", fnshm);
	dcell *in=NULL;
	while(!in){
	    if(!exist(fnlock)){//when fnlock exists, the data in fnshm is not good.
		info("Trying to read %s\n", fnshm);
		in=dcellread_mmap(fnshm);
	    }else{
		info("Will not read since %s exists\n", fnlock);
	    }
	    if(!in){
		info("Trying to create %s\n", fnshm);
		int fd=lock_file(fnlock, 0);//non blocking.
		if(fd>=0){//succeed to lock file. Another process is running.
		    cellarr *fc = cellarr_init(nlayer, "%s", fnshm); 
		    funsave(fc, data);
		    cellarr_close(fc);
		    remove(fnlock);
		    close(fd);
		}else{//some other process is working on the data.
		    //wait for the lock to release.
		    fd=open(fnlock, O_RDONLY);
		    //this lock will block.
		    if(fd>-1){
			flock(fd, LOCK_EX);
			flock(fd, LOCK_UN);
		    }
		    remove(fnlock);
		}
	    }
	}
	int nlayer2;
	screen=dcell2map(&nlayer2, in);
	assert(nlayer==nlayer2);
	dcellfree(in);
	free(fnshm);
    }else{
  	screen=calloc(nlayer,sizeof(map_t*));
	long nx = data->nx;
	long ny = data->ny;
	double dx = data->dx;
	for(int ilayer=0; ilayer<nlayer; ilayer++){
	    screen[ilayer]=mapnew(nx, ny, dx, NULL);
	}
	data->screen=screen;
	PINIT(data->mutex_ilayer);
	CALL(funmem, data, data->nthread);
    }
    return screen;
}
/**
 *   Generate vonkarman screens from turbulence statistics.
 */
map_t** vonkarman_screen(GENSCREEN_T *data){
    data->method=T_VONKARMAN;
    char fnspect[PATH_MAX];
    mymkdir("%s/.aos/spect/",HOME);
    snprintf(fnspect,PATH_MAX,"%s/.aos/spect/spect_%ldx%ld_dx1_%g_r0%g_L0%g.bin",
	     HOME,data->nx,data->ny,1./data->dx,data->r0,data->l0);
    if(exist(fnspect)){
	data->spect=dread("%s", fnspect);
    }else{
	info2("\nGenerating spect..."); TIC; tic;
	data->spect=turbpsd(data->nx,data->ny,data->dx,data->r0,data->l0,0.5); toc2("done");
	dwrite(data->spect,"%s",fnspect);
    }
    map_t **screen=create_screen(data, spect_screen_save, spect_screen_do);
    dfree(data->spect);
    return(screen);
}

/**
 *  Generate screens from PSD with power of 12/3 instead of 11/3.
 */
map_t** biharmonic_screen(GENSCREEN_T *data){
    data->method=T_BIHARMONIC;
    info2("\nGenerating spect..."); TIC; tic;
    data->spect=turbpsd_full(data->nx,data->ny,data->dx,data->r0,data->l0,-2,0.5); toc2("done");
    map_t **screen=create_screen(data, spect_screen_save, spect_screen_do);
    dfree(data->spect);
    return(screen);
}
/**
 * Generate one screen at a time and save to file
 */
static void fractal_screen_save(cellarr *fc, GENSCREEN_T *data){
    long nx=data->nx;
    long ny=data->ny;
    dmat *dm = dnew(data->nx, data->ny);
    for(int ilayer=0; ilayer<data->nlayer; ilayer++){
	drandn(dm, 1, data->rstat);
	double r0i=data->r0*pow(data->wt[ilayer], -3./5.);
	fractal(dm->p, nx, ny, data->dx, r0i, data->l0, data->ninit);
	remove_piston(dm->p, nx*ny);
	cellarr_dmat(fc, dm);
    }
    dfree(dm);
}
static void fractal_screen_do(GENSCREEN_T *data){
    rand_t *rstat=data->rstat;
    map_t** screen=data->screen;
    const double *wt=data->wt;
    long nx=screen[0]->nx;
    long ny=screen[0]->ny;
 repeat:
    LOCK(data->mutex_ilayer);
    int ilayer=data->ilayer;
    data->ilayer++;
    if(ilayer>=data->nlayer){
	UNLOCK(data->mutex_ilayer);
	return;
    }
    drandn((dmat*)screen[ilayer], 1, rstat);
    UNLOCK(data->mutex_ilayer);
    double r0i=data->r0*pow(wt[ilayer], -3./5.);
    //info("r0i=%g\n", r0i);
    fractal(screen[ilayer]->p, nx, ny, screen[0]->dx, r0i, data->l0, data->ninit);
    remove_piston(screen[ilayer]->p, nx*ny);
    goto repeat;
}

/**
 * Generate Fractal screens. Not good statistics.
 */

map_t **fractal_screen(GENSCREEN_T *data){
    data->method=T_FRACTAL;
    return create_screen(data, fractal_screen_save, fractal_screen_do);
}

/**
 *  Compute the covariance for separation of r, and put the values in cov. In
 *  kolmogorov spectrum, the variance are defined as half of the structure
 *  function between two points separated by rmax.
 */

dmat* turbcov(dmat *r, double rmax, double r0, double L0){
    double tg1=tgamma(11./6) * pow(24./5 * tgamma(6./5.), 5./6.)
	* pow(2 * M_PI/0.5e-6, -2) / pow(M_PI, 8./3.);
    double vkcoeff  = tg1 / pow(2, 5./6.);    
    double vkcoeff0 = tg1 * tgamma(5./6.) / 2 ;//for variance
    dmat *cov=dnew(r->nx, r->ny);
    long n=r->nx*r->ny;
    if(isinf(L0)){//kolmogorov.
	const double power=5./3.;
	double coeff=6.88*pow(2*M_PI/0.5e-6, -2) * pow(r0, -power);
	double sigma2=0.5*coeff*pow(rmax, power);
	for(long i=0; i<n; i++){
	    cov->p[i]=sigma2-0.5*coeff*pow(r->p[i], power);
	}
    }else{//von karman.
	const double f0=1./L0;
	const double r0f0p=pow(r0*f0, -5./3.);
	double ri, rk, rip, rkp;
	double r2pif0;	
	for(long i=0; i<n; i++){
	    if(fabs(r->p[i])<EPS){
		cov->p[i]=vkcoeff0*r0f0p;
	    }else{
		r2pif0=r->p[i]*2*M_PI*f0;
		bessik(r2pif0, 5./6., &ri, &rk, &rip, &rkp);
		cov->p[i]=vkcoeff*r0f0p*pow(r2pif0, 5./6.)*rk;
	    }
	}
    }
    return cov;
}

/**
 * Compute the turbulence spectrum at size nx*ny, with spacing dx. Notice that
 * the zero frequency component is in the corner psd->p[0].
 */

dmat *turbpsd_full(long nx,      /**<The size*/
		   long ny,      /**<The size*/
		   double dx,    /**<The sampling of spatial coordinate.*/
		   double r0,    /**<The Fried parameter*/
		   double L0,    /**<The outer scale*/
		   double slope, /**<should be -11/6 for von karman or kolmogorov
				    screens, or -2 for biharmonic screen (just
				    testing only).*/
		   double power  /**< optionally do a power of psd.*/
		   ){
    if(nx & 1 || ny & 1){
	warning("Screen is odd size.");
    }
    slope*=power;
    const double dfx=1./(nx*dx);
    const double dfy=1./(ny*dx);
    const double dfx2=dfx*dfx;
    const double L02=pow(L0,-2);
    const double scrnstr=pow(0.0229*pow(r0,-5./3.)*pow((0.5e-6)/(2.*M_PI),2)*(dfx*dfy),power);
    const int nx2=nx/2;
    const int ny2=ny/2;
    dmat *psd=dnew(nx,ny);
    for(int i=0;i<ny;i++){
	double r2y=pow((i<ny2?i:i-ny)*dfy,2);// to avoid fft shifting.
	double *psd1=psd->p+i*nx;
	for(int j=0;j<nx2;j++){
	    double r2x=j*j*dfx2;
	    psd1[j] = pow(r2x+r2y+L02,slope)*scrnstr;
	}
	for(int j=nx2;j<nx;j++){
	    double r2x=(j-nx)*(j-nx)*dfx2;
	    psd1[j] = pow(r2x+r2y+L02,slope)*scrnstr;
	}
    }
    return psd;
}

