#include "kernel.h"
/**
   A few kernels.
*/


__global__ void set_do(float *a, float alpha, int n){
    const int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	a[i]=alpha;
    }
}

__global__ void show_do(float *a, int nx, int ny){
    const int stepx=blockDim.x * gridDim.x;
    const int stepy=blockDim.y * gridDim.y;
    for(int iy=blockIdx.y * blockDim.y + threadIdx.y; iy<ny; iy+=stepy){
	for(int ix=blockIdx.x * blockDim.x + threadIdx.x; ix<nx; ix+=stepx){
	    printf("a(%d,%d)=%g\n", ix, iy, a[ix+iy*nx]);
	}
    }
}
__global__ void add_ptt_do(float *restrict opd, float (*restrict loc)[2], 
			   int n, float pis, float tx, float ty){
    const int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	opd[i]+=pis+loc[i][0]*tx+loc[i][1]*ty;
    }
}

__global__ void add_ngsmod_do(float *restrict opd, float (*restrict loc)[2], int n, 
			      float m0, float m1, float m2, float m3, float m4,
			      float thetax, float thetay, float scale, float ht, float MCC_fcp, float alpha
			      ){
    float scale1=1.f-scale;
    const int step=blockDim.x * gridDim.x;
    for(int i=blockIdx.x * blockDim.x + threadIdx.x; i<n; i+=step){
	float x=loc[i][0];
	float y=loc[i][1];
	float xy=x*y;
	float x2=x*x;
	float y2=y*y;
	opd[i]+= alpha*(+x*m0
			+y*m1
			+m2*((x2+y2-MCC_fcp)*scale1-2*scale*ht*(thetax*x+thetay*y))
			+m3*((x2-y2)*scale1 - 2*scale*ht*(thetax*x-thetay*y))
			+m4*(xy*scale1-scale*ht*(thetay*x+thetax*y)));
    }
}
