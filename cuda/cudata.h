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
#ifndef AOS_CUDA_DATA_H
#define AOS_CUDA_DATA_H
#include <map>
#include "types.h"
#include "perf.h"
#include "wfs.h"
extern int NGPU;
extern int MAXGPU;
typedef Real ATYPE;
typedef Real GTYPE;
namespace cuda_recon{
class curecon_t;
class curecon_geom;
}
struct cuperf_t;
class cupowfs_t;
class cuwfs_t;

typedef struct cudata_t{ 
    int igpu;
    static int recongpu;
    static cuarray<int> evlgpu;
    static cuarray<int> wfsgpu;
    std::map<uint64_t, void*> memhash;/*For reuse constant GPU memory*/
    std::map<void *, int> memcount; /*Store count of reused memory*/
    void *memcache;/*For reuse temp array for type conversion.*/
    long nmemcache;
    pthread_mutex_t memmutex;
    /**<for accphi */
    void *reserve;   /**<Reserve some memory in GPU*/
    cumapcell atm;   /**<atmosphere: array of cumap_t */
    static dmat *atmscale; /**<Scaling of atmosphere due to r0 variation*/
    cumapcell dmreal;/**<DM: array of cumap_t */
    cumapcell dmproj;/**<DM: array of cumap_t */
    int nps; /**<number of phase screens*/
    /*for perfevl */
    cuperf_t perf;
    /*for wfsgrad */
    cuarray<cupowfs_t>powfs;
    static cuarray<cuwfs_t>wfs;//must be static so we can address all elements.
    /*for recon */
    cuda_recon::curecon_t *recon;
    /*for moao*/
    cuarray<cumapcell> dm_wfs;
    cuarray<cumapcell> dm_evl;
    /*for mvm*/
    curmat mvm_m;/*the control matrix*/
    ATYPE *mvm_a; /*contains act result from mvm_m*mvm_g*/
    ATYPE **mvm_a2;/*contains act copied from other gpus for sum*/
    GTYPE *mvm_g;/*the gradients copied from gpu*/
    stream_t mvm_stream;
    cudata_t():nmemcache(0),memcache(NULL),reserve(0),nps(0),powfs(0),recon(0)
	      ,mvm_a(0),mvm_a2(0),mvm_g(0){
	pthread_mutex_init(&memmutex, 0);
    }
    ~cudata_t(){
	free(memcache); memcache=0;
    }
}cudata_t;
class cudata_t;
#ifdef __APPLE__
extern pthread_key_t cudata_key;
inline cudata_t* _cudata(){
    return (cudata_t*)pthread_getspecific(cudata_key);
}
#define cudata _cudata()
#else
extern __thread cudata_t *cudata;
#endif
extern cudata_t **cudata_all;/*use pointer array to avoid misuse. */

void gpu_print_mem(const char *msg);
long gpu_get_mem(void);
/**
   switch to the next GPU and update the pointer.
*/
inline void gpu_set(int igpu){
    extern cuarray<int> GPUS;
    if(igpu>=NGPU){
	error("Invalid igpu=%d", igpu);
    }
    cudaSetDevice(GPUS[igpu]);
#ifdef __APPLE__
    pthread_setspecific(cudata_key, cudata_all[igpu]);
#else
    cudata=cudata_all[igpu];
#endif
    if(cudata->reserve){
	cudaFree(cudata->reserve);
	cudata->reserve=NULL;
    }
}
/**
   returns next available GPU. Useful for assigning GPUs to particular wfs, evl, etc.
*/
inline int gpu_next(float add=1){
    static float cur=-1;
    cur+=add;
    int ans=(int)ceil(cur);
    if(ans>=NGPU){
	ans-=NGPU;
	cur-=NGPU;
    }
    return ans;
}

#endif
