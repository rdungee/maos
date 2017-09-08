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
#include "solve.h"
namespace cuda_recon{
Real cucg_t::solve(curcell &xout, const curcell &xin, stream_t &stream){
    Real ans;
    cgtmp.count++;
    if((ans=gpu_pcg(xout, this, precond, xin, cgtmp,
		    warm_restart, maxit, stream))>1){
	cgtmp.count_fail++;
	warning2("CG %5d(%5d) does not converge. residual=%g. maxit=%d\n", 
		 cgtmp.count, cgtmp.count_fail, ans, maxit);
    }
    return ans;
}

void cumuv_t::Forward(curcell &out, Real beta, const curcell &in, Real alpha, stream_t &stream){
    if(!M) error("M Can not be empty\n");
    if(!out){
	out=curcell(nx, 1, nxs, (int*)NULL);
    }else{
	curscale(out.M(), beta, stream);
    }
    cuspmul(out.M().P(), M, in.M().P(), 1, 'n', alpha, stream);
    if(U && V){
	curmv(Vx.P(), 0, V, in.M().P(), 't', 1, stream);
	curmv(out.M().P(), 1, U, Vx.P(), 'n', -alpha, stream);
    }
}
void cumuv_t::Trans(curcell &out, Real beta, const curcell &in, Real alpha, stream_t &stream){
    if(!M) error("M Can not be empty\n");
    if(!out){
	out=curcell(ny, 1, nys, (int*)NULL);
    }else{
	curscale(out.M(), beta, stream);
    }
    
    curscale(out.M(), beta, stream);
    cuspmul(out.M().P(), M, in.M().P(), 1, 't', alpha, stream);
    if(U && V){
	curmv(Vx.P(), 0, U, in.M().P(), 't', 1, stream);
	curmv(out.M().P(), 1, V, Vx.P(), 'n', -alpha, stream);
    }
}
void cumuv_t::Init(const MUV_T *in){
    if(!in) return;
    if(M || !in->M) error("in.M() should not be NULL and M should be NULL\n");
    dspcell *inM=dspcell_cast(in->M);
    dsp *Mc=dspcell2sp(inM);
    dmat *Uc=dcell2m(in->U);
    dmat *Vc=dcell2m(in->V);
    nx=inM->nx;
    ny=inM->ny;
    nxs=new int[nx];
    nys=new int[ny];
    for(int i=0; i<nx; i++){
	nxs[i]=inM->p[i]->nx;
    }
    for(int i=0; i<ny; i++){
	nys[i]=inM->p[i*inM->nx]->ny;
    }
    M=cusp(Mc, 1);
    cp2gpu(U, Uc);
    cp2gpu(V, Vc);
    dspfree(Mc); dfree(Uc); dfree(Vc);
    Vx=curmat(V.Ny(), 1);
}

cusolve_sparse::cusolve_sparse(int _maxit, int _warm_restart, MUV_T *_R, MUV_T *_L)
    :cucg_t(_maxit, _warm_restart){
    CR.Init(_R);
    CL.Init(_L);
}
cusolve_cbs::cusolve_cbs(spchol *_C, dmat *_Up, dmat *_Vp){
    if(!_C){
	error("C cannot be empty\n");
    }
    chol_convert(_C, 0);
    Cl=cusp(_C->Cl, 0);
    cp2gpu(Cp, _C->Cp, _C->Cl->nx, 1);
    if(_Up){
	cp2gpu(Up, _Up);
	cp2gpu(Vp, _Vp);
    }
}
Real cusolve_cbs::solve(curcell &xout, const curcell &xin, stream_t &stream){
    if(!xout) xout=xin.New();
    if(Cl.Type()==SP_CSC){
	chol_solve(xout.M().P(), xin.M().P(), stream);
    }else{
	error("To implemente\n");
    }
    if(Up){
	if(!Vr){
	    Vr=curmat(Vp.Ny(), 1);
	}
	curmv(Vr.P(), 0, Vp, xin.M().P(), 't', -1, stream);
	curmv(xout.M().P(), 1, Up, Vr.P(), 'n', 1, stream);
    }
    return 0;
}
/*solve in place*/
static __global__ void cuchol_solve_lower_do(Real *restrict y, Real *Cx, int *Cp, int *Ci, int n){
    int id=threadIdx.x;
    int nd=blockDim.x;
    extern __shared__ Real sb[];
    __shared__ Real val;
    /*first solve L\y*/
    
    for(int icol=0; icol<n; icol++){
	if(id==0){
	    y[icol]/=Cx[Cp[icol]];//divide by diagonal.
	    val=-y[icol];
	}
	__syncthreads();//this is necessary!
	for(int irow=Cp[icol]+1+id; irow<Cp[icol+1]; irow+=nd){
	    y[Ci[irow]]+=val*Cx[irow];
	}
	__syncthreads();
    }
    /*Next solve L'\y. Use reduction algorithm instead of atomic add.*/
    for(int icol=n-1; icol>-1; icol--){
	sb[id]=0;
	for(int irow=Cp[icol]+1+id; irow<Cp[icol+1]; irow+=nd){
	    sb[id]+=Cx[irow]*y[Ci[irow]];
	}
	for(int step=(blockDim.x>>1);step>0;step>>=1){
	    __syncthreads();
	    if(id<step){
		sb[id]+=sb[id+step];
	    }
	}
	if(id==0){
	    y[icol]=(y[icol]-sb[0])/Cx[Cp[icol]];
	}
	__syncthreads();//this is necessary!
    }
}
void cusolve_cbs::chol_solve(Real *out, const Real *in, stream_t &stream){
    if(!Cl || !Cp) error("Invalid\n");
    int n=Cl.Nx();
    if(!y){
	y=curmat(Cl.Nx(), 1);
    }
    perm_f_do<<<DIM(n, 256),0,stream>>>(y.P(), in, Cp.P(), n);
    //only 1 block for synchronization. //todo: improve the implementation.
    const int NTH=256;
    cuchol_solve_lower_do<<<1,NTH, NTH*sizeof(Real),stream>>>(y.P(), Cl.Px(), Cl.Pp(), Cl.Pi(), n); 
    perm_i_do<<<DIM(n, 256),0,stream>>>(out, y.P(), Cp.P(), n);
    cudaStreamSynchronize(stream);
}
}
