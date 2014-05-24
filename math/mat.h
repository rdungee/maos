/*
  Copyright 2009-2013 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#ifndef AOS_LIB_MAT_H
#define AOS_LIB_MAT_H
#ifndef AOS_LIB_MATH_H
#error "Don't include this file directly"
#endif

#include "random.h"
#ifdef USE_SINGLE
#define cabs2(A)     (powf(crealf(A),2)+powf(cimagf(A),2))
#else
#define cabs2(A)     (pow(creal(A),2)+pow(cimag(A),2))
#endif
#ifdef __clang__
/*clang doesnot accept the gcc version in C++ mode*/
#define PALL(T,A,pp) typedef T pp##_ptr[(A)->nx]; pp##_ptr *pp=(pp##_ptr*)(A)->p
#else 
/*clang version caused <anonymous> must be uninitailized error in cuda code.*/
#define PALL(T,A,pp) T (*pp)[(A)->nx]=(T(*)[(A)->nx])(A)->p
#endif
#define PDMAT(M,P)   PALL(double,M,P)
#define PDCELL(M,P)  PALL(dmat*,M,P)
#define dfree(A)     ({dfree_do((A),0);(A)=NULL;})
#define dcp2(A,B)    memcpy(A->p,B->p,sizeof(double)*A->nx*A->ny)
#define dcellfree(A) ({dcellfree_do(A);A=NULL;})
#define dcellfreearr(A,n) ({for(int in=0; A&&in<n; in++){dcellfree(A[in]);};free(A);A=NULL;})
#define dzero(A)     if(A) memset((A)->p, 0, (A)->nx*(A)->ny*sizeof(double))
#define dhash(A,key) hashlittle(A->p, A->nx*A->ny*sizeof(double), key)

#define PSMAT(M,P)   PALL(float,M,P)
#define PSCELL(M,P)  PALL(smat*,M,P)
#define sfree(A)     ({sfree_do((A),0);(A)=NULL;})
#define scp2(A,B)    memcpy(A->p,B->p,sizeof(float)*A->nx*A->ny)
#define scellfree(A) ({scellfree_do(A);A=NULL;})
#define scellfreearr(A,n) ({for(int in=0; A&&in<n; in++){scellfree(A[in]);};free(A);A=NULL;})
#define szero(A) if(A) memset((A)->p, 0, (A)->nx*(A)->ny*sizeof(float))
#define shash(A,key) hashlittle(A->p, A->nx*A->ny*sizeof(float), key)

#define PCMAT(M,P)   PALL(dcomplex,M,P)
#define PCCELL(M,P)  PALL(cmat*,M,P)
#define cfree(A)     ({cfree_do(A,0);A=NULL;})
#define ccellfree(A) ({ccellfree_do(A);A=NULL;})
#define ccellfreearr(A,n) ({for(int in=0; A&&in<n; in++){ccellfree(A[in]);};free(A);A=NULL;})
#define czero(A)     if(A) memset((A)->p, 0, (A)->nx*(A)->ny*sizeof(dcomplex))
#define chash(A,key) hashlittle(A->p, A->nx*A->ny*sizeof(dcomplex), key)

#define PZMAT(M,P)   PALL(fcomplex,M,P)
#define PZCELL(M,P)  PALL(zmat*,M,P) 
#define zfree(A)     ({zfree_do(A,0);A=NULL;})
#define zcellfree(A) ({zcellfree_do(A);A=NULL;})
#define zzero(A)     if(A) memset((A)->p, 0, (A)->nx*(A)->ny*sizeof(fcomplex))
#define zhash(A,key) hashlittle(A->p, A->nx*A->ny*sizeof(fcomplex), key)

#define AOS_MAT_DEF(X,XR,Y,T,R)					\
X(mat) *X(new_ref)(long nx, long ny, T *p) CHECK_UNUSED_RESULT; \
X(mat) *X(new_data)(long nx, long ny, T *p) CHECK_UNUSED_RESULT; \
X(mat) *X(new)(long nx, long ny) CHECK_UNUSED_RESULT;\
void X(init)(X(mat)**A, long nx, long ny) ; \
void X(free_keepdata)(X(mat) *A);\
void X(free_do)(X(mat) *A, int keepdata);\
void X(resize)(X(mat) *A, long nx, long ny);\
X(mat) *X(ref)(const X(mat) *in) CHECK_UNUSED_RESULT;\
X(mat) *X(ref_reshape)(X(mat) *in, long nx, long ny) CHECK_UNUSED_RESULT;\
X(mat) *X(refcols)(X(mat) *in, long icol, long ncol) CHECK_UNUSED_RESULT;\
X(mat) *X(sub)(const X(mat) *in, long sx, long nx, long sy, long ny) CHECK_UNUSED_RESULT;\
int X(isnan)(const X(mat)*A);\
X(mat) *X(cat)(const X(mat) *in1, const X(mat) *in2, int dim) CHECK_UNUSED_RESULT;\
void X(arrfree)(X(mat) **As, int n);\
X(mat) *X(dup)(const X(mat) *in) CHECK_UNUSED_RESULT;\
void X(cp)(X(mat) **out0, const X(mat) *in);\
X(mat) *X(trans)(const X(mat) *A) CHECK_UNUSED_RESULT;\
void X(set)(X(mat) *A, const T val);\
void X(maxmin)(const T *restrict p, long N, R *max, R *min);\
R X(max)(const X(mat) *A) CHECK_UNUSED_RESULT;\
R X(maxabs)(const X(mat) *A) CHECK_UNUSED_RESULT;\
R X(min)(const X(mat) *A) CHECK_UNUSED_RESULT;\
R X(norm2)(const X(mat) *in) CHECK_UNUSED_RESULT;\
void X(randu)(X(mat) *A, const T mean, rand_t *rstat);\
void X(randn)(X(mat) *A, const T sigma, rand_t *rstat);\
void X(show)(const X(mat) *A, const char *format,...) CHECK_ARG(2);\
void X(scale)(X(mat) *A, R w);\
T X(sum)(const X(mat) *A) CHECK_UNUSED_RESULT;\
T X(trace)(const X(mat) *A) CHECK_UNUSED_RESULT;\
void X(add)(X(mat) **B0, T bc,const X(mat) *A, const T ac);\
void X(adds)(X(mat*)A, const T ac);\
T X(inn)(const X(mat)*A, const X(mat) *B);			\
T X(wdot)(const T *a, const X(mat) *w, const T *b) CHECK_UNUSED_RESULT;\
T X(wdot2)(const T *a, const X(mat) *w, const T *b) CHECK_UNUSED_RESULT;\
T X(wdot3)(const T *a, const X(mat) *w, const T *b) CHECK_UNUSED_RESULT;\
void X(cwm)(X(mat) *B, const X(mat) *A);\
void X(cwm3)(X(mat) *restrict A, const X(mat) *restrict B, const X(mat) *restrict C);\
void X(cwmcol)(X(mat) *restrict A, const X(mat) *restrict B);\
void X(cwm3col)(X(mat) *restrict A,const X(mat) *restrict W,const X(mat) *restrict B);\
void X(cwmrow)(X(mat) *restrict A, const X(mat) *restrict B);\
void X(cwmcol2)(X(mat) *restrict A, \
		const T *restrict B1, const R wt1,	\
		const T *restrict B2, const R wt2);	\
void X(cwmrow2)(X(mat) *restrict A,			\
		const T *restrict B1, const R wt1,	\
		const T *restrict B2, const R wt2);	\
void X(cwdiv)(X(mat) *B, const X(mat) *A, T value);			\
void X(mulvec)(T *restrict y, const X(mat) * restrict A, const T *restrict x, const T alpha);\
void X(mm)(X(mat)**C0, const T beta, const X(mat) *A, const X(mat) *B, const char trans[2], const T alpha); \
void X(invspd_inplace)(X(mat) *A);\
X(mat)* X(invspd)(const X(mat) *A) CHECK_UNUSED_RESULT;\
void X(inv_inplace)(X(mat)*A);\
X(mat)* X(inv)(const X(mat) *A) CHECK_UNUSED_RESULT;\
X(mat) *X(chol)(const X(mat) *A) CHECK_UNUSED_RESULT; \
X(mat) *X(mcc)(const X(mat) *A, const X(mat) *wt) CHECK_UNUSED_RESULT;\
X(mat) *X(imcc)(const X(mat) *A, const X(mat) *wt) CHECK_UNUSED_RESULT;\
X(mat) *X(tmcc)(const X(mat) *A, const X(mat) *wt) CHECK_UNUSED_RESULT;\
X(mat) *X(pinv)(const X(mat) *A, const X(mat) *wt, const X(sp) *Wsp) CHECK_UNUSED_RESULT;\
T X(diff)(const X(mat) *A, const X(mat) *B) CHECK_UNUSED_RESULT;\
void X(circle)(X(mat) *A, R cx, R cy, R dx, R dy, R r, T val); \
void X(circle_mul)(X(mat) *A, R cx, R cy, R dx, R dy, R r, T val);\
void X(circle_symbolic)(X(mat) *A, R cx, R cy, R dx, R dy, R r);\
void X(fftshift)(X(mat) *A);\
void X(cpcorner2center)(X(mat) *A, const X(mat)*B);\
void X(shift)(X(mat) **B0, const X(mat) *A, int sx, int sy);\
void X(rotvec)(X(mat) *A, const R theta);\
void X(rotvect)(X(mat) *A, const R theta);\
void X(rotvecnn)(X(mat) **B0, const X(mat) *A, R theta);\
void X(mulvec3)(T *y, const X(mat) *A, const T *x);\
void X(cog)(R *grad,const X(mat) *i0,R offsetx, R offsety, R thres, R bkgrnd);\
void X(shift2center)(X(mat) *A, R offsetx, R offsety);\
int X(clip)(X(mat) *A, R min, R max);\
void X(gramschmidt)(X(mat) *Mod, R *amp);	\
void X(muldiag)(X(mat) *A, const X(mat) *s);\
void X(muldiag2)(X(mat) *A, const X(mat) *s);\
void X(cwpow)(X(mat) *A, R power);\
void X(cwexp)(X(mat) *A, R alpha);\
void X(cwpow_thres)(X(mat) *A, R power, R thres);		\
void X(svd)(X(mat) **U, XR(mat) **Sdiag, X(mat) **VT, const X(mat) *A); \
void X(svd_pow)(X(mat) *A, R power, R thres);  \
void X(expm)(X(mat) **out, R alpha, X(mat) *A, R beta); \
void X(polyval)(X(mat) *A, XR(mat)*p);\
void X(addI)(X(mat) *A, T val);\
void X(tikcr)(X(mat) *A, T thres);\
void X(mulsp)(X(mat) **yout, const X(mat) *x, const X(sp) *A, const T alpha);\
X(mat)* X(logspace)(R emin, R emax, long n) CHECK_UNUSED_RESULT;\
X(mat)* X(linspace)(R min, R dx, long n) CHECK_UNUSED_RESULT;\
X(mat)* X(interp1)(const X(mat) *xin, const X(mat) *yin, const X(mat) *xnew) CHECK_UNUSED_RESULT;\
X(mat)* X(interp1_2)(const X(mat) *xyin, const X(mat) *xnew) CHECK_UNUSED_RESULT; \
X(mat)* X(interp1linear)(const X(mat) *xin, const X(mat) *yin, const X(mat) *xnew) CHECK_UNUSED_RESULT;\
X(mat)* X(interp1log)(const X(mat) *xin, const X(mat) *yin, const X(mat) *xnew) CHECK_UNUSED_RESULT;\
void X(blend)(X(mat) *restrict A, X(mat) *restrict B, int overlap);\
void X(histfill)(X(mat) **out, const X(mat)* A, R center, R spacing, int n);\
X(mat) *X(spline_prep)(X(mat) *x, X(mat) *y);\
X(mat)* X(spline_eval)(X(mat) *coeff, X(mat)* x, X(mat)*xnew);\
X(mat)* X(spline)(X(mat) *x,X(mat) *y,X(mat) *xnew);\
void X(cwlog10)(X(mat) *A);\
void X(cwlog)(X(mat) *A);\
void X(embed)(X(mat) *restrict A, const X(mat) *restrict B, const R theta); \
long X(fwhm)(X(mat) *A);\
void X(sort)(X(mat) *A, int ascend);\
X(mat) *X(enc)(X(mat) *A, X(mat) *dvec, int type, int nthread);\
typedef T (*X(minsearch_fun))(T *x, void *info);			\
int X(minsearch)(T *x, T *scale, int nmod, T ftol, X(minsearch_fun) fun, void *info);\
void X(bessik)(T x, T xnu, T *ri, T *rk, T *rip, T *rkp);

/*The following are only useful for cmat */
#define AOS_CMAT_DEF(X,XR,Y,T,R)					\
void X(cwmc)(X(mat) *restrict A, const X(mat) *restrict B, const R alpha); \
void X(cwmd)(X(mat) *restrict A, const XR(mat) *restrict B, const R alpha); \
void X(embed_wvf)(X(mat) *restrict A, const R *opd, const R *amp,	\
		  const int nopdx, const int nopdy,			\
		  const R wvl, const R theta);				\
void X(embedc)(X(mat) *restrict A, const X(mat) *restrict B, const R theta,CEMBED flag); \
void X(embedd)(X(mat) *restrict A, XR(mat) *restrict B, const R theta);	\
void X(embedscaleout)(X(mat) *restrict A, const X(mat) * in,		\
		      R xoutscale,R youtscale,				\
		      const R theta, CEMBED flag);			\
void X(cpcorner)(X(mat) *A, const X(mat) *restrict B, CEMBED flag);	\
void X(abstoreal)(X(mat) *A);						\
void X(abs2toreal)(X(mat) *A);						\
void X(cpd)(X(mat)**restrict A, const XR(mat) *restrict B);		\
void X(real2d)(XR(mat)**restrict A0, R alpha,const X(mat) *restrict B, R beta);	\
void X(abs22d)(XR(mat)**restrict A0, R alpha,const X(mat) *restrict B, R beta);	\
void X(cp)(X(mat)**restrict A0, const X(mat) *restrict B);		\
void X(tilt2)(X(mat) *otf, X(mat) *otfin, R sx, R sy, int pinct);	\
void X(tilt)(X(mat) *otf, R sx, R sy, int pinct);			\
void X(cogreal)(R *grad,const X(mat) *i0,R offsetx,			\
		R offsety,R thres, R bkgrnd);				\
void X(cogabs)(R *grad,const X(mat) *i0,R offsetx,			\
	       R offsety,R thres, R bkgrnd);				\
void X(inv_inplace)(X(mat)*A);						\
void X(invspd_inplace)(X(mat) *A);					\
void X(mulvec)(T *restrict y, const X(mat) * restrict A,		\
	       const T *restrict x, const T alpha);


#endif