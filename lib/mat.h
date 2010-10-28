/*
  Copyright 2009,2010 Lianqi Wang <lianqiw@gmail.com> <lianqiw@tmt.org>
  
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

#ifndef AOS_MAT_H
#define AOS_MAT_H
#include "random.h"
#include "type.h"
void dembed(dmat *restrict A, dmat *restrict B, const double theta);

#define AOS_MAT_DEF(X,Y,T)\
X(mat) *X(new_ref)(T *p, long nx, long ny) CHECK_UNUSED_RESULT;\
X(mat) *X(new_data)(T *p, long nx, long ny) CHECK_UNUSED_RESULT;\
X(mat) *X(new)(long nx, long ny) CHECK_UNUSED_RESULT;\
void X(free_keepdata)(X(mat) *A);\
void X(free_do)(X(mat) *A, int keepdata);\
void X(resize)(X(mat) *A, long nx, long ny);\
X(mat) *X(ref)(X(mat) *in) CHECK_UNUSED_RESULT;\
X(mat) *X(ref_reshape)(X(mat) *in, int nx, int ny) CHECK_UNUSED_RESULT;\
X(mat) *X(refcols)(X(mat) *in, long icol, long ncol) CHECK_UNUSED_RESULT;\
X(mat) *X(sub)(const X(mat) *in, long sx, long nx, long sy, long ny) CHECK_UNUSED_RESULT;\
X(mat) *X(cat)(const X(mat) *in1, const X(mat) *in2, int dim) CHECK_UNUSED_RESULT;\
void X(arrfree)(X(mat) **As, int n);\
X(mat) *X(dup)(const X(mat) *in) CHECK_UNUSED_RESULT;\
void X(cp)(X(mat) **out0, const X(mat) *in);\
X(mat) *X(trans)(const X(mat) *A) CHECK_UNUSED_RESULT;\
void X(set)(X(mat) *A, const T val);\
void X(zero)(X(mat) *out);\
double X(max)(const X(mat) *A) CHECK_UNUSED_RESULT;\
double X(min)(const X(mat) *A) CHECK_UNUSED_RESULT;\
double X(norm2)(const X(mat) *in) CHECK_UNUSED_RESULT;\
void X(randu)(X(mat) *A, const T mean, struct_rand *rstat);\
void X(randn)(X(mat) *A, const T sigma, struct_rand *rstat);\
void X(show)(const X(mat) *A, const char *format,...) CHECK_ARG(2);\
void X(scale)(X(mat) *A, T w);\
T X(sum)(const X(mat) *A) CHECK_UNUSED_RESULT;\
void X(add)(X(mat) **B0, T bc,const X(mat) *A, const T ac);\
T X(inn)(const X(mat)*A, const X(mat) *B);			\
T X(wdot)(const T *a, const X(mat) *w, const T *b) CHECK_UNUSED_RESULT;\
T X(wdot2)(const T *a, const X(mat) *w, const T *b) CHECK_UNUSED_RESULT;\
T X(wdot3)(const T *a, const X(mat) *w, const T *b) CHECK_UNUSED_RESULT;\
void X(cwm)(X(mat) *B, const X(mat) *A);\
void X(mulvec)(T *restrict y, const X(mat) * restrict A, const T *restrict x, const T alpha);\
void X(mm)(X(mat)**C0, const X(mat) *A, const X(mat) *B, const char trans[2], const T alpha);\
void X(invspd_inplace)(X(mat) *A);\
X(mat)* X(invspd)(const X(mat) *A) CHECK_UNUSED_RESULT;\
void X(inv_inplace)(X(mat)*A);\
X(mat)* X(inv)(const X(mat) *A) CHECK_UNUSED_RESULT;\
X(mat) *X(mcc)(const X(mat) *A, const X(mat) *wt) CHECK_UNUSED_RESULT;\
X(mat) *X(imcc)(const X(mat) *A, const X(mat) *wt) CHECK_UNUSED_RESULT;\
X(mat) *X(tmcc)(const X(mat) *A, const X(mat) *wt) CHECK_UNUSED_RESULT;\
X(mat) *X(pinv)(const X(mat) *A, const X(mat) *wt, const X(sp) *Wsp) CHECK_UNUSED_RESULT;\
T X(diff)(const X(mat) *A, const X(mat) *B) CHECK_UNUSED_RESULT;\
void X(circle)(X(mat) *A, double cx, double cy, double r, T val);\
void X(circle_symbolic)(X(mat) *A, double cx, double cy, double r);\
void X(fftshift)(X(mat) *A);\
void X(cpcorner2center)(X(mat) *A, const X(mat)*B);\
void X(shift)(X(mat) **B0, const X(mat) *A, int sx, int sy);\
void X(rotvec)(X(mat) *A, const double theta);\
void X(rotvecnn)(X(mat) **B0, const X(mat) *A, double theta);\
void X(mulvec3)(T *y, const X(mat) *A, const T *x);\
void X(cog)(double *grad,const X(mat) *i0,double offsetx, double offsety, double thres, double bkgrnd);\
void X(shift2center)(X(mat) *A, double offsetx, double offsety);\
int X(clip)(X(mat) *A, double min, double max);\
void X(gramschmidt)(X(mat) *Mod, double *amp);	\
void X(muldiag)(X(mat) *A, X(mat) *s);\
void X(cwpow)(X(mat) *A, double power);\
void X(svd)(dmat **Sdiag, X(mat) **U, X(mat) **VT, const X(mat) *A);\
void X(svd_pow)(X(mat) *A, double power);\
void X(addI)(X(mat) *A, T val);\
void X(tikcr)(X(mat) *A, T thres);\
void X(mulsp)(X(mat) **yout, const X(mat) *x, const X(sp) *A, const T alpha);\
X(mat)* X(logspace)(double emin, double emax, long n) CHECK_UNUSED_RESULT;\
X(mat)* X(linspace)(double emin, double emax, long n) CHECK_UNUSED_RESULT;\
X(mat)* X(interp1)(dmat *xin, dmat *yin, dmat *xnew) CHECK_UNUSED_RESULT;\
X(mat)* X(interp1log)(dmat *xin, dmat *yin, dmat *xnew) CHECK_UNUSED_RESULT;\
void X(histfill)(dmat **out, const X(mat)* A, double center, double spacing, int n);\
X(mat) *X(spline_prep)(X(mat) *x, X(mat) *y);\
X(mat)* X(spline_eval)(X(mat) *coeff, X(mat)* x, X(mat)*xnew);\
X(mat)* X(spline)(X(mat) *x,X(mat) *y,X(mat) *xnew);\
X(cell)* X(bspline_prep)(X(mat)*x, X(mat)*y, X(mat) *z);\
X(mat) *X(bspline_eval)(X(cell)*coeff, X(mat) *x, X(mat) *y, X(mat) *xnew, X(mat) *ynew);

#endif
