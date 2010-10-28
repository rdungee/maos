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

#ifndef AOS_RECON_UTILS_H
#define AOS_RECON_UTILS_H
#include "maos.h"

void apply_L2(dcell **xout, spcell *L2, const dcell *xin, double alpha, int nthread);
void apply_invpsd(dcell **xout, INVPSD_T *extra,
		  const dcell *xin, double alpha);
void TTFR(dcell* x, const dcell *TTF, const dcell *PTTF);
void applyW(dcell *xin, const dsp *W0, 
	      const dmat *W1, const double *wt);
dcell* calcWmcc(const dcell *A, const dcell *B, const dsp *W0, 
		const dmat *W1, const dmat *wt);

spcell *act_slaving(LOC_T **aloc, spcell *HA, dmat *W1, dcell *NW);
#endif
