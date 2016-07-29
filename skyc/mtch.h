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

#ifndef SKYC_MTCH_H
#define SKYC_MTCH_H
#include "skyc.h"
#include "types.h"
void psf2i0gxgy(dmat *i0, dmat *gx, dmat *gy, dmat *psf, DTF_S *dtf);
void genmtch(dcell **mtche, dmat **sanea,
	     dcell *i0, dcell *gx, dcell *gy, double pixtheta, 
	     double rne, double bkgrnd, int cr);
#endif
