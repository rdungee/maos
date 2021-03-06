/*
  Copyright 2009-2018 Lianqi Wang <lianqiw-at-tmt-dot-org>
  
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

#include "common.h"
#ifndef AOS_SIM_UTILS_H
#define AOS_SIM_UTILS_H
void atm2xloc(dcell **opdx, const SIM_T *simu);
void sim_update_etf(SIM_T *simu);
void seeding(SIM_T *simu);
SIM_T* init_simu(const PARMS_T *parms,POWFS_T *powfs, APER_T *aper,RECON_T *recon,int iseed);
void free_simu(SIM_T *simu);
void print_progress(const SIM_T *simu);
void save_skyc(POWFS_T *powfs, RECON_T *recon, const PARMS_T *parms);
void genatm(SIM_T *simu);
void setup_recon_HXW_predict(SIM_T *simu);
#endif
