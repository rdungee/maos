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

#ifndef AOS_MATARRH_TYPE_H
#define AOS_MATARRH_TYPE_H
#include "numtype.h"
#include "array.h"
/**
   \file type.h Defines the math data types like dmat, cmat, dcell, ccell,
   dsp, csp data types.

   Don't use ulong for dimensions because subtracting a bigger ulong from a
   smaller ulong overflows.  */

typedef enum CEMBED{
    C_FULL,
    C_ABS2,
    C_REAL,
    C_ABS,
    C_LITERAL
}CEMBED;


typedef Mat<double> dmat;
typedef Mat<float> smat;
typedef Mat<dcomplex> cmat;
typedef Mat<fcomplex> zmat;
typedef Mat<long> lmat;
typedef Sparse<double> dsp;
typedef Sparse<float>  ssp;
typedef Sparse<dcomplex> csp;
typedef Sparse<fcomplex> zsp;


/**
   OPD or Amplitude map defined on square/rectangular grids. with equal spacing
   on x/y. 
*/
class map_t:public dmat{
public:
    double ox=0;      /**<Origin in x*/
    double oy=0;      /**<Origin in y*/
    double dx=0;      /**<Sampling along x*/
    double dy=0;      /**<Sampling along y*/
    double h=0;       /**<Heigh conjugation of this surface*/
    double vx=0;      /**Wind velocity. Useful for atmospheric grid*/
    double vy=0;      /**Wind velocity. Useful for atmospheric grid*/
    double iac=0;     /**<Inter-actuator coupling. >0: use cubic influence function*/
public:
    map_t(){}
    map_t(Int nxi, Int nyi, Real *pi, Real dxi, Real dyi)
	:dmat(nxi, nyi, pi, 0),ox(-nxi/2*dxi),oy(-nyi/2*dyi),dx(dxi),dy(dyi){}
    map_t(const map_t&in):dmat(in),ox(in.ox),oy(in.oy),dx(in.dx),dy(in.dy),h(in.h),vx(in.vx),vy(in.vy),iac(in.iac){}
};

/**
   Map with different x/y sampling. 
*/
class rmap_t:public dmat{
public:
    double ox=0;      /**<Origin in x*/
    double oy=0;      /**<Origin in y*/
    double dx=0;      /**<Sampling along x (first dimension)*/
    double dy=0;      /**<Sampling along y (second dimension)*/
    double txdeg=0;   /**<the x tilt angle in degree wrt beam (90 is prep), */
    double tydeg=0;   /**<the y tilt angle in degree wrt beam (90 is prep), */
    double ftel=0;    /**<Effective focal length of the telescope*/
    double fexit=0;   /**<The distance between the exit pupil and the focus*/
    double fsurf=0;   /**<The distance between the tilted surface (M3) and the focus*/
public:
    rmap_t(){}
    rmap_t(Int nxi, Int nyi, Real *pi, Real dxi, Real dyi)
	:dmat(nxi, nyi, pi, 0),ox(-nxi/2*dxi),oy(-nyi/2*dyi),dx(dxi),dy(dyi){}
    rmap_t(const rmap_t&in):dmat(in),ox(in.ox),oy(in.oy),dx(in.dx),dy(in.dy),txdeg(in.txdeg),tydeg(in.tydeg),ftel(in.ftel),fexit(in.fexit),fsurf(in.fsurf){}
};

/**
   Store starting x,y for each col
*/
struct locstatcol_t{
    double xstart; /**<starting x of this column*/
    double ystart; /**<starting y of this column*/
    long   pos;    /**<starting index of this column*/
};

/**
   Stores array of locstatcol_t

*/
struct locstat_t{
    Array<locstatcol_t> cols; /**<Information about each column*/
    double dx;          /**<Sampling of the grid along x*/
    double dy;          /**<Sampling of the grid along y*/
    double xmin;        /**<Minimum x*/
    double ymin;        /**<Minimum y*/
    long   ncol;        /**<Number of consecutive columns found*/
    long   nx;       /**<Size for embedding*/
    long   ny;       /**<Size for embedding*/
};
/**
   Struct for coordinates like plocs, xloc, aloc etc.
*/
class loc_t:protected dmat{
public:
    pcompat<double*>locx;  /**<x coordinates of each point*/
    pcompat<double*>locy;  /**<y coordinates of each point*/
    pcompat<long>   nloc;   /**<number of points*/
    locstat_t stat;/**<points to column statistics*/
    map_t map;    /**<point to the map used for identifying neihboring points.*/
    double dx;     /**<Sampling along x*/
    double dy;     /**<Sampling along y*/
    double ht=0;     /**<Conjugation height of the loc grid.*/
    double iac=0;    /**<Inter-actuator coupling. >0: use cubic influence function for ray tracing*/
    int npad=0;      /*padding when create map*/
public:
    loc_t(Int nloci, Real dxi, Real dyi):dmat(nloci, 2),dx(dxi),dy(dyi){
	id=M_LOC64;
	locx.SetP(P());
	locy.SetP(P()+Nx());
	nloc.SetP(Nx());
    }
};
/**
   low left point of each subaperture.
   
   don't change the leading 5 elements. so that pts_t can be used as loc_t.
*/
class pts_t:protected loc_t{
    uint32_t id;
    double *origx; /**<The x origin of each subaperture*/
    double *origy; /**<The y origin of each subaperture*/
    long nsa;      /**<number of subapertures.*/
    union{
	double dsa;    /**<side length of subaperture*/
	double dsax;   /**<side length of subaperture*/
    };
    double dsay;   /**<side length of subaperture*/
    double dummy1; /**<Place holder*/
    double dummy2;  /**<Place holder*/
    locstat_t *stat;/**<padding so that we can be casted to loc_t*/
    map_t *map;    /**<treat pts_t as loc_t and compute the MAP*/
    int npad;      /*padding when create map*/
    int nx;        /**<number of cols per subaperture*/
    int ny;        /**<number of rows per subaperture*/
    double dx;     /**<sampling of points in each subaperture*/
    double dy;     /**<sampling of points in each subaperture. dy=dx normally required.*/
    
};

typedef Cell<cmat> ccell;
typedef Cell<zmat> zcell;
typedef Cell<dmat> dcell;
typedef Cell<smat> scell;
typedef Cell<lmat> lcell;

typedef Cell<ccell> cccell;
typedef Cell<zcell> zccell;
typedef Cell<dcell> dccell;
typedef Cell<scell> sccell;
typedef Cell<lcell> lccell;

typedef Cell<cccell> ccccell;
typedef Cell<zccell> zcccell;
typedef Cell<dccell> dcccell;
typedef Cell<sccell> scccell;
typedef Cell<lccell> lcccell;

typedef Cell<map_t> mapcell;
typedef Cell<rmap_t> rmapcell;
typedef Cell<loc_t> loccell;
typedef Cell<mapcell> mapccell;
typedef Cell<rmapcell> rmapccell;
typedef Cell<loccell> locccell;

/*
typedef CELLARR(cmat*) ccell;
typedef CELLARR(zmat*) zcell;
typedef CELLARR(dmat*) dcell;
typedef CELLARR(smat*) scell;
typedef CELLARR(lmat*) lcell;
*/
typedef Cell<csp> cspcell;
typedef Cell<zsp> zspcell;
typedef Cell<dsp> dspcell;
typedef Cell<ssp> sspcell;

/*
typedef CELLARR(dsp*) dspcell;
typedef CELLARR(ssp*) sspcell;
typedef CELLARR(csp*) cspcell;
typedef CELLARR(zsp*) zspcell;
*/
/*
typedef CELLARR(ccell*) cccell;
typedef CELLARR(zcell*) zccell;
typedef CELLARR(dcell*) dccell;
typedef CELLARR(scell*) sccell;
typedef CELLARR(lcell*) iccell;

typedef CELLARR(cccell*) ccccell;
typedef CELLARR(zccell*) zcccell;
typedef CELLARR(dccell*) dcccell;
typedef CELLARR(sccell*) scccell;
typedef CELLARR(iccell*) icccell;

typedef CELLARR(map_t*) mapcell;
typedef CELLARR(rmap_t*) rmapcell;

typedef CELLARR(loc_t*) loccell;

typedef CELLARR(mapcell*) mapccell;
typedef CELLARR(rmapcell*) rmapccell;
typedef CELLARR(loccell*) locccell;


typedef struct cell{
    ARR(struct cell*);
    struct cell *m;
    }cell;*/
#undef ARR
#undef CELLARR
#undef MATARR

/*A method to simulate operator overloading for indexing arrys*/
#if DEBUG
INLINE void assert_1d(long i, long nx, long ny){
    if(i<0 || i>=nx*ny){
	error("%ld is out of range for (%ld,%ld) array\n", i, nx, ny);
    }
}
INLINE void assert_2d(long ix, long iy, long nx, long ny){
    if(ix<0 || ix>=nx || iy<0 || iy>=ny){
	error("(%ld,%ld) is out of range for (%ld,%ld) array\n", ix, iy, nx, ny);
    }
}
#define IND1(A,i) ((A)->p[assert_1d((i), (A)->nx, (A)->ny),(i)])
#define IND2(A,ix,iy) ((A)->p[assert_2d((ix), (iy), (A)->nx, (A)->ny),(ix)+(A)->nx*(iy)])
//#define PIND1(A,i) ((A)->p+(assert_1d(i, (A)->nx, (A)->ny),(i)))
//#define PIND2(A,ix,iy) ((A)->p+(assert_2d((ix), (iy), (A)->nx, (A)->ny),(ix)+(A)->nx*(iy)))
#else
#define IND1(A,i) ((A)->p[(i)])
#define IND2(A,ix,iy) ((A)->p[(ix)+(A)->nx*(iy)])
#endif
#define PIND1(A,i) ((A)->p+(i))
#define PIND2(A,ix,iy) ((A)->p+(ix)+(A)->nx*(iy))

#define IND0(A) error("Invalid use. Use IND(A,i) or IND(A,ix,iy)\n");
#define PIND0(A) error("Invalid use. Use PIND(A,i) or PIND(A,ix,iy)\n");
#define IND_GET(_0,_1,_2,_3,NAME,...) NAME
#define IND(...) IND_GET(_0,__VA_ARGS__,IND2,IND1,IND0,IND0)(__VA_ARGS__)
#define PIND(...) IND_GET(_0,__VA_ARGS__,PIND2,PIND1,PIND0,PIND0)(__VA_ARGS__)
#define PCOL(A,iy) &((A)->p[(iy)*(A)->nx])

#endif
