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

#ifndef AOS_TYPE_H
#define AOS_TYPE_H
#include "common.h"
#include "misc.h"
#if defined(DLONG)
typedef long spint; //we use 32 bit in sparse index.
#else
typedef int spint;
#endif
/**
   \file lib/type.h Defines the math data types like dmat, cmat, dcell, ccell,
   dsp, csp data types.

   Don't use ulong for dimensions because subtracting a bigger ulong from a
   smaller ulong overflows.  */


/*
  Separate definition of struct with typedef. Put definition of struct in a
private file and typedef in a public interface to hide the struct from the user
*/

/**
   a double matrix object contains 2-d array of double numbers
 */
typedef struct dmat{
    double *p;  /**< the pointer to allocated memory. */
    long nx;    /**< number of rows */
    long ny;    /**< number of columns */
    char *header;/**<Optionally contain a header*/
    long *nref; /**< reference count */
    struct mmap_t *mmap;/**< not NULL if mmaped. unmap the mmaped memory*/
}dmat;

/**
   a double complex matrix object contains 2-d arrays of
   double complex numbers. */
typedef struct cmat{
    dcomplex *p; /**< the pointer to allocated memory. */
    long nx;     /**< number of rows */
    long ny;     /**< number of columns */
    char *header;/**<Optionally contain a header*/
    long *nref;  /**< reference count */
    struct mmap_t *mmap;/**< not NULL if mmaped. unmap the mmaped memory*/
    struct fft_t *fft;
}cmat;

/**
   an 2-d block matrix of cmat.
 */
typedef struct ccell{
    cmat **p;   /**< Contains an array of pointers to cmat. */
    long nx;    /**< number of rows */
    long ny;    /**< number of columns */
    char *header;/**<Optionally contain a header*/
    struct mmap_t *mmap; /**<not NULL if mmaped. use unmap the mmaped memory*/
}ccell;
/**
   an 2-d block matrix of dmat.
 */
typedef struct dcell{
    dmat **p;   /**< Contains an array of pointers to dmat. */
    long nx;    /**< number of rows */
    long ny;    /**< number of columns */
    char *header;/**<Optionally contain a header*/
    struct mmap_t *mmap; /**< not NULL if mmaped. unmap the mmaped memory*/
}dcell;

typedef enum CEMBED{
    C_FULL,
    C_ABS2,
    C_REAL,
    C_ABS,
    C_LITERAL
}CEMBED;
/**
   a sparse array of double numbers stored in
   compressed column format, i.e. MATLAB format */
typedef struct dsp{
    long nzmax ;/**< maximum number of entries */
    long m ;	/**< number of rows */
    long n ;	/**< number of columns */
    spint *p ;   /**< column pointers (size n+1) or col indlces (size nzmax) when nz!=-1 */
    spint *i ;   /**< row indices, size nzmax */
    double *x ;	/**< numerical values, size nzmax */
    long nz ;   /**< number of entries in triplet matrix, -1 for compressed-col */
    long *nref; /**< reference counting like dmat */
}dsp;
/**
   an 2-d array of sparse.
 */
typedef struct spcell{
    dsp **p;    /**< Contains an array of pointers to dsp. */
    long nx;    /**< number of rows */
    long ny;    /**< number of columns */
}spcell;
/**
   a sparse array of double complex numbers stored in
   compressed column format */
typedef struct csp{
    long nzmax ;/**< maximum number of entries */
    long m ;	/**< number of rows */
    long n ;	/**< number of columns */
    spint *p ;   /**< column pointers (size n+1) or col indlces (size nzmax) */
    spint *i ;   /**< row indices, size nzmax */
    dcomplex *x;/**< numerical values, size nzmax */
    long nz ;   /**< # of entries in triplet matrix, -1 for compressed-col */
    long *nref; /**< reference counting like cmat*/
}csp;
/**
   an 2-d array of csp.
 */
typedef struct cspcell{
    csp **p;    /**< Contains an array of pointers to csp. */
    long nx;    /**< number of rows */
    long ny;    /**< number of columns */
}cspcell;


/**
   OPD or Amplitude map defined on square/rectangular grids. with equal spacing
   on x/y. The first 5 elements are all from dmat. Do not change them.
*/
typedef struct map_t{
    //The OPD, takes the same form of dmat so can be casted.
    double *p;  /**< the pointer to allocated memory. */
    long nx;    /**< number of rows */
    long ny;    /**< number of columns */
    char *header;/**<Optionally contain a header*/
    long *nref; /**< reference count */
    struct mmap_t *mmap;/**< not NULL if mmaped. unmap the mmaped memory*/
    double ox;      /**<Origin in x*/
    double oy;      /**<Origin in y*/
    double dx;      /**<Sampling along x and y*/
    double h;       /**<Heigh conjugation of this surface*/
    double vx;      /**Wind velocity. Useful for atmospheric grid*/
    double vy;      /**Wind velocity. Useful for atmospheric grid*/
#if USE_POSIX_SHM == 1
    //long shm;       /**Records the length of the memory mmaped of positive. -1 means it is part of shared shm.*/
#endif
} map_t;

/**
   Map with different x/y sampling.
*/
typedef struct rectmap_t{
    double *p;      /**<The OPD*/
    long nx;        /**<Number of points along x*/
    long ny;        /**<Number of points along y*/
    char *header;/**<Optionally contain a header*/
    long *nref; /**< reference count */
    long type;  /**< specify whether this is allocated or in shared memory. */
    double ox;      /**<Origin in x*/
    double oy;      /**<Origin in y*/
    double dx;      /**<Sampling along x (first dimension)*/
    double dy;      /**<Sampling along y (second dimension)*/
    double txdeg;   /**<the x tilt angle in degree wrt beam (90 is prep), */
    double tydeg;   /**<the y tilt angle in degree wrt beam (90 is prep), */
    double ftel;    /**<Effective focal length of the telescope*/
    double fexit;   /**<The distance between the exit pupil and the focus*/
    double fsurf;   /**<The distance between the tilted surface (M3) and the focus*/
}rectmap_t;

/**
   map of locs
*/
typedef struct locmap_t{
    long *p;       /**<The map, of size nx*ny*/
    double ox;     /**<Origin of the map along x*/
    double oy;     /**<Origin of the map along y*/
    int nx;        /**<Number of points along x*/
    int ny;        /**<Number of points along y*/
    int npad;      /**<Padding along the boundary. just for checking. no need in computation.*/
}locmap_t;

/**
   Struct for coordinates like plocs, xloc, aloc etc.
*/
typedef struct loc_t{
    double *locx;  /**< x coordinates of each point*/
    double *locy;  /**< y coordinates of each point*/
    long   nloc;   /**< number of points*/
    double dx;     /**< Sampling*/
    locmap_t *map; /**< point to the map used for identifying neihboring points.*/
    struct locstat_t *stat;/**<points to column statistics*/
}loc_t;

/**
   Store starting x,y for each col
*/
typedef struct locstatcol_t{
    double xstart; /**<starting x of this column*/
    double ystart; /**<starting y of this column*/
    long   pos;    /**<starting index of this column*/
}locstatcol_t;

/**
   Stores array of locstatcol_t
*/
typedef struct locstat_t{
    locstatcol_t *cols; /**<Information about each column*/
    double dx;          /**<Sampling of the grid*/
    double xmin;        /**<Minimum x*/
    double ymin;        /**<Minimum y*/
    long   ncol;        /**<Number of consecutive columns found*/
    long   nrow;        /**<Maximum number of rows*/
}locstat_t;
/**
   low left point of each subaperture.
   
   don't change the leading 5 elements. so that pts_t can be
   cast to loc_t

   2009-12-08: Bug found: I set the origx, origy to the lower
   left point in the subaperture instead of the true
   edge. This is good for raytracing. But I used it as the
   saloc orig for mkg, which is really wrong.  use
   powfs->saloc as the orig of subapertures. PTS just for ray
   tracing.
*/
typedef struct pts_t{
    double *origx; /**<The x origin of each subaperture*/
    double *origy; /**<The y origin of each subaperture*/
    long nsa;      /**<number of subapertures.*/
    double dsa;    /**<side length of subaperture*/
    locmap_t *map; /**<treat pts_t as loc_t and compute the MAP*/
    double *area;  /**<subaperture area, sum(amp^2)*/
    double dx;     /**<sampling of points in each subaperture*/
    int nx;        /**<number of col and row of points per subaperture*/
}pts_t;

/**
   Cell of anything
*/
typedef struct cell{
    dmat **p;   /**< Contains an array of pointers to dmat. */
    long nx;    /**< number of rows */
    long ny;    /**< number of columns */
    char *header;/**<Optionally contain a header*/
    struct mmap_t *mmap; /**< not NULL if mmaped. unmap the mmaped memory*/
    struct fft_t *fft; /**<Use a pair of dmat to do FFT in split storage.*/
}cell;

#define AOS_CMAT(A) c##A
#define AOS_CSPARSE(A) c##A
#define AOS_DMAT(A) d##A
#define AOS_SPARSE(A) A

#endif
