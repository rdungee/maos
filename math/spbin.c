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
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>

#include <fcntl.h>
#include <unistd.h>
#include "type.h"
#include "spbin.h"
#include "defs.h"
/*
   Function to write sparse matrix data into file pointed using a file
   pointer. Generally used by library developer.  We do not convert data during
   saving, but rather do the conversion during reading.
*/
void Y(spwritedata)(file_t *fp, const X(sp) *sp){
    header_t header={M_SPT, 0, 0, NULL};
    if(sp && sp->nzmax){
	header.nx=sp->m;
	header.ny=sp->n;
	header.str=sp->header;
    }
    write_header(&header, fp);
    if(sp && sp->nzmax){
	Y(spsort)((X(sp)*)sp);/*sort the matrix to have the right order */
	uint64_t nzmax;
	nzmax=sp->p[sp->n];/*don't use sp->nzmax, which maybe larger than actual */
	zfwritelarr(fp, 1, &nzmax);
	zfwrite(sp->p, sizeof(spint), sp->n+1, fp);
	zfwrite(sp->i, sizeof(spint), nzmax, fp);
	zfwrite(sp->x ,sizeof(T),nzmax,fp);  
    }
}
/**
   Function to read sparse matrix data from file pointer into memory. Used by
   library developer.
  */
X(sp) *Y(spreaddata)(file_t *fp, header_t *header){
    header_t header2;
    if(!header){
	header=&header2;
	read_header(header, fp);
    }
    long m=header->nx;
    long n=header->ny;
    uint64_t nzmax;
    X(sp) *out=NULL;
    if(m!=0 && n!=0){
	uint32_t magic2=0;
	switch(header->magic){
	case M_SPT64:
	    magic2=M_INT64;
	    break;
	case M_SPT32:
	    magic2=M_INT32;
	    break;
	default:
	    error("This is not a valid sparse matrix file. magic=%x\n", header->magic);
	}
	zfread(&nzmax,sizeof(uint64_t),1,fp);
	if(nzmax!=0){
	    out=Y(spnew)(m,n,nzmax);
	    out->header=header->str; header->str=NULL;
	    readspintdata(fp, magic2, out->p, n+1);
	    readspintdata(fp, magic2, out->i, nzmax);
	    zfread(out->x, sizeof(T), nzmax, fp);
	}
    }
    free(header->str);
    return out;
}

/**
   User callable function to write sparse matrix into file. 

   Usage: spwrite(A,"A.bin");
*/
void Y(spwrite)(const X(sp) *sp, const char *format,...){
    format2fn;
    /*write the sparse matrix to file to later load from matlab */
    file_t *fp=zfopen(fn,"wb");
    Y(spwritedata)(fp, sp);
    /*don't worry about the warning of 0x401ee45 in valgrind. That is the IO  */
    zfclose(fp);
}
/**
   User callable function to write cell array of sparse matrix into file. 

   Usage: spcellwrite(A,"A.bin"); */
void Y(spcellwrite)(const Y(spcell) *spc, const char *format,...){
    format2fn;
    file_t *fp=zfopen(fn,"wb");
    header_t header={MCC_ANY, 0, 0, NULL};
    if(spc){
	header.nx=spc->nx;
	header.ny=spc->ny;
	header.str=spc->header;
    }
    write_header(&header, fp);
    if(spc){
	for(unsigned long iy=0; iy<spc->ny; iy++){
	    for(unsigned long ix=0; ix<spc->nx; ix++){
		Y(spwritedata)(fp, spc->p[ix+iy*spc->nx]);
	    }
	}
    }
    zfclose(fp);
}
/**
   User callable function to read sparse metrix from file. 

   Usage: A=spread("A.bin");*/
X(sp)* Y(spread)(const char *format,...){
    format2fn;
    file_t *fp=zfopen(fn,"rb");
    X(sp) *out=Y(spreaddata)(fp, 0);
    zfeof(fp);
    zfclose(fp);
    return out;
}
/**
   User callable function to read cell array of sparse matrix
   from file. Usage: A=spcellread("A.bin");
 */
Y(spcell) *Y(spcellread)(const char *format,...){
    format2fn;
    file_t *fp=zfopen(fn,"rb");
    header_t header;
    read_header(&header, fp);
    header_t *headerc=NULL;
    long nx, ny;
    if(!iscell(header.magic)){
	headerc=&header;
	warning("%s is not a sparse cell file. want %d, got %u\n", fn, MCC_ANY, header.magic);
	nx=1; ny=1;
    }else{
	free(header.str);
	nx=header.nx;
	ny=header.ny;
    }
    Y(spcell) *out;
    out=Y(spcellnew)(nx,ny);
    out->header=header.str; header.str=NULL;
    for(unsigned long ix=0; ix<nx*ny; ix++){
	out->p[ix]=Y(spreaddata)(fp, headerc);
    }
    zfeof(fp);
    zfclose(fp);
    return out;
}