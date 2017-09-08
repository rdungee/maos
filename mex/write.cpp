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
#ifdef __INTEL_COMPILER
#undef _GNU_SOURCE /*avoid compiling problem*/
#endif



/*compile with 
mex write.c -largeArrayDims
usage write(filename, data);
*/
#include "io.h"
static char *mx2str(const mxArray *header){
    char *str=NULL;
    if(header && mxIsChar(header)){
	int nheader=mxGetM(header)*mxGetN(header)+1;
	str=(char*)malloc(nheader);
	mxGetString(header, str, nheader);
	//convert matlab \n (2 chars) to C \n (1 char with a space)
	const char *str2=str+strlen(str)-1;
	for(char *h=str; h<str2; h++){
	    if(h[0]=='\\' && h[1]=='n'){
		h[0]=' ';
		h[1]='\n';
	    }
	}
    }
    return str;
}
static void writedata(file_t *fp, int type, const mxArray *arr, const mxArray *header){
    char *str=mx2str(header);
    header_t header2={0,0,0,0};
   
    //uint64_t m,n;
    if(!arr){
	header2.ndim=0;
	header2.dims=0;
    }else{
	header2.ndim=mxGetNumberOfDimensions(arr);
	header2.dims=(mwSize*)mxGetDimensions(arr);
	if(header2.ndim==0){
	    header2.ntot=0;
	}else{
	    header2.ntot=1;
	    for(uint64_t id=0; id<header2.ndim; id++){
		header2.ntot*=header2.dims[id];
	    }
	}
	if(!header2.ntot){
	    arr=NULL;
	}
    }
    if(arr && mxIsCell(arr)){
	header2.magic=MCC_ANY;
	int issparse=0;
	mxArray *in;
	int type2;
	if(mxGetNumberOfElements(arr)==0){
	    in=NULL;
	    issparse=0;
	    type2=M_DBL;
	}else{
	    for(size_t ix=0; ix<mxGetNumberOfElements(arr); ix++){
		in=mxGetCell(arr,ix);
		if(in){
		    if(mxIsSparse(in)){
			issparse=1;
			if(mxIsComplex(in)){
			    type2=MAT_CSP;
			}else{
			    type2=MAT_SP;
			}
		    }else{
			issparse=0;
			if(mxIsComplex(in)){
			    if(mxIsSingle(in)){
				type2=M_ZMP;
			    }else{
				type2=M_CMP;
			    }
			}else{
			    if(mxIsSingle(in)){
				type2=M_FLT;
			    }else{
				type2=M_DBL;
			    }
			}
		    }
		    break;
		}
	    }
	}
	if(!in){//all cell empty
	    issparse=0;
	    type2=M_DBL;
	}
	//don't write global header.
	write_header(&header2, fp);
	for(size_t ix=0; ix<mxGetNumberOfElements(arr); ix++){
	    in=mxGetCell(arr, ix);
	    if(in && !mxIsEmpty(in) && mxIsSparse(in) !=issparse)
		error("can only save cell array of all sparse or all dense");
	    if(header && mxIsCell(header)){
		writedata(fp, type2, in, mxGetCell(header, ix));
	    }else{
		writedata(fp, type2, in, header);
	    }
	}
    }else{/*not cell.*/
	if(type == MAT_SP || ((arr) && mxIsSparse(arr))){
	    if(sizeof(mwIndex)==4){
		header2.magic=M_SP32;
	    }else{
		header2.magic=M_SP64;
	    }
	    write_header(&header2, fp);
	    if(header2.ntot){
		mwIndex *Jc=mxGetJc(arr);
		long n=header2.dims[1];
		long nzmax=Jc[n];
		zfwrite(&nzmax, sizeof(double), 1, fp);
		zfwrite(mxGetJc(arr), sizeof(mwIndex), n+1, fp);
		zfwrite(mxGetIr(arr), sizeof(mwIndex), nzmax, fp);
		zfwrite(mxGetPr(arr), sizeof(double), nzmax, fp);
	    }
	}else if(type == MAT_CSP || ((arr) && mxIsSparse(arr))){
	    if(sizeof(mwIndex)==4){
		header2.magic=M_CSP32;
	    }else{
		header2.magic=M_CSP64;
	    }
	    write_header(&header2, fp);
	    if(header2.ntot){
		mwIndex *Jc=mxGetJc(arr);
		long n=header2.dims[1];
		long nzmax=Jc[n];
		zfwrite(&nzmax, sizeof(double), 1, fp);
		zfwrite(mxGetJc(arr), sizeof(mwIndex), n+1, fp);
		zfwrite(mxGetIr(arr), sizeof(mwIndex), nzmax, fp);
		zfwrite_dcomplex(mxGetPr(arr),mxGetPi(arr),nzmax,fp);
	    }
	}else if(type == M_DBL || ((arr) && mxIsDouble(arr) && !mxIsComplex(arr))){
	    header2.magic=M_DBL;
	    write_header(&header2, fp);
	    if(header2.ntot){
		zfwrite(mxGetPr(arr), sizeof(double), header2.ntot, fp);
	    }  
	}else if(type == M_FLT || ((arr) && mxIsSingle(arr) && !mxIsComplex(arr))){
	    header2.magic=M_FLT;
	    write_header(&header2, fp);
	    if(header2.ntot){
		zfwrite(mxGetPr(arr), sizeof(float), header2.ntot, fp);
	    }
	}else if(type == M_CMP || ((arr)&& mxIsDouble(arr) && mxIsComplex(arr))){
	    header2.magic=M_CMP;
	    write_header(&header2, fp);
	    if(header2.ntot){
		zfwrite_dcomplex(mxGetPr(arr),mxGetPi(arr), header2.ntot, fp);
	    }
	}else if(type == M_ZMP || ((arr)&& mxIsSingle(arr) && mxIsComplex(arr))){
	    header2.magic=M_ZMP;
	    write_header(&header2, fp);
	    if(header2.ntot){
		zfwrite_fcomplex((float*)mxGetPr(arr),(float*)mxGetPi(arr), header2.ntot, fp);
	    }
	}else{
	    error("Unrecognized data type");
	}
    }
    free(str);
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
    (void)nlhs;
    (void)plhs;
    file_t *fp;
    int ifn=1;
    const mxArray *header=NULL;
    if(nrhs==2){/*data and file*/
	ifn=1;
    }else if(nrhs==3){/*data, header, and file*/
	ifn=2;
	header=prhs[1];
    }else{
	mexErrMsgTxt("Usage: write(var,'file') or write(var, header, 'file')\n");
    }
    int nlen=mxGetM(prhs[ifn])*mxGetN(prhs[ifn])+1;
    char *fn=(char*)malloc(nlen);
    mxGetString(prhs[ifn],fn,nlen);
    fp=zfopen(fn,"wb");
    if(!fp){
	mexErrMsgTxt("Error writing file.\n");
	return;
    }
    free(fn);
    writedata(fp, 0, prhs[0], header);
    zfclose(fp);
}
