#ifdef __INTEL_COMPILER
#undef _GNU_SOURCE /*avoid compiling problem*/
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
/*compile with 
mex write.c -largeArrayDims
usage write(filename, data);
*/
#include "io.h"
static void write_header2(file_t *fp, const mxArray *header){
    if(mxIsChar(header)){
	int nheader=mxGetM(header)*mxGetN(header)+1;
	char *header2=malloc(nheader);
	mxGetString(header, header2, nheader);
	write_header(header2,fp);
	free(header2);
    }
}
static void writedata(file_t *fp, int type, const mxArray *arr, const mxArray *header){
    uint32_t magic;
    uint64_t m,n;
    if(!arr){
	m=0; n=0;
    }else{
	m=mxGetM(arr);
	n=mxGetN(arr);
    }
    if(arr && mxIsCell(arr)){
	int issparse=0;
	mxArray *in;
	int type2;
	long ix;
	if(mxGetNumberOfElements(arr)==0){
	    in=NULL;
	    issparse=0;
	    magic=MC_DBL;
	    type2=M_DBL;
	}else{
	    in=mxGetCell(arr,0);
	    if(mxIsSparse(in)){
		issparse=1;
		if(mxIsComplex(in)){
		    magic=MC_CSP;
		    type2=MAT_CSP;
		}else{
		    magic=MC_SP;
		    type2=MAT_SP;
		}
	    }else{
		issparse=0;
		if(mxIsComplex(in)){
		    magic=MC_CMP;
		    type2=M_CMP;
		}else{
		    magic=MC_DBL;
		    type2=M_DBL;
		}
	    }
	}
	if(header && !mxIsCell(header)){
	    write_header2(fp, header);/*header for the cell.*/
	}
	zfwrite(&magic, sizeof(uint32_t), 1, fp);
	zfwrite(&m, sizeof(uint64_t), 1, fp);
	zfwrite(&n, sizeof(uint64_t), 1, fp);
	for(ix=0; ix<mxGetNumberOfElements(arr); ix++){
	    in=mxGetCell(arr, ix);
	    if(in && !mxIsEmpty(in) && mxIsSparse(in) !=issparse)
		error("can only save cell array of all sparse or all dense");
	    if(header && mxIsCell(header)){
		writedata(fp, type2, in, mxGetCell(header, ix));
	    }else{
		writedata(fp, type2, in, NULL);
	    }
	}
    }else{/*not cell.*/
	if(header){
	    write_header2(fp, header);
	}
	if(type == MAT_SP || ((arr) && mxIsSparse(arr))){
	    if(sizeof(mwIndex)==4){
		magic=M_SP32;
	    }else{
		magic=M_SP64;
	    }
	    zfwrite(&magic, sizeof(uint32_t), 1, fp);
	    zfwrite(&m, sizeof(uint64_t), 1, fp);
	    zfwrite(&n, sizeof(uint64_t), 1, fp);
	    if(m!=0 && n!=0){
		mwIndex *Jc=mxGetJc(arr);
		long nzmax=Jc[n];
		zfwrite(&nzmax, sizeof(double), 1, fp);
		zfwrite(mxGetJc(arr), sizeof(mwIndex), n+1, fp);
		zfwrite(mxGetIr(arr), sizeof(mwIndex), nzmax, fp);
		zfwrite(mxGetPr(arr), sizeof(double), nzmax, fp);
	    }
	}else if(type == MAT_CSP || ((arr) && mxIsSparse(arr))){
	    if(sizeof(mwIndex)==4){
		magic=M_CSP32;
	    }else{
		magic=M_CSP64;
	    }
	    zfwrite(&magic, sizeof(uint32_t), 1, fp);
	    zfwrite(&m, sizeof(uint64_t), 1, fp);
	    zfwrite(&n, sizeof(uint64_t), 1, fp);
	    if(m!=0 && n!=0){
		mwIndex *Jc=mxGetJc(arr);
		long nzmax=Jc[n];
		zfwrite(&nzmax, sizeof(double), 1, fp);
		zfwrite(mxGetJc(arr), sizeof(mwIndex), n+1, fp);
		zfwrite(mxGetIr(arr), sizeof(mwIndex), nzmax, fp);
		zfwrite_complex(mxGetPr(arr),mxGetPi(arr),nzmax,fp);
	    }
	}else if(type == M_DBL || ((arr)&& mxIsDouble(arr))){
	    magic=M_DBL;
	    zfwrite(&magic, sizeof(uint32_t), 1, fp);
	    zfwrite(&m, sizeof(uint64_t), 1, fp);
	    zfwrite(&n, sizeof(uint64_t), 1, fp);
	    if(m!=0 && n!=0){
		zfwrite(mxGetPr(arr), sizeof(double), m*n, fp);
	    }  
	}else if(type == M_CMP || ((arr)&& mxIsDouble(arr))){
	    magic=M_CMP;
	    zfwrite(&magic, sizeof(uint32_t), 1, fp);
	    zfwrite(&m, sizeof(uint64_t), 1, fp);
	    zfwrite(&n, sizeof(uint64_t), 1, fp);
	    if(m!=0 && n!=0){
		zfwrite_complex(mxGetPr(arr),mxGetPi(arr), m*n, fp);
	    }
	}else{
	    error("Unrecognized data type");
	}
    }
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
	error("Usage: write(a,'a') or write(a, header, 'a')\n");
    }
    int nlen=mxGetM(prhs[ifn])*mxGetN(prhs[ifn])+1;
    char *fn=malloc(nlen);
    mxGetString(prhs[ifn],fn,nlen);
    fp=openfile(fn,"wb");
    free(fn);
    writedata(fp, 0, prhs[0], header);
    zfclose(fp);
}
