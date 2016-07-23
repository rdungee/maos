#pragma once
#include "numtype.h"
template <typename T>
class Dev{
    T p_;
public:
    void *operator new[](size_t size){
	void *p=::operator new[](size);
	memset(p, 0, size);//zeroing memory.
	return p;
    }
};
#if defined(DLONG)
///wrap x86 cmpxchg assembly code into a function. 8 Byte integer
inline Int cmpxchg(Int *ptr, Int old, Int newval){
    volatile Int *__ptr = (volatile Int *)(ptr);	
    Int __ret;                                     
    __asm__ volatile( "lock; cmpxchgq %2,%1"
		      : "=a" (__ret), "+m" (*__ptr) 
		      : "r" (newval), "0" (old)                     
		      : "memory");				 
    return __ret;
}
#else
///wrap x86 cmpxchg assembly code Into a function. 4 byte Integer
inline Int cmpxchg(Int *ptr, Int old, Int newval){
    volatile Int *__ptr = (volatile Int *)(ptr);	
    Int __ret;                                     
    __asm__ volatile( "lock; cmpxchg %2,%1"
		      : "=a" (__ret), "+m" (*__ptr) 
		      : "r" (newval), "0" (old)                     
		      : "memory");				 
    return __ret;
}
#endif
///Add val to ptr and return new value atomically using x86 assembly code.
inline Int AtomicAdd(Int *ptr, Int val){
    Int old;
    do{
	old=*ptr;
    }while(cmpxchg(ptr, old, old+val)!=old);
    return old+val;
}
///Magic number denoting the data type. Specialized for each type.
template <typename T>
struct Magic{
};
template <>
struct Magic<double>{
    enum { magic=M_DBL};
};
template <>
struct Magic<float>{
    enum { magic=M_FLT};
};
template <>
struct Magic<dcomplex>{
    enum { magic=M_CMP};
};
template <>
struct Magic<fcomplex>{
    enum { magic=M_ZMP};
};
template <>
struct Magic<int>{//always 32 bit int
    enum { magic=M_INT32};
};
template <>
struct Magic<long int>{//always 64 bit int. long is not.
    enum { magic=M_INT64};
};

//For backward compatibility without exposure class internal data.
template <typename T>
class pcompat{
private:
    T p_=0;
public:
    void SetP(T pi){
	p_=pi;
    }
    operator T() const{
	return p_;
    }
    T operator()() const {
	return p_;
    }
};

/**
   A base class for 2-d arrays. Not to be used directly.
 */
class TwoDim{
public:
private:
    Int nx_;       /**<Number of rows (usually x) (fast changing)*/
    Int ny_;       /**<Number of columns (usually y) (slow changing)*/
protected:
    explicit TwoDim(Int nxi=0, Int nyi=1):nx_(nxi),ny_(nxi?nyi:0),nx(nxi),ny(nxi?nyi:0){ };
public:
    std::string header; /**<header for optional data and description.*/
    virtual void MakeHeader(){}
    ///Return the number of rows.    
    Int Nx() const noexcept {return nx_;}
    ///Return the number of columns.
    Int Ny() const noexcept {return ny_;}
    ///Return the total number of elements
    Int N() const noexcept {return nx_*ny_;}
    ///Check whether the array is empty or not.

    ///Do not use operator bool(), which automatically convert array to Int.
    bool operator !()const noexcept{
	return nx_==0 || ny_==0;
    }
    void SetSize(Int nxi, Int nyi){
	nx=nx_=nxi;
	ny=ny_=nyi;
    }
    void SetSize(const TwoDim& in){
	nx=nx_=in.nx_;
	ny=ny_=in.ny_;
    }
    ///Temporary, for drop in compatibility.
public:
    Int nx;
    Int ny;
};

//Holding memory address and reference count.
template <typename T>
class RefP{
    T  *p_=0;
    Int *c_=0;
    void Unref_(){
	if(c_){
	    if(!AtomicAdd(c_, -1)){
		delete [] p_;
		delete c_;
	    }
	}
    }
    void Ref_(){
	if(c_){
	    AtomicAdd(c_, 1);
	}	
    }
public:
    RefP(T *pi=0):p_(pi){
	if(p_){
	    c_=new Int(1);
	}
    }
    ~RefP(){
	Unref_();
    }
    RefP(const RefP& in):p_(in.p_),c_(in.c_){
	Ref_();
    }
    RefP &operator=(const RefP& in){
	if(c_){
	    Unref_();
	}
	p_=in.p_;
	c_=in.c_;
	Ref_();
	return *this;
    }
    Int Nref()const{
	return c_?c_[0]:0;
    }
};
/**
   Array is responsible for memory allocation, reference counting, and automatic memory deallocation.

   It implements a smart pointer for automatic reference counting and free the
   memory when the last reference clears. By recording p0, we can deallocate
   the original pointer even when only a subarray is referenced in the new object.

   When nref_ is empty, we do not own the pointer and therefore won't free it.
*/
template <typename T>
class Array:public TwoDim{
private:
    T *p_=0;     /**<Records the original pointer from memory allocation.*/
    RefP<T> nref_; /**<Count the reference and records the original pointer*/
    Int ldx_=0;     /**<Leading dimension of memory allocation*/
    //Int offset_=0;  /**<Pointer offset*/
    mutable bool protect_=0;/**<Prevent overriding of p_*/
public:
    T* const &p;///Temporary, for drop in compatibility.
    struct fft_t *fft=0;
    struct mmap_t *mmap=0;
protected:
    ///Decrease counter and free memory. Keep other information intact.
    void CheckProtect(Int check_ref){
	if(protect_||(check_ref && Nref()>1)){
	    fatal("Cannot reassign/resize protected array\n");
	}
    }
    ///Use pin or allocate memory 
    void InitMem(T* pin){
	protect_=0;
	if(Nx()>0 && Ny()>0){
	    if(!pin){
		ldx_=Nx();
		p_=(T*)new Dev<T>[Nx()*Ny()];
	    }else{
		p_=pin;
		warn(p_ << "foreign");
	    }
	}
	nref_=RefP<T>(p_);
    }
    ///Reference memory, with optional offset.
    void RefMem(const Array &in, Int offset=0){
	//Prevent freeing memory during self assignment.
	RefP<T> tmp=in.nref_;
	p_=in.p_+offset;
	ldx_=in.ldx_;
	nref_=tmp;
	protect_=0;
	fft=0;
	mmap=0;
    }
public:

    ///Allocate memory and create reference counter [if size>0].
    explicit Array(Int nxi=0, int nyi=1, T*pin=0, int foreign=0)
	:TwoDim(nxi, nyi),p_(pin),ldx_(nx),p(p_){
	if(foreign){
	    if(!pin) fatal("When borrowing, p must not be NULL.\n");
	}else{
	    InitMem(pin);
	}
    }
   
    ///Copy constructor that does reference counting.
    Array(const Array&in) noexcept
	:TwoDim(in),p_(in.p_),nref_(in.nref_),ldx_(in.ldx_),p(p_){
    }
    ///Copy assignment operator.
    Array &operator=(const Array &in){
	CheckProtect(0);
	RefMem(in);
	SetSize(in);
    	return *this;	
    }
  
    ///Create a subarray starting at (ix0, iy0) with dimension (nx, ny). Memory is referenced.
    void Sub(const Array&in, Int ix0, Int nxi, Int iy0, Int nyi){
	//Reference the memory with offset
	RefMem(in, ix0+iy0*in.ldx_);
	//Adjust the size
	SetSize(nxi?nxi:(in.Nx()-ix0), nyi?nyi:(in.Ny()-iy0));
	//Set protection of input.
	in.protect_=1;
	assert((!nxi || !nyi) || (ix0>=0 && iy0>=0 && ix0<in.Nx() && iy0<in.Ny() && (ix0+Nx()<=in.Nx() ) && (iy0+Ny()<=in.Ny())));
    }
    ///Reference the input array. Also inherits the size. Use RefMem() to avoid copying sizing info.
    void Ref(const Array&in){
	//Reference the memory
	RefMem(in);
	//Copy the size.
	SetSize(in);
    }
    void Reshape(const Int nxi, const Int nyi){
	if(Nx()!=Ldx() && Ny()!=1){
	    fatal("Cannot reshape referenced sub array\n");
	}else if(Nx()*Ny()!=nxi*nyi){
	    fatal("Cannot array of resize"<<Nx()<<"x"<<Ny()<<" to "<<nxi<<"x"<<nyi);
	}else{
	    SetSize(nxi, nyi);
	    ldx_=nxi;
	}
    }
    ///Resize an array. New memory is allocated because T may be a class object
    ///type. Objects are copied by reference.
    Int Resize(const Int nxi, const Int nyi){
	CheckProtect(1);//This checks *nref_==1, which infers ldx_==nx_;
	//If total number of elements remain the same. Just change ldx_ and Nx();
	if(Nx()*Ny()!=nxi*nyi){
	    Array<T> pnew(nxi, nyi);
	    Int nymin=MIN(nyi, Ny());
	    Int nxmin=MIN(nxi, Nx());
	    for(Int iy=0; iy<nymin; iy++){
		for(Int ix=0; ix<nxmin; ix++){
		    pnew(ix, iy)=(*this)(ix, iy);
		}
	    }
	    *this=pnew;
	}else{
	    SetSize(nxi, nyi);
	    ldx_=nxi;
	}
	return 1;
    }
    
    ///Reinitialize with same size
    void New(){
	InitMem(0);
    }
    ///Reinit with different size
    void New(Int nxi, Int nyi){
	SetSize(nxi, nyi);
	InitMem(0);
    }
    void New(const TwoDim &size){New(size.Nx(), size.Ny());}
  
    ///Copy Content from one array to another.
    virtual ~Array(){
    }
    
public:
    //Implicit pointer retrival (for backward compatibility).
    operator const T*() const noexcept{ return p_; }
    operator T*() noexcept{ return p_; }
    //Explicit pointer retrieval
    T*P() noexcept {return (T*)p_;}
    const T*P() const noexcept {return (T*)p_;}
    //Explicit pointer retrieval
    T *operator()() noexcept {return P();}
    const T *operator()() const noexcept {return P();}
    //Retrieves original pointer (calloc output)
    const T*P0() const noexcept {return (T*)p_;}
    
    //Boundary checking is performed in debug mode. use -DNDEBUG to disable
    //boundary checking. causes 3 times slow down

    ///Retrieves pointer to elements.
    T* P(Int ix, Int iy){
	assert(ix>=0 && ix<Nx() && iy >=0 && iy<Ny());
	return P()+ix+iy*ldx_;
    }
    const T* P(Int ix, Int iy)const{
	assert(ix>=0 && ix<Nx() && iy >=0 && iy<Ny());
	return P()+ix+iy*ldx_;
    }
    
    ///Use to retrieve elements in pointer like fashion. For const arrays. Operation is invalid if ldx>nx. 
    const T &operator [](Int i)const{
	assert(ldx_==Nx() && i >=0 && i < Nx()*Ny());//Sanity check in debug mode.
	return P()[i];
    }
    
    ///Use to retrieve elements in 1d function call fashion. For const arrays.
    const T &operator ()(Int i)const noexcept{
	return operator[](i);
    }
    
    ///Use to retrieve elemnets in 2d function call fashion. For const arrays. 
    const T &operator ()(Int ix, Int iy)const{
	//Performs boundary check in debug mode.
	assert(ix>=0 && ix<Nx() && iy >=0 && iy<Ny());
	return P()[ix+iy*ldx_];
    }
  
    //Non const versions uses const version equivalent.
    ///Use to retrieve elements in pointer like fashion
    T &operator [](Int ix){
	return const_cast<T&>(static_cast<const Array&>(*this)[ix]);
    }
    ///Use to retrieve elements in 1d function call fashion
    T &operator ()(Int ix){
	return const_cast<T&>(static_cast<const Array&>(*this)(ix));
    }
    ///Use to retrieve elemnets in 2d function call fashion
    T &operator ()(Int ix, Int iy){
	return const_cast<T&>(static_cast<const Array&>(*this)(ix, iy));
    }
    ///Insert n elements into a vector array. Resize it simultaneously if resize is true
    void Insert(Int ix, Int nelem, bool resize=true){
	if(!nelem) return;
	if(resize){
	    if(Ny()<=1){
		Resize(Nx()+nelem, 1);
	    }else if(Nx()<=1){
		Resize(1, Ny()+nelem);
	    }else{
		error("Can only insert into a vector\n");
	    }
	}
	for(Int i=N()-nelem-1; i>=ix; i--){
	    p_[i+nelem]=p_[i];
	}
	for(Int i=ix; i<ix+nelem; i++){
	    p_[i]=0;
	}
    }
    ///Returns the value of reference counter
    Int Nref()const noexcept{ return nref_.Nref(); }
    ///Return the leading dimension of underlying memory.
    Int Ldx() const noexcept {return ldx_;}

    //Returns column pointer.
    T *Col(Int icol) {assert(icol>=0 && icol<Ny());return P()+icol*ldx_;}
    const T *Col(Int icol) const {assert(icol>=0 && icol<Ny()); return P()+icol*ldx_;}
    
    ///Use to check whether two arrays points to the same memory location and is the same dimension.
    bool operator==(const Array&B)const noexcept{
	return (Nx()==B.Nx() && Ny()==B.Ny() && p_==B.p_);
    }
    ///Inverse of ==.
    bool operator!=(const Array&B)const noexcept{
	return !operator==(B);
    }
};
///Display Data
template <typename T>
ostream& operator<<(ostream &os, const Array<T> &in){
    os<<in.Nx()<<"x"<<in.Ny()<<": ("<<typeid(T).name()<<")"<<endl;
    for(Int ix=0; ix<MIN(10,in.Nx()); ix++){
	for(Int iy=0; iy<MIN(10,in.Ny()); iy++){
	    os<<in(ix,iy)<<" ";
	}
	if(in.Ny()>10){
	    os<<"...";
	}
	if(in.Ny()){
	    os<<endl;
	}
    }
    if(in.Nx()>10){
	os<<"..."<<endl;
    }
    if(in.Nx())
	os<<endl;
    return os;
}

/**
   Dense matrix of fundmental numerical data types.
*/

template <typename T>
class Mat:public Array<T>{
public:
    using Array<T>::Array;//Inherits constructor (C++11)
    using Array<T>::CheckProtect;
    using Array<T>::P;
    using Array<T>::Col;
    using Array<T>::Ldx;
    using Array<T>::Nx;
    using Array<T>::Ny;
    using Array<T>::N;
public:
  
};

typedef Mat<Real> RMat;
typedef Mat<Comp> CMat;
typedef Mat<Int>  IMat;

/**
   Dense matrix of objects.
*/
template <typename T>
class Cell:public Array<T>{
private:
    T m_;  /**<Stores continuous data of the underlying matrix for block matrix. May be empty*/
public:
    using Array<T>::P;
    using Array<T>::Ldx;
    using Array<T>::Nx;
    using Array<T>::Ny;
    using Array<T>::N;

public:
    ///Allocate an array
    Cell(Int nxi=0, Int nyi=1):Array<T>(nxi, nyi){}
    virtual ~Cell(){
    }
       
    ///Return the number of rows of each element.
    IMat Nxs() const{
	IMat size=IMat(Nx(), Ny());
	for(Int iy=0; iy<Ny(); iy++){
	    for(Int ix=0; ix<Nx(); ix++){
		size(ix,iy)=(*this)(ix, iy).Nx();
	    }
	}
	return size;
    }

    ///Return the number of cols of each element.
    IMat Nys() const{
	IMat size=IMat(Nx(), Ny());
	for(Int iy=0; iy<Ny(); iy++){
	    for(Int ix=0; ix<Nx(); ix++){
		size(ix,iy)=(*this)(ix, iy).Ny();
	    }
	}
	return size;
    }
    
  
    ///Create a cell as a block matrix. Each sub-cell references the big
    ///array. Works with both vector and matrix All elements are of the same
    ///size.
    void NewDeep(Int nxi, Int nyi, Int nnx, Int nny){
	this->New(nxi,nyi);
	m_.New(nxi*nnx, nyi*nny);
	for(Int iy=0; iy<nyi; iy++){
	    for(Int ix=0; ix<nxi; ix++){
		(*this)(ix, iy).Sub(m_, ix*nnx, nnx, iy*nny, nny);
	    }
	}
    }
    ///Create a cell as a block matrix. Each sub-cell references the big
    ///array. Works with both vector and matrix All sub-cells can be different
    ///size, but must match in shape like block matrix.
    void NewDeep(const Array<Int> &nnx, const Array<Int> &nny){
	if(nnx.Nx()!=nny.Nx() || nny.Ny()!=nny.Ny()){
	    fatal("nnx and nny mismatch");
	}
	const Int nxi=nnx.Nx();
	const Int nyi=nny.Ny();
	this->New(nxi, nyi);
	Int totx=0, toty=0;
	IMat nnx2(nxi,1), nny2(nyi,1);
	for(Int iy=0; iy<nyi; iy++){
	    for(Int ix=0; ix<nxi; ix++){
		if(iy==0){//first column. compute row size
		    totx+=nnx(ix,0);
		    nnx2[ix]=nnx(ix,0);
		}
		if(ix==0){//first row. compute column size.
		    toty+=nny(0, iy);
		    nny2[iy]=nny(0, iy);
		}
		if((nnx(ix, iy)!=0 && nny(ix,iy)!=nny2[iy])
		   || (nny(ix,iy)!=0 && nnx(ix, iy)!=nnx2[ix])){ 
		    debug(nnx);
		    debug(nny);
		    fatal("cell ("<<ix<<","<<iy<<") should be ("<<nnx(ix,0)<<","<<nny(0,iy)<<")"
			  <<" but is ("<<nnx(ix, iy)<<","<<nny(ix, iy)<<")"<<endl);
		}
	    }
	}
	m_.New(totx, toty);
	toty=0;
	for(Int iy=0; iy<nyi; iy++){
	    totx=0;
	    for(Int ix=0; ix<nxi; ix++){
		(*this)(ix, iy).Sub(m_, totx, nnx2[ix], toty, nny2[iy]);
		totx+=nnx2[ix];
	    }
	    toty+=nny2[iy];
	}
    }
    
    ///Create another array with the same size
    void NewDeep(const Cell &in){
	NewDeep(in.Nxs(), in.Nys());
    }

    //Use base class destructor, copy/move constructor and copy/move assignment operators.
    
    ///Return the continuous array.
    const T& M()const{
	return m_;
    }
    
    ///Return the continuous array as a new reference. 
    T& M(){
	return m_;
    }
};
//Compresed Column Storage Sparse matrix. Always properly ordered.
template <typename T>
class Sparse:public TwoDim{
private:
    Mat<Int> p_; ///<Column pointers (size ny+1)
    Mat<Int> i_; ///<Row indices, size nzmax
    Mat<T> x_;   ///<Numerical values, size nzmax
protected:
    void SetCompat();//Set compatibility flags.
public:
    explicit Sparse(Int nxi=0, Int nyi=0, Int nzm=0)
	:TwoDim(nxi, nyi), p_(nyi+1, 1), i_(nzm, 1),x_(nzm, 1){
	SetCompat();
    }
    explicit Sparse(Int nxi, Int nyi, Int nzm, Int *pi, Int *ii, T *xi)
	:TwoDim(nxi,nyi), p_(nyi+1, 1, pi), i_(nzm, 1, ii), x_(nzm, 1, xi){
	SetCompat();
    }
    Sparse(const Sparse &in)noexcept:TwoDim(in),p_(in.p_),i_(in.i_),x_(in.x_){
	SetCompat();
    }

    Sparse &operator=(const Sparse &in){
	if(p_.P0() && p_.P0()!=in.p_.P0()){
	    SetSize(in);
	    p_=in.p_;
	    i_=in.i_;
	    x_=in.x_;
	}
	SetCompat();
    }
    ~Sparse(){};
    Int Resize(Int nzm){
	int ans=i_.Resize(nzm,1) && x_.Resize(nzm,1);
	SetCompat();
	return ans;
    }
    Int Resize(Int nxi, Int nyi, Int nzm){
	SetSize(nxi, nyi);
	return p_.Resize(nyi+1,1) && Resize(nzm);
    }
    Mat<Int>&P() {return p_;}
    const Mat<Int>&P()const {return p_;}
    Int P(Int ii) const {return p_(ii);}

    Mat<Int>&I() {return i_;}
    const Mat<Int>&I()const {return i_;}
    Int I(Int ii) const {return i_(ii);}
    
    Mat<T>& X() {return x_;}
    const Mat<T>& X()const {return x_;}
    T X(Int ii) const {return x_(ii);}
    
    Int Nzmax() {return i_.N();}
    //Slow. Only use for debugging.
    T operator()(Int ix, Int iy) const{
	assert(ix>=0 && ix<Nx() && iy>=0 && iy<Ny());
	for(Int ic=P(iy); ic<P(iy+1); ic++){
	    if(I(ic)==ix){
		return X(ic);
	    }
	}
	return 0;
    }
    //Slow. Only use for debugging.
    void Set(Int ix, Int iy, T val) {
	assert(ix>=0 && ix<Nx() && iy>=0 && iy<Ny());
	Int found=0;
	Int ic;
	for(ic=p_(iy); ic<p_(iy+1); ic++){
	    if(i_(ic)==ix){//Already exist.
		x_(ic)=val;
		found=1;
		break;
	    }
	}
	if(!found){
	    i_.Insert(ic, 1, i_.Ny()<=p_(Ny()));
	    x_.Insert(ic, 1, x_.Ny()<=p_(Ny()));
	    for(Int ii=iy+1; ii<Ny()+1; ii++){
		p_(ii)++;
	    }
	    i_(ic)=ix;
	    x_(ic)=val;
	}
    }

public: //For backward compatibility
    pcompat<Int*> p;
    pcompat<Int*> i;
    pcompat<T*> x;
    pcompat<Int> nzmax;
    
};

///Display Data
template <typename T>
ostream& operator<<(ostream &os, const Sparse<T> &in){
    os<<in.Nx()<<"x"<<in.Ny()<<":"<<endl;
    for(Int iy=0; iy<in.Ny(); iy++){
	for(Int ic=in.P(iy); ic<in.P(iy+1); ic++){
	    os<<"("<<in.I(ic)<<","<<iy<<"): "<<in.X(ic)<<endl;
	}
    }
    return os;
}

///trait class to itentify the underline numerical type. 
template <typename T>
struct num_type{
    typedef T type;
};

///partial specialization for trait class
template <typename T>
struct num_type<Array<T>>{
    typedef typename num_type<T>::type type;
};


///partial specialization for trait class 
template <typename T>
struct num_type<Mat<T>>{
    typedef typename num_type<T>::type type;
};


///partial specialization for trait class
template <typename T>
struct num_type<Cell<T>>{
    typedef typename num_type<T>::type type;
};
template <>
struct Magic<Sparse<double>>{
    enum {magic=M_DSP};
};
template <>
struct Magic<Sparse<float>>{
    enum {magic=M_SSP};
};
template <>
struct Magic<Sparse<dcomplex>>{
    enum {magic=M_CSP};
};
template <>
struct Magic<Sparse<fcomplex>>{
    enum {magic=M_ZSP};
};

template <typename T>
struct Magic<Array<T>>{
    enum { magic=MCC_ANY};
};

template <typename T>
struct Magic<Mat<T>>{
    enum { magic=MCC_ANY};
};
template <typename T>
struct Magic<Cell<T>>{
    enum { magic=MCC_ANY};
};

template <typename T>
struct Magic<Array<T>*>{
    enum { magic=MCC_ANY};
};

template <typename T>
struct Magic<Mat<T>*>{
    enum { magic=MCC_ANY};
};
template <typename T>
struct Magic<Cell<T>*>{
    enum { magic=MCC_ANY};
};
typedef Cell<RMat> RCell;
typedef Cell<CMat> CCell;
typedef Cell<IMat> ICell;

template <typename T>
void Sparse<T>::SetCompat(){
    p.SetP(p_());
    i.SetP(i_());
    x.SetP(x_());
    nzmax.SetP(Nzmax());
}

template <typename T>
void WriteHeader(const Array<T>&A, File &fp){
    Header header;
    header.magic=Magic<T>::magic;
    header.nx=A.Nx();
    header.ny=A.Ny();
    header.str=A.header;
    fp.WriteHeader(header);
}

template <typename T>
void WriteData(const Mat<T>& A, File &fp){
    WriteHeader(A, fp);
    if(A.Ldx()==A.Nx()){
	fp.Write(A.P(), sizeof(T), A.Nx()*A.Ny());
    }else{
	for(Int iy=0; iy<A.Ny(); iy++){
	    fp.Write(A.Col(iy), sizeof(T), A.Nx());
	}
    }
}

template <typename T>
void ReadData(Mat<T>& A, File &fp){
    Header header;
    fp.ReadHeader(header);
    if(A.Nx()!=(Int)header.nx || A.Ny()!=(Int)header.ny){
	A=Mat<T>(header.nx, header.ny);
    }
    A.header=header.str;
    if(A.Ldx()==A.Nx()){
	fp.Read(A.P(), sizeof(T), A.Nx()*A.Ny());
    }else{
	for(Int iy=0; iy<A.Ny(); iy++){
	    fp.Read(A.Col(iy), sizeof(T), A.Nx());
	}
    }
}

template <typename T>
void WriteData(const Cell<T>& A, File &fp){
    WriteHeader(A, fp);
    for(Int iy=0; iy<A.Ny(); iy++){
	for(Int ix=0; ix<A.Nx(); ix++){
	    WriteData(A(ix, iy), fp);
	}
    }
}
template <typename T>
void ReadData(Cell<T>& A, File &fp){
    Header header;
    fp.ReadHeader(header);
    if(A.Nx()!=(Int)header.nx || A.Ny()!=(Int)header.ny){
	A=Cell<T>(header.nx, header.ny);
    }
    for(Int iy=0; iy<A.Ny(); iy++){
	for(Int ix=0; ix<A.Nx(); ix++){
	    ReadData(A(ix, iy), fp);
	}
    }
}

///Save to file. 
template <typename T>
void WriteBin(const T&A, const char *format, ...){
    format2fn;
    File fp(fn, "wb");
    WriteData(A, fp);
}
//Read from file.
template <typename T>
void ReadBin(T&A, const char *format, ...){
    format2fn;
    File fp(fn, "rb");
    ReadData(A, fp);
}


///Check whether two arrys match in size.
inline void CheckMatch(const TwoDim &A, const TwoDim &B){
    if(A.Nx()!=B.Nx() || A.Ny() != B.Ny()){
	fatal("Mismatch: "<<A.Nx()<<"x"<<A.Ny()<<" != "<<B.Nx()<<"x"<<B.Ny());
    }
}


///Copy from one array to another
template <typename T>
void Copy(Mat<T> &out, const Mat<T> &in){
    if(out.P0()==in.P0()){
	fatal("Cannot copy overlapped region\n");
    }
    if(!out) out.New(in);
    else CheckMatch(out, in);
    if(out.Ldx() == out.Nx() && in.Ldx()==in.Nx()){
	memcpy(out.P(), in.P(), sizeof(T)*in.Nx()*in.Ny());
    }else{
	for(Int icol=0; icol<in.Ny(); icol++){
	    memcpy(out.Col(icol), in.Col(icol), sizeof(T)*in.Nx());
	}
    }	
}
///Copy from one array to another
template <typename T>
void CopyOverlap(Mat<T> &out, const Mat<T> &in){
    if(!out){
	Copy(out, in);
    }else{
	if(out.P0()==in.P0()){
	    fatal("Cannot copy overlapped region\n");
	}
	Int nymin=MIN(out.Ny(), in.Ny());
	Int nxmin=MIN(out.Ny(), in.Nx());
	//Column size match
	if(out.Ldx() == out.Nx() && in.Ldx()==in.Nx() && in.Nx()==out.Nx()){
	    memcpy(out.P(), in.P(), sizeof(T)*in.Nx()*nymin);
	}else{
	    for(Int icol=0; icol<nymin; icol++){
		memcpy(out.Col(icol), in.Col(icol), sizeof(T)*nxmin);
	    }
	}
    }	
}
///Deep Copy from one cell to another
template <typename T>
void CopyDeep(Cell<T> &out, const Cell<T> &in){
    if(out.P0()==in.P0()){
	fatal("Cannot copy overlapped region\n");
    }
    if(!out) out.NewDeep(in);
    else CheckMatch(out, in);
    for(Int iy=0; iy<in.Ny(); iy++){
	for(Int ix=0; ix<in.Nx(); ix++){
	    Copy(out(ix, iy), in(ix, iy));
	}
    }
}

///Shallow (reference) Copy from one cell to another
template <typename T>
void Copy(Cell<T> &out, const Cell<T> &in){
    if(out.P0()==in.P0()){
	fatal("Cannot copy overlapped region\n");
    }
    if(!out) out.New(in);
    else CheckMatch(out, in);
    for(Int iy=0; iy<in.Ny(); iy++){
	for(Int ix=0; ix<in.Nx(); ix++){
	    out(ix, iy)=in(ix, iy);
	}
    }
}

///Shallow (reference) Copy from one cell to another in overlapping region only.
template <typename T>
void CopyOverlap(Cell<T> &out, const Cell<T> &in){
    if(!out){
	Copy(out, in);
    }else{
	if(out.P0()==in.P0()){
	    fatal("Cannot copy overlapped region\n");
	}
	Int nymin=MIN(out.Ny(), in.Ny());
	Int nxmin=MIN(out.Ny(), in.Nx());
	for(Int iy=0; iy<nymin; iy++){
	    for(Int ix=0; ix<nxmin; ix++){
		out(ix, iy)=in(ix, iy);
	    }
	}
    }
}

///Modifying element by element
template <typename T, typename F>
void Apply(Mat<T>& A, F fun){
    for(Int iy=0; iy<A.Ny(); iy++){
	for(Int ix=0; ix<A.Nx(); ix++){
	    fun(A(ix, iy));
	}
    }
}

///Reduction operation for all elements without modifying
template <typename T, typename F>
void Apply(const Mat<T>& A, F fun){
    for(Int iy=0; iy<A.Ny(); iy++){
	for(Int ix=0; ix<A.Nx(); ix++){
	    fun(A(ix, iy));
	}
    }
}
///Modifying only diagonal elements
template <typename T, typename F>
void ApplyDiag(Mat<T>& A, F fun){
    if(A.Ny()!=A.Nx()) fatal("Not square:"<<A.Nx()<<"x"<<A.Ny()<<endl);
    for(Int ix=0; ix<A.Nx(); ix++){
	fun(A(ix, ix));
    }
}
///Reduction operation for only diagonal elements
template <typename T, typename F>
void ApplyDiag(const Mat<T>& A, F fun){
    if(A.Ny()!=A.Nx()) fatal("Not square:"<<A.Nx()<<"x"<<A.Ny()<<endl);
    for(Int ix=0; ix<A.Nx(); ix++){
	fun(A(ix, ix));
    }
}
//Two arrays of same size
///Modifying first array element by element
template <typename T, typename F>
void Apply2(Mat<T>& A, const Mat<T>& B, F fun){
    if(!B) return; //no-op
    else if(!A){
	A=Mat<T>(B.Nx(), B.Ny());
    }else{
	CheckMatch(A,B);
    }
    for(Int iy=0; iy<B.Ny(); iy++){
	for(Int ix=0; ix<B.Nx(); ix++){
	    fun(A(ix, iy), B(ix, iy));
	}
    }
}
///Reduction operation for all elements
template <typename T, typename F>
void Apply2(const Mat<T>& A, const Mat<T>& B, F fun){
    CheckMatch(A,B);
    for(Int iy=0; iy<B.Ny(); iy++){
	for(Int ix=0; ix<B.Nx(); ix++){
	    fun(A(ix, iy), B(ix, iy));
	}
    }
}

//Three arrays of same size
///Modifying element by element
template <typename T, typename F>
void Apply3(Mat<T>& A, const Mat<T>& B, const Mat<T>& C, F fun){
    if(!B && !C) return; //no-op
    else if(!A){
	A=Mat<T>(B.Nx(), B.Ny());
    }else{
	CheckMatch(A,B);
    }
    CheckMatch(B,C);
    for(Int iy=0; iy<B.Ny(); iy++){
	for(Int ix=0; ix<B.Nx(); ix++){
	    fun(A(ix, iy), B(ix, iy),C(ix,iy));
	}
    }
}

///Reduction operation for all elements
template <typename T, typename F>
void Apply3(const Mat<T>& A, const Mat<T>& B, const Mat<T>& C, F fun){
    CheckMatch(A,B);
    CheckMatch(B,C);
    for(Int iy=0; iy<B.Ny(); iy++){
	for(Int ix=0; ix<B.Nx(); ix++){
	    fun(A(ix, iy), B(ix, iy), C(ix, iy));
	}
    }
}

///Modifying element by element
template <typename T, typename F>
void Apply(Cell<T>& A, F fun){
    for(Int iy=0; iy<A.Ny(); iy++){
	for(Int ix=0; ix<A.Nx(); ix++){
	    Apply(A(ix, iy), fun);
	}
    }
}
///Reduction operation for all elements
template <typename T, typename F>
void Apply(const Cell<T>& A, F fun){
    for(Int iy=0; iy<A.Ny(); iy++){
	for(Int ix=0; ix<A.Nx(); ix++){
	    Apply(A(ix, iy), fun);
	}
    }
}
///Modifying only diagonal elements
template <typename T, typename F>
void ApplyDiag(Cell<T>& A, F fun){
    if(A.Ny()!=A.Nx()) fatal("Not square");
    for(Int ix=0; ix<A.Nx(); ix++){
	Apply(A(ix, ix), fun);
    }
}
///Reduction operation for only diagonal elements
template <typename T, typename F>
void ApplyDiag(const Cell<T>& A, F fun){
    if(A.Ny()!=A.Nx()) fatal("Not square");
    for(Int ix=0; ix<A.Nx(); ix++){
	Apply(A(ix, ix), fun);
    }
}
//Two arrays of same size. Use different name for easy error diagnosis during compiling.
///Modifying element by element
template <typename T, typename F>
void Apply2(Cell<T>& A, const Cell<T>& B, F fun){
    if(!A){
	A=Cell<T>(B.Nx(), B.Ny());
    }else{
	CheckMatch(A,B);
    }
    for(Int iy=0; iy<A.Ny(); iy++){
	for(Int ix=0; ix<A.Nx(); ix++){
	    Apply2(A(ix, iy), B(ix, iy), fun);
	}
    }
}

///Reduction operation for all elements
template <typename T, typename F>
void Apply2(const Cell<T>& A, const Cell<T>& B, F fun){
    CheckMatch(A,B);
    for(Int iy=0; iy<A.Ny(); iy++){
	for(Int ix=0; ix<A.Nx(); ix++){
	    Apply2(A(ix, iy), B(ix, iy),  fun);
	}
    }
}     
  
//Three arrays of same size
///Modifying element by element
template <typename T, typename F>
void Apply3(Cell<T>& A, const Cell<T>& B, const Cell<T>& C, F fun){
    if(!A){
	A=Cell<T>(B.Nx(), B.Ny());
    }else{
	CheckMatch(A,B);
    }
    CheckMatch(B,C);
    for(Int iy=0; iy<B.Ny(); iy++){
	for(Int ix=0; ix<B.Nx(); ix++){
	    Apply3(A(ix, iy), B(ix, iy),C(ix,iy), fun);
	}
    }
}

///Reduction operation for all elements
template <typename T, typename F>
void Apply3(const Cell<T>& A, const Cell<T>& B, const Cell<T>& C, F fun){
    CheckMatch(A,B);
    CheckMatch(B,C);
    for(Int iy=0; iy<B.Ny(); iy++){
	for(Int ix=0; ix<B.Nx(); ix++){
	    Apply3(A(ix, iy), B(ix, iy), C(ix, iy), fun);
	}
    }
}

///Set each element to alpha.
template <typename T>
void Set(T& A, typename num_type<T>::type alpha){
    Apply(A, [=](typename num_type<T>::type &val){val=alpha;});
}
///Scale each element by alpha.
template <typename T>
void Scale(T& A, typename num_type<T>::type alpha){
    typedef typename num_type<T>::type S;
    if(alpha==(S)0){
	Set(A, 0);
    }else{
	Apply(A, [=](S &val){val*=alpha;});
    }
}
///Sum all the elements
template <typename T>
auto Sum(const T& A)->typename num_type<T>::type{
    typedef typename num_type<T>::type S;
    S res=0;
    Apply(A, [&](S val){res+=val;});
    return res;
}
///Sum after applying a function to each elements.
template <typename T, typename F>
auto Sum(const T& A, F fun)->typename num_type<T>::type{
    typedef typename num_type<T>::type S;
    S res=0;
    Apply(A, [&](S val){res+=fun(val);});
    return res;
}
///Sum all the elements after taking absolute value
template <typename T>
auto SumAbs(const T& A)->typename num_type<T>::type{
    //return Sum(A, abs);
    typedef typename num_type<T>::type S;
    return Sum(A, [](S val){return fabs(val);});
}				 
///Sum all the diagonal elements
template <typename T>
auto SumDiag(const T& A)->typename num_type<T>::type{
    typedef typename num_type<T>::type S;
    S res=0;
    ApplyDiag(A, [&](S val)->void{res+=val;});
    return res;
}

///Sum after multiplying elements by elements.
template <typename T>
auto Dot(const T& A, const T& B)->typename num_type<T>::type{
    typedef typename num_type<T>::type S;
    S res=0;
    Apply2(A, B, [&](S val, S val2){res+=val*val2;});
    return res;
}

///Add a scalar to every diagonal element.
template <typename T>
void AddDiag(T& A, typename num_type<T>::type alpha){
    typedef typename num_type<T>::type S;
    ApplyDiag(A, [=](S& val1){val1+=alpha;});
}

///Add a scalar to every element.
template <typename T>
void Add(T& A, typename num_type<T>::type alpha){
    typedef typename num_type<T>::type S;
    Apply(A, [=](S& val1){val1+=alpha;});
}

///Component wise A=A*alpha+B*beta.
template <typename T>
void Add(T& A, typename num_type<T>::type alpha, const T& B, typename num_type<T>::type beta){
    typedef typename num_type<T>::type S;
    Apply2(A, B, [=](S& val1, S val2){val1=val1*alpha+val2*beta;});
}

///Component wise A=A*B.
template <typename T>
void Cwm(T& A, const T& B){
    typedef typename num_type<T>::type S;
    Apply2(A, B, [=](S& val1, S val2){val1=val1*val2;});
}

///Component wise A=A*B*beta.
template <typename T>
void Cwm(T& A, const T& B, typename num_type<T>::type beta){
    typedef typename num_type<T>::type S;
    Apply2(A, B, [=](S& val1, S val2){val1=val1*val2*beta;});
}
