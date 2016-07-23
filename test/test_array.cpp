#include <tgmath.h>
#include "catch.h"
#include "../math/array.h"
using std::abs;
#define REQUIRE_EQ(A,B) REQUIRE(abs(A-B)<=abs(A+B)*1e-14)
TEST_CASE("array", "[array]"){
    const Int N=6;
    //Test constructors
    typedef Array<double> DArr;
    SECTION("empty object"){
	DArr a0;
	REQUIRE(a0.P()==0);
	REQUIRE(!a0);
	REQUIRE(a0.N()==0);
	REQUIRE(a0.Nref()==0);
	a0=a0;//Self assignment
    }
    //Create a new object
    DArr aa(N,N);
    SECTION("Basic"){
	REQUIRE(aa.Nref()==1);
	REQUIRE(aa.Nx()==N);
	REQUIRE(aa.Ny()==N);
	REQUIRE(aa.Ldx()==N);
	REQUIRE(aa.N()==N*N);
	REQUIRE(aa.P()==aa.p);
    }
    SECTION("weak reference"){
	//Weak reference from aa.
	DArr ab(N,N,aa.P(),1);
	REQUIRE(ab.P()==aa.P());
	REQUIRE(aa.Nref()==1);
	REQUIRE(ab.Nref()==0);
	ab=ab;//Self assignment with weak reference
	REQUIRE(ab.Nref()==0);
	DArr ac=ab;//Copy constructor
	REQUIRE(ac.Nref()==0);
	ab=aa; //Assignment operator
	REQUIRE(ab.Nref()==aa.Nref());
	REQUIRE(ab.Nref()==2);
    }
    SECTION("strong reference"){
	aa=aa;//Self assignment with strong reference
	REQUIRE(aa.Nref()==1);
	//Copy constructor
	DArr ac(aa);
	REQUIRE(aa.P()==ac.P());
	REQUIRE(aa.Nref()==2);
	ac.Ref(ac);
	ac.Ref(aa);//Self reference. Nothing changes
	REQUIRE(aa.P()==ac.P());
	REQUIRE(aa.Nref()==2);
    }
    SECTION("Sub array"){
	DArr ab;
	ab.Sub(aa, 2, 4, 3, 3);
	REQUIRE(ab.P()==aa.P(2,3));
	REQUIRE(ab.Nx()==4);
	REQUIRE(ab.Ny()==3);
	REQUIRE(ab.Ldx()==aa.Ldx());
	REQUIRE(ab.Nref()==2);
	REQUIRE_THROWS(aa=ab);
	ab=aa;
	REQUIRE(ab.P()==aa.P());
	REQUIRE(aa.Nref()==2);
	ab=DArr();
	REQUIRE(aa.Nref()==1);
    }
    SECTION("Reshape"){
	aa.Reshape(N*2, N/2);
	REQUIRE(aa.N()==N*N);
	REQUIRE(aa.Nx()==N*2);
	REQUIRE(aa.Nref()==1);
    }
     SECTION("NEW"){
	REQUIRE(aa.Nref()==1);
	aa.New();
	REQUIRE(aa.N()==N*N);
	aa.New(N*2,N*2);
	REQUIRE(aa.N()==N*N*4);
	REQUIRE(aa(0,0)==0);
     }
     
     SECTION("Resize"){
	aa(0,0)=1;
	aa(N-1,1)=2;
	aa.Resize(N*2,N);
	REQUIRE(aa(0,0)==1);
	REQUIRE(aa(1,0)==0);
	REQUIRE(aa(N-1,1)==2);
	REQUIRE(aa(N,1)==0);
     }
}
TEST_CASE("mat", "[mat]"){
    const Int N=4;
    RMat aa(N,N);
    Set(aa, -2); //Set all elements to 2.
    SECTION("Display"){
	debug(aa);
	REQUIRE_EQ(aa(1,2),-2);
    }
    SECTION("Resize"){
	aa.Resize(N*2, N*2);
	REQUIRE(aa(N-1,N-1)==-2);
	REQUIRE(aa(N,N)==0);
    }
    SECTION("Sum"){
	REQUIRE(Sum(aa)==-2*N*N);
	REQUIRE(SumAbs(aa)==2*N*N);
	REQUIRE(SumDiag(aa)==-2*N);
    }
    SECTION("Dot"){
	REQUIRE(Dot(aa,aa)==4*N*N);
    }
    SECTION("AddDiag"){
	AddDiag(aa, 4);
	REQUIRE(SumDiag(aa)==2*N);
	REQUIRE(Sum(aa)==4*N-2*N*N);
    }
    SECTION("ADD"){
	Add(aa,4);
	REQUIRE(Sum(aa)==2*N*N);
	RMat bb(N,N);
	Set(bb, 1);
	Add(aa, 0.5, bb, 2);
	REQUIRE(Sum(aa)==N*N*3);
    }
    SECTION("Cwm"){
	RMat bb(N,N);
	Set(bb, 3);
	Cwm(aa, bb);
	REQUIRE(Sum(aa)==-6*N*N);
	Cwm(aa, bb, 2);
	REQUIRE(Sum(aa)==-36*N*N);
    }
    SECTION("Insert"){
	RMat ab(6,1);
	Set(ab,2);
	ab.Insert(2, 2);//Insert and expand
	REQUIRE(ab.N()==8);
	REQUIRE(Sum(ab)==12);
	REQUIRE(ab(3)==0);
	REQUIRE(ab(4)==2);
	ab.Insert(2, 2, 0);//Insert without expand.
	REQUIRE(Sum(ab)==8);
	REQUIRE(ab.N()==8);
    }
}

TEST_CASE("cell","[cell]"){
    const int N=6;
    SECTION("Basic"){
	RCell ac;
	RCell bc(N,N*2);
	REQUIRE(bc.N()==N*N*2);
	REQUIRE(bc.Nx()==N);
	REQUIRE(bc.Ny()==N*2);
	REQUIRE(!bc(1,2));

	RMat a;
	a.New(bc);
	REQUIRE(a.N()==2*N*N);
    }
    SECTION("NewDeep 1"){
	RCell ac;
	ac.NewDeep(N,N*2,3,2);
	REQUIRE(ac.M().Nx()==N*3);
	REQUIRE(ac.M().Ny()==N*4);
	RCell bc=ac;
	REQUIRE(ac.M().P()==bc.M().P());
	REQUIRE_THROWS(ac.M()=RMat(N,N));
	RMat aa=ac.M();
	REQUIRE(aa.Nx()==N*3);
    }
    SECTION("New Deep 2"){
	RCell ac;
	Array<Int>nnx(N, N);
	Array<Int>nny(N, N);
	ac.NewDeep(nnx, nny);
	assert(ac.M().N()==0);
	nnx(0,0)=1; nny(0,0)=1;
	nnx(1,0)=2; nny(1,0)=1;
	nnx(0,1)=1; nny(0,1)=3;
	nnx(1,1)=2; nny(1,1)=3;
	ac.NewDeep(nnx, nny);
	Set(ac, 2);
	REQUIRE(ac.M().Nx()==3);
	REQUIRE(ac.M().Ny()==4);
	REQUIRE(&ac.M()(1,1)==&ac(1,1)(0,0));
	REQUIRE(!ac(2,2));
	RCell bc;
	bc.NewDeep(ac);
	REQUIRE(bc.P()!=ac.P());
	REQUIRE(bc.M().N()==ac.M().N());
	REQUIRE(&bc.M()(1,1)==&bc(1,1)(0,0));

	bc.NewDeep(ac.Nxs(), ac.Nys());
	REQUIRE(&bc.M()(1,1)==&bc(1,1)(0,0));
    }
    SECTION("Math"){
	RCell ac;
	ac.NewDeep(N,N*2,3,2);	
	Set(ac, 1);
	REQUIRE(ac(2,2)(1,1)==1);
	Scale(ac, 2);
	REQUIRE(ac(2,2)(1,1)==2);
	REQUIRE(Sum(ac)==2*ac.N()*6);
	RCell bc(ac);
	Scale(bc, 2);//scales both
	REQUIRE(Dot(ac, bc)==16*ac.N()*6);
	bc.NewDeep(bc);
	Set(bc, 3);
	REQUIRE(Dot(ac, bc)==12*ac.N()*6);
    }
}
TEST_CASE("sparse", "[sparse]"){
    const int N=6;
    typedef Sparse<Real> Rsp;
    Rsp as(N,N);
    REQUIRE(as.Nzmax()==0);
    REQUIRE(as(2,4)==0);
    as.Set(2,4, 10);
    REQUIRE(as(2,4)==10);
    as.Set(2,4, 4);
    REQUIRE(as(2,4)==4);
    debug(as);
}
