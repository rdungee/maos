#include <tgmath.h>
#include "catch.h"
#include "../math/array.h"
using std::abs;
#define REQUIRE_EQ(A,B) REQUIRE(abs(A-B)<=abs(A+B)*1e-14)
TEST_CASE("array", "[array]"){
    const Int N=16;
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
    REQUIRE(aa.Nx()==N);
    REQUIRE(aa.Ny()==N);
    REQUIRE(aa.Ldx()==N);
    REQUIRE(aa.N()==N*N);
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
    SECTION("Resize"){
	aa(0,0)=1;
	aa(N-1,1)=2;
	aa.Resize(N*2,N);
	REQUIRE(aa(0,0)==1);
	REQUIRE(aa(1,0)==0);
	REQUIRE(aa(N-1,1)==2);
	REQUIRE(aa(N,1)==0);
    }
    SECTION("NEW"){
	aa.New();
	REQUIRE(aa.N()==N*N);
	aa.New(N*2,N*2);
	REQUIRE(aa.N()==N*N*4);
	REQUIRE(aa(0,0)==0);
    }
}
TEST_CASE("mat", "[mat]"){
    const Int N=4;
    RMat aa(N,N);
    Set(aa, -2); //Set all elements to 2.
    REQUIRE_EQ(aa(1,2),-2);
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
}
