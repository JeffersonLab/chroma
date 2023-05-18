#include "chromabase.h"
#include "handle.h"
#include <cmath>

#include "util/gauge/reunit.h"
#include "gtest/gtest.h"

#include "actions/ferm/fermstates/simple_fermstate.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/linop/clover_term_qdp_w.h"
#include "actions/ferm/linop/exp_clover_term_qdp_w.h"

using namespace Chroma;
using namespace QDP;

template<typename TestType>
class ExpCloverFixtureT : public TestType {
public:
	using T = LatticeFermion;
	using Q = multi1d<LatticeColorMatrix>;
	using P = multi1d<LatticeColorMatrix>;

	void SetUp() 
    {
	  	u.resize(Nd);
	    for(int mu=0; mu < Nd; ++mu) {
	    	gaussian(u[mu]);
	    	reunit(u[mu]);
			// u[mu] = 1;
	    }
	  	multi1d<int> bcs(4);
		bcs[0] = bcs[1] = bcs[2] = 1;
		bcs[3] = -1;
		simpleFermState = new SimpleFermState<T,P,Q>(bcs, u);

		CloverFermActParams p;
		p.Mass = Real(Mass);
		p.clovCoeffR = 1;
		p.clovCoeffT = 1;
		p.u0 = 1;
		p.anisoParam.anisoP = false;
		p.anisoParam.t_dir = 3;
		p.anisoParam.xi_0  = Real(1);
		p.twisted_m = 0;
		p.twisted_m_usedP = false;

		clov.create( simpleFermState, p );
		eclov.create( simpleFermState, p );
	}

	static constexpr double Mass = 0.1;

	void TearDown() {}

	Q u;

	Handle< FermState<T,P,Q> > simpleFermState;
	QDPCloverTerm clov;
	QDPExpCloverTerm<> eclov;
	
};

class ExpClovFixture: public ExpCloverFixtureT<::testing::Test> {};

TEST_F(ExpClovFixture, CheckOp)
{
    LatticeFermion src, res, res_exp, dummy, diff;
	gaussian(src);
	res = zero;
	res_exp = zero;
	

	// We will be going for an exponential in the end:
	// 
	// so: 
	//  exp(x) = (diag mass)[ 1 + E + 1/2 E^2 + .... ]
	// 
	//  First test: (diag mass)[ 1 + E ] = regular clover term.
	for( int cb=0; cb < 2; ++cb ) {
		clov.apply( res, src, PLUS, cb);
		eclov.applyUnexp( dummy, src, PLUS, cb);
		res_exp[rb[cb]] = src + dummy;
		res_exp[rb[cb]] *= Real( Nd + Mass );

		diff[rb[cb]] = res_exp - res;
		Double normdiff = sqrt( norm2(diff,rb[cb])/norm2(src,rb[cb]) );
		QDPIO::cout << "Diff (" << cb << ") = " << normdiff << "\n";

		ASSERT_LT( toDouble(normdiff), 1.0e-14);
	}
}

TEST_F(ExpClovFixture, CheckRefOp)
{
    LatticeFermion src, res, diff;
	src = zero;
	src.elem(0).elem(0).elem(0).real() = 1;
	eclov.fillRefDiag(0.5);

	res = zero;	
	eclov.applyRef(res,src,PLUS,2);

	double r2 = res.elem(0).elem(0).elem(0).real();

	res = zero;	
	eclov.applyRef(res,src,PLUS,10);

	double r10 = res.elem(0).elem(0).elem(0).real();

	res = zero;	
	eclov.applyRef(res,src,PLUS,13);

	double r13 = res.elem(0).elem(0).elem(0).real();

	double ref = exp(0.5);
	QDPIO::cout << "N=2   res="<< r2 << " exp(0.5)=" << ref << " abs. diff=" << abs(ref-r2) << "\n";
	QDPIO::cout << "N=10   res="<< r10 << " exp(0.5)=" << ref << " abs. diff=" << abs(ref-r10) << "\n";
	QDPIO::cout << "N=13   res="<< r13 << " exp(0.5)=" << ref << " abs. diff=" << abs(ref-r13) << "\n";
	
	ASSERT_LT( abs(ref-r13), 1.0e-15);
}



TEST_F(ExpClovFixture, CheckApplyPower0)
{
    LatticeFermion src, res, diff;
	gaussian(src);
	res = zero;
	
	
	// A^0 = I 
	for( int cb=0; cb < 2; ++cb ) {
		
		eclov.applyPower(res,src, PLUS, cb, 0);
	
		diff[rb[cb]] = res - src;
		Double normdiff = sqrt( norm2(diff,rb[cb])/norm2(src,rb[cb]) );
		QDPIO::cout << "Diff (" << cb << ") = " << normdiff << "\n";

		ASSERT_LT( toDouble(normdiff), 1.0e-14);
	}
}

TEST_F(ExpClovFixture, CheckApplyPower1)
{
    LatticeFermion src, res, res2, diff;
	gaussian(src);
	res = zero;
	res2 = zero;
	
	for( int cb=0; cb < 2; ++cb ) {
		eclov.applyUnexp(res, src, PLUS, cb);
		eclov.applyPower(res2, src, PLUS, cb, 1);

		diff[rb[cb]] = res2 - res;
		Double normdiff = sqrt( norm2(diff,rb[cb])/norm2(src,rb[cb]) );
		QDPIO::cout << "Diff (" << cb << ") = " << normdiff << "\n";

		ASSERT_LT( toDouble(normdiff), 1.0e-14);
	}
}

TEST_F(ExpClovFixture, CheckApplyPower2)
{
    LatticeFermion src, res, res2, diff;
	gaussian(src);
	res = zero;
	res2 = zero;
	
	for( int cb=0; cb < 2; ++cb ) {
		eclov.applyUnexp(res, src, PLUS, cb);
		eclov.applyUnexp(res2, res, PLUS, cb);
		eclov.applyPower(res, src, PLUS, cb, 2);

		diff[rb[cb]] = res2 - res;
		Double normdiff = sqrt( norm2(diff,rb[cb])/norm2(src,rb[cb]) );
		QDPIO::cout << "Diff (" << cb << ") = " << normdiff << "\n";

		ASSERT_LT( toDouble(normdiff), 1.0e-14);
	}
}

TEST_F(ExpClovFixture, CheckApplyPower3)
{
    LatticeFermion src, res, res2, diff;
	gaussian(src);
	res = zero;
	res2 = zero;
	
	for( int cb=0; cb < 2; ++cb ) {
		eclov.applyUnexp(res, src, PLUS, cb);
		eclov.applyUnexp(res2, res, PLUS, cb);
		eclov.applyUnexp(res, res2, PLUS, cb);
		eclov.applyPower(res2, src, PLUS, cb, 3);

		diff[rb[cb]] = res2 - res;
		Double normdiff = sqrt( norm2(diff,rb[cb])/norm2(src,rb[cb]) );
		QDPIO::cout << "Diff (" << cb << ") = " << normdiff << "\n";

		ASSERT_LT( toDouble(normdiff), 1.0e-14);
	}
}

TEST_F(ExpClovFixture, CheckApplyPower4)
{
    LatticeFermion src, res, res2, diff;
	gaussian(src);
	res = zero;
	res2 = zero;
	
	for( int cb=0; cb < 2; ++cb ) {
		eclov.applyUnexp(res, src, PLUS, cb);
		eclov.applyUnexp(res2, res, PLUS, cb);
		eclov.applyUnexp(res, res2, PLUS, cb);
		eclov.applyUnexp(res2, res, PLUS, cb);
		eclov.applyPower(res, src, PLUS, cb, 4);

		diff[rb[cb]] = res2 - res;
		Double normdiff = sqrt( norm2(diff,rb[cb])/norm2(src,rb[cb]) );
		QDPIO::cout << "Diff (" << cb << ") = " << normdiff << "\n";

		ASSERT_LT( toDouble(normdiff), 1.0e-14);
	}
}

TEST_F(ExpClovFixture, CheckApplyPower5)
{
    LatticeFermion src, res, res2, diff;
	gaussian(src);
	res = zero;
	res2 = zero;
	
	for( int cb=0; cb < 2; ++cb ) {
		eclov.applyUnexp(res, src, PLUS, cb);
		eclov.applyUnexp(res2, res, PLUS, cb);
		eclov.applyUnexp(res, res2, PLUS, cb);
		eclov.applyUnexp(res2, res, PLUS, cb);
		eclov.applyUnexp(res, res2, PLUS, cb);
		eclov.applyPower(res2, src, PLUS, cb, 5);

		diff[rb[cb]] = res2 - res;
		Double normdiff = sqrt( norm2(diff,rb[cb])/norm2(src,rb[cb]) );
		QDPIO::cout << "Diff (" << cb << ") = " << normdiff << "\n";

		ASSERT_LT( toDouble(normdiff), 1.0e-14);
	}
}

TEST_F(ExpClovFixture, CheckApplyPower7)
{
    LatticeFermion src, res, res2, diff;
	gaussian(src);
	res = zero;
	res2 = zero;
	
	for( int cb=0; cb < 2; ++cb ) {
		eclov.applyUnexp(res, src, PLUS, cb);
		eclov.applyUnexp(res2, res, PLUS, cb);
		eclov.applyUnexp(res, res2, PLUS, cb);
		eclov.applyUnexp(res2, res, PLUS, cb);
		eclov.applyUnexp(res, res2, PLUS, cb);
		eclov.applyUnexp(res2, res, PLUS, cb);
		eclov.applyUnexp(res, res2, PLUS, cb);
		eclov.applyPower(res2, src, PLUS, cb, 7);

		diff[rb[cb]] = res2 - res;
		Double normdiff = sqrt( norm2(diff,rb[cb])/norm2(src,rb[cb]) );
		QDPIO::cout << "Diff (" << cb << ") = " << normdiff << "\n";

		ASSERT_LT( toDouble(normdiff), 1.0e-14);
	}
}
TEST_F(ExpClovFixture, CheckExp)
{
    LatticeFermion src, res, res2, dummy, diff;
	gaussian(src);
	res = zero;
	res2 = zero;
	

	// We will be going for an exponential in the end:
	// 
	// so: 
	//  exp(x) = (diag mass)[ 1 + E + 1/2 E^2 + .... ]
	// 
	//  First test: (diag mass)[ 1 + E ] = regular clover term.
	eclov.applyRef(res, src, PLUS,28);
	for( int cb=0; cb < 2; ++cb ) {
		eclov.apply(res2,src, PLUS, cb);
	}
		
	diff = res2 - res;
	Double normdiff = sqrt( norm2(diff)/norm2(src) );
	QDPIO::cout << "Diff  = " << normdiff << "\n";

	ASSERT_LT( toDouble(normdiff), 1.0e-14);
	
}

TEST_F(ExpClovFixture, CheckExpInv)
{
    LatticeFermion src, res, res2, dummy, diff;
	gaussian(src);
	res = zero;
	res2 = zero;
	
	LatticeReal tr_Minv;
	
	//eclov.choles(0);
	//eclov.choles(1);

	// We will be going for an exponential in the end:
	// 
	// so: 
	//  exp(x) = (diag mass)[ 1 + E + 1/2 E^2 + .... ]
	// 
	//  First test: (diag mass)[ 1 + E ] = regular clover term.
	for( int cb=0; cb < 2; ++cb ) {
		eclov.apply(res,src, PLUS, cb);
		eclov.applyInv(res2,res, PLUS, cb);
	}
		
	diff = res2 - src;
	Double normdiff = sqrt( norm2(diff)/norm2(src) );
	QDPIO::cout << "Diff  = " << normdiff << "\n";

	ASSERT_LT( toDouble(normdiff), 1.0e-14);
	
}
