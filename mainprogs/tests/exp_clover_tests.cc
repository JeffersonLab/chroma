#include "chromabase.h"
#include "handle.h"

#include "util/gauge/reunit.h"
#include "gtest/gtest.h"


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
	    }
	}


	void TearDown() {}

	Q u;
};

class ExpClovFixture: public ExpCloverFixtureT<::testing::Test> {};

TEST_F(ExpClovFixture, CheckOp)
{
    QDPIO::cout << "Fred\n";
}

