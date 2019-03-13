#include "chromabase.h"
#include "handle.h"
#include "seoprec_linop.h"
#include "eoprec_linop.h"
#include "eoprec_wilstype_fermact_w.h"
#include "seoprec_wilstype_fermact_w.h"
#include "actions/ferm/linop/eoprec_clover_linop_w.h"
#include "actions/ferm/linop/seoprec_clover_linop_w.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "util/gauge/reunit.h"
#include "gtest/gtest.h"

#include "io/xml_group_reader.h"
#include "multi_rhs_xml.h"

#include "actions/ferm/invert/syssolver_mrhs_proxy_params.h"
using namespace Chroma;
using namespace QDP;
using namespace MultiRHSTesting;

template<typename TestType>
class MultiRHSFixtureT : public TestType {
public:
	using T = LatticeFermion;
	using Q = multi1d<LatticeColorMatrix>;
	using P = multi1d<LatticeColorMatrix>;

	using S_asymm_T = EvenOddPrecWilsonTypeFermAct<T,P,Q>;
	using S_symm_T =  SymEvenOddPrecWilsonTypeFermAct<T,P,Q>;
	using LinOpSymm_T = SymEvenOddPrecCloverLinOp;
	using LinOpAsymm_T = EvenOddPrecCloverLinOp;

	void SetUp() {
	  u.resize(Nd);
	    for(int mu=0; mu < Nd; ++mu) {
	    	gaussian(u[mu]);
	    	reunit(u[mu]);
	    }
	    std::istringstream input_asymm(fermact_xml_asymm);
	    std::istringstream input_symm(fermact_xml_symm);

	    XMLReader xml_in_asymm(input_asymm);
	    XMLReader xml_in_symm(input_symm);


	    S_asymm = dynamic_cast<S_asymm_T*>(TheFermionActionFactory::Instance().createObject("CLOVER",
	    											   xml_in_asymm,
	    											   "FermionAction"));

	    S_symm = dynamic_cast<S_symm_T*>(TheFermionActionFactory::Instance().createObject("SEOPREC_CLOVER",
	        											   xml_in_symm,
	        											   "FermionAction"));
	}


	void TearDown() {}

	Q u;
	Handle<S_symm_T> S_symm;
	Handle<S_asymm_T> S_asymm;
};

class MultiRHSFixture : public MultiRHSFixtureT<::testing::Test> {};
// class QPropTest : public SymmFixtureT<::testing::TestWithParam<std::string>>{};

TEST_F(MultiRHSFixture, CheckParamReader)
{
	std::istringstream inv_param_stream(inv_param_multi_rhs_proxy_xml);
	XMLReader inv_param_xml(inv_param_stream);

	try {
	SysSolverMRHSProxyParams param(inv_param_xml,"InvertParam");

	QDPIO::cout << "param has BlockSize = " << param.BlockSize << std::endl;
	QDPIO::cout << "param has SubSolver ID =" << param.SubInverterXML.id << std::endl;
	QDPIO::cout << "param has SubSolver Path =" << param.SubInverterXML.path << std::endl;
	QDPIO::cout << "param has SubSolver XML = " << param.SubInverterXML.xml << std::endl;

		ASSERT_EQ( param.BlockSize, 4 );
	}
	catch( const std::string e ) {
		QDPIO::cout << "Caught exception " << e << std::endl;
		FAIL();
	}

 }

