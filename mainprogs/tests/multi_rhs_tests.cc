#include "chromabase.h"
#include "handle.h"
#include "linearop.h"
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
#include "actions/ferm/invert/syssolver_mrhs_proxy.h"
#include "actions/ferm/invert/syssolver_linop_mrhs_factory.h"
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

		state = S_asymm->createState(u);
		M_asymm =dynamic_cast<LinOpAsymm_T *>(S_asymm->linOp(state));
		M_symm =dynamic_cast<LinOpSymm_T *>(S_symm->linOp(state));

	}


	void TearDown() {}

	Q u;
	Handle<S_symm_T> S_symm;
	Handle<S_asymm_T> S_asymm;
	Handle<FermState<T,P,Q> > state;
	Handle<LinOpAsymm_T> M_asymm;
	Handle<LinOpSymm_T> M_symm;
};

class MultiRHSFixture : public MultiRHSFixtureT<::testing::Test> {};
class MHRSSolverProxyTest : public MultiRHSFixtureT<::testing::TestWithParam<std::string>>{};

TEST_F(MultiRHSFixture, CheckParamReader)
{
	std::istringstream inv_param_stream(inv_param_multi_rhs_proxy_cg_xml);
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

TEST_F(MultiRHSFixture, CheckMRHSOperator)
{
	int num_rhs =4;
	Handle< LinearOperatorArray<T> > M_multi(S_symm->linOpMRHS(state,num_rhs));
	auto& s = M_multi->subset();
	ASSERT_EQ( M_multi->size(), 4);

	multi1d<T> psi(num_rhs);
	multi1d<T> chi(num_rhs);


	for(int i=0;i < num_rhs; ++i) gaussian(psi[i]);
	(*M_multi)(chi,psi,PLUS);

	for(int i=0; i < num_rhs; ++i) {
		T tmp = zero;
		(*M_symm)( tmp, psi[i], PLUS );
		tmp[ s ] -= chi[i];
		Double diff = sqrt(norm2(tmp,s));
		Double diff_chi = sqrt(norm2(chi[i],s));
		Double rel_diff = diff/diff_chi;
		QDPIO::cout << "i= "<<i << " Diff= " << diff << " Rel diff= " << rel_diff << std::endl;
		ASSERT_LT( toDouble(rel_diff), 1.0e-8);
	}
 }

TEST_F(MultiRHSFixture, CheckLinOpMRHSWProxy)
{
	std::istringstream inv_param_stream(inv_param_multi_rhs_proxy_cg_xml);
	XMLReader inv_param_xml(inv_param_stream);

	SysSolverMRHSProxyParams param(inv_param_xml,"InvertParam");

	// Smart pointers including our handle cannot cast covariantly.\
	// This is a Handle<> mimic of static_cast_ptr
	LinOpMRHSSysSolverProxy<T,P,Q> the_solver(param, S_symm.cast_static<FermAct4D<T,P,Q>>(), state);

	ASSERT_EQ( the_solver.size(), param.BlockSize);

}

TEST_P(MHRSSolverProxyTest, CheckLinOpMRHSWProxyWorks)
{
	std::istringstream inv_param_stream(GetParam());
	XMLReader inv_param_xml(inv_param_stream);

	SysSolverMRHSProxyParams param(inv_param_xml,"InvertParam");

	// Smart pointers including our handle cannot cast covariantly.\
	// This is a Handle<> mimic of static_cast_ptr
	LinOpMRHSSysSolverProxy<T,P,Q> the_solver(param, S_symm.cast_static<FermAct4D<T,P,Q>>(), state);

	const int N = the_solver.size();
	const Subset& s = the_solver.subset();
	multi1d<T> chi(N);
	multi1d<T> psi(N);

	for(int i=0; i < N; ++i) {
		chi[i]=zero;
		psi[i]=zero;
		gaussian(chi[i], s);
	}

	// Solve all poles at once
	SystemSolverResultsMRHS_t res = the_solver(psi, chi);

	ASSERT_EQ( res.resid.size(), the_solver.size());
	ASSERT_EQ( res.n_count.size(), the_solver.size());

	for(int i=0; i < N; ++i) {
		T tmp = zero;
		(*M_symm)( tmp, psi[i], PLUS );
		tmp[ s ] -= chi[i];
		Double diff = sqrt(norm2(tmp,s));
		Double diff_chi = sqrt(norm2(chi[i],s));
		Double rel_diff = diff/diff_chi;
		QDPIO::cout << "i= "<<i << " Diff= " << diff << " Rel diff= " << rel_diff << std::endl;
		ASSERT_LT( toDouble(rel_diff), 1.0e-8);
	}
}

TEST_P(MHRSSolverProxyTest, CheckMdagMMRHSWProxyWorks)
{
	std::istringstream inv_param_stream(GetParam());
	XMLReader inv_param_xml(inv_param_stream);

	SysSolverMRHSProxyParams param(inv_param_xml,"InvertParam");

	// Smart pointers including our handle cannot cast covariantly.\
	// This is a Handle<> mimic of static_cast_ptr
	MdagMMRHSSysSolverProxy<T,P,Q> the_solver(param, S_symm.cast_static<FermAct4D<T,P,Q>>(), state);

	const int N = the_solver.size();
	const Subset& s = the_solver.subset();
	multi1d<T> chi(N);
	multi1d<T> psi(N);

	for(int i=0; i < N; ++i) {
		chi[i]=zero;
		psi[i]=zero;
		gaussian(chi[i], s);
	}

	// Solve all poles at once
	SystemSolverResultsMRHS_t res = the_solver(psi, chi);

	ASSERT_EQ( res.resid.size(), the_solver.size());
	ASSERT_EQ( res.n_count.size(), the_solver.size());

	for(int i=0; i < N; ++i) {
		T tmp = zero;
		T tmp2 = zero;
		(*M_symm)( tmp, psi[i], PLUS );
		(*M_symm)( tmp2, tmp, MINUS );
		tmp2[ s ] -= chi[i];
		Double diff = sqrt(norm2(tmp2,s));
		Double diff_chi = sqrt(norm2(chi[i],s));
		Double rel_diff = diff/diff_chi;
		QDPIO::cout << "i= "<<i << " Diff= " << diff << " Rel diff= " << rel_diff << std::endl;
		ASSERT_LT( toDouble(rel_diff), 1.0e-8);
	}
}

INSTANTIATE_TEST_CASE_P(MRHSSyssolverProxy,
						MHRSSolverProxyTest,
                        ::testing::Values(inv_param_multi_rhs_proxy_cg_xml,
                        		inv_param_multi_rhs_proxy_bicgstab_xml));

#ifdef BUILD_QUDA
INSTANTIATE_TEST_CASE_P(MRHSQUDASyssolverProxy,
						MHRSSolverProxyTest,
                        ::testing::Values(inv_param_multi_rhs_proxy_quda_bicgstab_xml,
                        		inv_param_multi_rhs_proxy_quda_multigrid_xml));
#endif

TEST_F(MultiRHSFixture, CheckLinOpMRHSProxyFectoryCreateion)
{
	std::istringstream inv_param_stream(inv_param_multi_rhs_proxy_cg_xml);
	XMLReader inv_param_xml(inv_param_stream);

	Handle<LinOpMRHSSystemSolver<LatticeFermion>> the_solver =
			TheLinOpFermMRHSSystemSolverFactory::Instance().createObject("MULTI_RHS_PROXY_INVERTER",
					inv_param_xml,"InvertParam",S_symm.cast_static<FermAct4D<T,P,Q>>(), state);

	const int N = the_solver->size();
	const Subset& s = the_solver->subset();
	multi1d<T> chi(N);
	multi1d<T> psi(N);

	for(int i=0; i < N; ++i) {
		chi[i]=zero;
		psi[i]=zero;
		gaussian(chi[i], s);
	}

	// Solve all poles at once
	SystemSolverResultsMRHS_t res =(*the_solver)(psi, chi);

	ASSERT_EQ( res.resid.size(), N);
	ASSERT_EQ( res.n_count.size(), N);

	for(int i=0; i < N; ++i) {
		T tmp = zero;
		(*M_symm)( tmp, psi[i], PLUS );
		tmp[ s ] -= chi[i];
		Double diff = sqrt(norm2(tmp,s));
		Double diff_chi = sqrt(norm2(chi[i],s));
		Double rel_diff = diff/diff_chi;
		QDPIO::cout << "i= "<<i << " Diff= " << diff << " Rel diff= " << rel_diff << std::endl;
		ASSERT_LT( toDouble(rel_diff), 1.0e-8);
	}
}
