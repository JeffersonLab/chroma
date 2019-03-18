#include "chromabase.h"
#include "handle.h"
#include "linearop.h"
#include "seoprec_linop.h"
#include "eoprec_linop.h"
#include "actions/ferm/fermacts/eoprec_clover_fermact_w.h"
#include "actions/ferm/fermacts/seoprec_clover_fermact_w.h"
#include "actions/ferm/linop/eoprec_clover_linop_w.h"
#include "actions/ferm/linop/seoprec_clover_linop_w.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "util/gauge/reunit.h"
#include "gtest/gtest.h"

#include "io/xml_group_reader.h"
#include "multi_rhs_xml.h"

#include "actions/ferm/invert/syssolver_mrhs_proxy_params.h"
#include "actions/ferm/invert/syssolver_mrhs_twisted_params.h"
#include "actions/ferm/invert/syssolver_mrhs_proxy.h"
#include "actions/ferm/invert/syssolver_mrhs_twisted_proxy.h"
#include "actions/ferm/invert/syssolver_linop_mrhs_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_mrhs_factory.h"

#include "actions/ferm/linop/shifted_linop_w.h"
#include "actions/ferm/linop/multi_twist_linop_w.h"

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
	using S_symm_T =  SymEvenOddPrecLogDetWilsonTypeFermAct<T,P,Q>;
	using LinOpSymm_T = SymEvenOddPrecLogDetLinearOperator<T,P,Q>;
	using LinOpAsymm_T = EvenOddPrecLinearOperator<T,P,Q>;

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
		M_asymm =S_asymm->linOp(state);
		M_symm =S_symm->linOp(state);

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
class MHRSSolverTwistedSeoprecProxyTest : public MultiRHSFixtureT<::testing::TestWithParam<std::string>>{};
class MHRSSolverTwistedEoprecProxyTest : public MultiRHSFixtureT<::testing::TestWithParam<std::string>>{};
class MHRSSolverTwistedSeoprecProxyFactoryTest : public MultiRHSFixtureT<::testing::TestWithParam<std::string>>{};


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

TEST_F(MultiRHSFixture, CheckTwistedParamReader)
{
	std::istringstream inv_param_stream(inv_param_multi_rhs_twisted_proxy_seoprec_cg_xml);
	XMLReader inv_param_xml(inv_param_stream);

	try {
	SysSolverMRHSTwistedParams param(inv_param_xml,"InvertParam");

	QDPIO::cout << "param has BlockSize = " << param.BlockSize << std::endl;
	QDPIO::cout << "param has Twists = {";
	for(int i=0; i < param.Twists.size(); ++i) {
		QDPIO::cout << " " << param.Twists[i];
	}
	QDPIO::cout << " }" << std::endl;

	QDPIO::cout << "param has SubSolver ID =" << param.SubInverterXML.id << std::endl;
	QDPIO::cout << "param has SubSolver Path =" << param.SubInverterXML.path << std::endl;
	QDPIO::cout << "param has SubSolver XML = " << param.SubInverterXML.xml << std::endl;

		ASSERT_EQ( param.BlockSize, 4 );
	}
	catch( const std::string & e ) {
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

TEST_F(MultiRHSFixture, CheckMRHSTwistedOperator)
{

	// Base operator
	Handle< LinOpSymm_T > M_base( S_symm->linOp(state));

	// Setup Twisted Params
	std::istringstream inv_param_stream(inv_param_multi_rhs_twisted_proxy_seoprec_cg_xml);
	XMLReader inv_param_xml(inv_param_stream);
	SysSolverMRHSTwistedParams param(inv_param_xml,"InvertParam");

	// Create a MultiTwist Op
	SymEvenOddPrecLogDetMultiTwistLinOp<T,P,Q> M_multi(M_base, param.Twists);

	// Get the subset
	auto& s = M_multi.subset();

	// Get the Size
	ASSERT_EQ( M_multi.size(), param.BlockSize);
	const int N = M_multi.size();

	multi1d<T> psi(N);
	multi1d<T> chi(N);


	for(int i=0;i < N; ++i) {
		psi[i] = zero;
		gaussian(psi[i],s);
	}

	(M_multi)(chi,psi,PLUS);

	for(int i=0; i < N; ++i) {
		T tmp = zero;
		SymEvenOddPrecLogDetTwistedShiftedLinOp<T,P,Q> theShiftedOp((*M_base), param.Twists[i]);

		(theShiftedOp)( tmp, psi[i], PLUS );
		tmp[ s ] -= chi[i];
		Double diff = sqrt(norm2(tmp,s));
		Double diff_chi = sqrt(norm2(chi[i],s));
		Double rel_diff = diff/diff_chi;
		QDPIO::cout << "i= "<<i << " Diff= " << diff << " Rel diff= " << rel_diff << std::endl;
		ASSERT_LT( toDouble(rel_diff), 1.0e-8);
	}
 }

TEST_F(MultiRHSFixture, CheckAsymmMRHSTwistedOperator)
{

	// Base operator
	Handle< LinOpAsymm_T > M_base( S_asymm->linOp(state));

	// Setup Twisted Params
	std::istringstream inv_param_stream(inv_param_multi_rhs_twisted_proxy_eoprec_cg_xml);
	XMLReader inv_param_xml(inv_param_stream);
	SysSolverMRHSTwistedParams param(inv_param_xml,"InvertParam");

	// Create a MultiTwist Op
	EvenOddPrecMultiTwistLinOp<T,P,Q> M_multi(M_base, param.Twists);

	// Get the subset
	auto& s = M_multi.subset();

	// Get the Size
	ASSERT_EQ( M_multi.size(), param.BlockSize);
	const int N = M_multi.size();

	multi1d<T> psi(N);
	multi1d<T> chi(N);


	for(int i=0;i < N; ++i) {
		psi[i] = zero;
		gaussian(psi[i],s);
	}

	(M_multi)(chi,psi,PLUS);

	for(int i=0; i < N; ++i) {
		T tmp = zero;
		EvenOddPrecTwistedShiftedLinOp<T,P,Q> theShiftedOp((*M_base), param.Twists[i]);

		(theShiftedOp)( tmp, psi[i], PLUS );
		tmp[ s ] -= chi[i];
		Double diff = sqrt(norm2(tmp,s));
		Double diff_chi = sqrt(norm2(chi[i],s));
		Double rel_diff = diff/diff_chi;
		QDPIO::cout << "i= "<<i << " Diff= " << diff << " Rel diff= " << rel_diff << std::endl;
		ASSERT_LT( toDouble(rel_diff), 1.0e-8);
	}
 }


TEST_F(MultiRHSFixture, CheckLinOpMRHSSysSolverProxy)
{
	std::istringstream inv_param_stream(inv_param_multi_rhs_proxy_cg_xml);
	XMLReader inv_param_xml(inv_param_stream);

	SysSolverMRHSProxyParams param(inv_param_xml,"InvertParam");

	// Smart pointers including our handle cannot cast covariantly.\
	// This is a Handle<> mimic of static_cast_ptr
	LinOpMRHSSysSolverProxy<T,P,Q> the_solver(param, S_symm.cast_static<FermAct4D<T,P,Q>>(), state);

	ASSERT_EQ( the_solver.size(), param.BlockSize);

}

TEST_P(MHRSSolverProxyTest, CheckLinOpMRHSProxyWorks)
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

TEST_P(MHRSSolverProxyTest, CheckMdagMMRHSProxyWorks)
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

TEST_F(MultiRHSFixture, CheckLinOpMRHSProxyFectoryCreation)
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

TEST_F(MultiRHSFixture, CheckMdagMMRHSProxyFectoryCreation)
{
	std::istringstream inv_param_stream(inv_param_multi_rhs_proxy_cg_xml);
	XMLReader inv_param_xml(inv_param_stream);

	Handle<MdagMMRHSSystemSolver<LatticeFermion>> the_solver =
			TheMdagMFermMRHSSystemSolverFactory::Instance().createObject("MULTI_RHS_PROXY_INVERTER",
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

TEST_P(MHRSSolverTwistedSeoprecProxyTest, CheckLinOpMRHSTwistedProxySolverWorks)
{
	std::istringstream inv_param_stream(GetParam());
	XMLReader inv_param_xml(inv_param_stream);

	SysSolverMRHSTwistedParams param(inv_param_xml,"InvertParam");

	// Smart pointers including our handle cannot cast covariantly.\
	// This is a Handle<> mimic of static_cast_ptr
	SymEvenOddPrecLogDetLinOpMRHSSysSolverTwistedProxy<T,P,Q>
		the_solver(param, S_symm.cast_static<FermAct4D<T,P,Q>>(), state);

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
		SymEvenOddPrecLogDetTwistedShiftedLinOp<T,P,Q> M_twisted((*M_symm), param.Twists[i]);
		(M_twisted)( tmp, psi[i], PLUS );
		tmp[ s ] -= chi[i];
		Double diff = sqrt(norm2(tmp,s));
		Double diff_chi = sqrt(norm2(chi[i],s));
		Double rel_diff = diff/diff_chi;
		QDPIO::cout << "i= "<<i << " Diff= " << diff << " Rel diff= " << rel_diff << std::endl;
		ASSERT_LT( toDouble(rel_diff), 1.0e-8);
	}
}

TEST_P(MHRSSolverTwistedSeoprecProxyTest, CheckMdagMMRHSTwistedProxySolverWorks)
{
	std::istringstream inv_param_stream(GetParam());
	XMLReader inv_param_xml(inv_param_stream);

	SysSolverMRHSTwistedParams param(inv_param_xml,"InvertParam");

	// Smart pointers including our handle cannot cast covariantly.\
	// This is a Handle<> mimic of static_cast_ptr
	SymEvenOddPrecLogDetMdagMMRHSSysSolverTwistedProxy<T,P,Q>
		the_solver(param, S_symm.cast_static<FermAct4D<T,P,Q>>(), state);

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
		SymEvenOddPrecLogDetTwistedShiftedLinOp<T,P,Q> M_twisted((*M_symm), param.Twists[i]);
		(M_twisted)( tmp, psi[i], PLUS );
		(M_twisted)( tmp2, tmp, MINUS );
		tmp2[ s ] -= chi[i];
		Double diff = sqrt(norm2(tmp2,s));
		Double diff_chi = sqrt(norm2(chi[i],s));
		Double rel_diff = diff/diff_chi;
		QDPIO::cout << "i= "<<i << " Diff= " << diff << " Rel diff= " << rel_diff << std::endl;
		ASSERT_LT( toDouble(rel_diff), 1.0e-8);
	}
}
#if 1
INSTANTIATE_TEST_CASE_P(MRHSSyssolverTwistedSeoprecProxy,
						MHRSSolverTwistedSeoprecProxyTest,
                        ::testing::Values(inv_param_multi_rhs_twisted_proxy_seoprec_cg_xml));
#endif

TEST_P(MHRSSolverTwistedEoprecProxyTest, CheckLinOpMRHSTwistedProxySolverWorks)
{
	std::istringstream inv_param_stream(GetParam());
	XMLReader inv_param_xml(inv_param_stream);

	SysSolverMRHSTwistedParams param(inv_param_xml,"InvertParam");

	// Smart pointers including our handle cannot cast covariantly.\
	// This is a Handle<> mimic of static_cast_ptr
	EvenOddPrecLinOpMRHSSysSolverTwistedProxy<T,P,Q>
		the_solver(param, S_asymm.cast_static<FermAct4D<T,P,Q>>(), state);

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
		EvenOddPrecTwistedShiftedLinOp<T,P,Q> M_twisted((*M_asymm), param.Twists[i]);
		(M_twisted)( tmp, psi[i], PLUS );
		tmp[ s ] -= chi[i];
		Double diff = sqrt(norm2(tmp,s));
		Double diff_chi = sqrt(norm2(chi[i],s));
		Double rel_diff = diff/diff_chi;
		QDPIO::cout << "i= "<<i << " Diff= " << diff << " Rel diff= " << rel_diff << std::endl;
		ASSERT_LT( toDouble(rel_diff), 1.0e-8);
	}
}

TEST_P(MHRSSolverTwistedEoprecProxyTest, CheckMdagMMRHSTwistedProxySolverWorks)
{
	std::istringstream inv_param_stream(GetParam());
	XMLReader inv_param_xml(inv_param_stream);

	SysSolverMRHSTwistedParams param(inv_param_xml,"InvertParam");

	// Smart pointers including our handle cannot cast covariantly.\
	// This is a Handle<> mimic of static_cast_ptr
	EvenOddPrecMdagMMRHSSysSolverTwistedProxy<T,P,Q>
		the_solver(param, S_asymm.cast_static<FermAct4D<T,P,Q>>(), state);

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
		EvenOddPrecTwistedShiftedLinOp<T,P,Q> M_twisted((*M_asymm), param.Twists[i]);
		(M_twisted)( tmp, psi[i], PLUS );
		(M_twisted)( tmp2, tmp, MINUS );
		tmp2[ s ] -= chi[i];
		Double diff = sqrt(norm2(tmp2,s));
		Double diff_chi = sqrt(norm2(chi[i],s));
		Double rel_diff = diff/diff_chi;
		QDPIO::cout << "i= "<<i << " Diff= " << diff << " Rel diff= " << rel_diff << std::endl;
		ASSERT_LT( toDouble(rel_diff), 1.0e-8);
	}
}
#if 1
INSTANTIATE_TEST_CASE_P(MRHSSyssolverTwistedEoprecProxy,
						MHRSSolverTwistedEoprecProxyTest,
                        ::testing::Values(inv_param_multi_rhs_twisted_proxy_eoprec_cg_xml));
#endif

TEST_P(MHRSSolverTwistedSeoprecProxyFactoryTest, CheckLinOpMRHSTwistedProxyFactoryCreate)
{
	std::istringstream inv_param_stream(GetParam());
	XMLReader inv_param_xml(inv_param_stream);

	// To grab the twists for checking
	SysSolverMRHSTwistedParams param(inv_param_xml,"InvertParam");

	Handle<LinOpMRHSSystemSolver<LatticeFermion>> the_solver(
			TheLinOpFermMRHSSystemSolverFactory::Instance().createObject("MULTI_RHS_SEOPREC_TWISTED_PROXY_INVERTER",
					inv_param_xml, "InvertParam", S_symm.cast_static<FermAct4D<T,P,Q>>(), state)
	);


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
	SystemSolverResultsMRHS_t res = (*the_solver)(psi, chi);

	ASSERT_EQ( res.resid.size(), N);
	ASSERT_EQ( res.n_count.size(), N);

	for(int i=0; i < N; ++i) {
		T tmp = zero;
		SymEvenOddPrecLogDetTwistedShiftedLinOp<T,P,Q> M_twisted((*M_symm), param.Twists[i]);
		(M_twisted)( tmp, psi[i], PLUS );
		tmp[ s ] -= chi[i];
		Double diff = sqrt(norm2(tmp,s));
		Double diff_chi = sqrt(norm2(chi[i],s));
		Double rel_diff = diff/diff_chi;
		QDPIO::cout << "i= "<<i << " Diff= " << diff << " Rel diff= " << rel_diff << std::endl;
		ASSERT_LT( toDouble(rel_diff), 1.0e-8);
	}
}

TEST_P(MHRSSolverTwistedSeoprecProxyFactoryTest, CheckMdagMMRHSTwistedProxyFactoryCreate)
{
	std::istringstream inv_param_stream(GetParam());
	XMLReader inv_param_xml(inv_param_stream);
	// To grab the twists for checking
	SysSolverMRHSTwistedParams param(inv_param_xml,"InvertParam");

	// Smart pointers including our handle cannot cast covariantly.\
	// This is a Handle<> mimic of static_cast_ptr
	Handle<MdagMMRHSSystemSolver<LatticeFermion>> the_solver(
			TheMdagMFermMRHSSystemSolverFactory::Instance().createObject("MULTI_RHS_SEOPREC_TWISTED_PROXY_INVERTER",
					inv_param_xml, "InvertParam", S_symm.cast_static<FermAct4D<T,P,Q>>(), state)
	);

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
	SystemSolverResultsMRHS_t res = (*the_solver)(psi, chi);

	ASSERT_EQ( res.resid.size(), N);
	ASSERT_EQ( res.n_count.size(), N);

	for(int i=0; i < N; ++i) {
		T tmp = zero;
		T tmp2 = zero;
		SymEvenOddPrecLogDetTwistedShiftedLinOp<T,P,Q> M_twisted((*M_symm), param.Twists[i]);
		(M_twisted)( tmp, psi[i], PLUS );
		(M_twisted)( tmp2, tmp, MINUS );
		tmp2[ s ] -= chi[i];
		Double diff = sqrt(norm2(tmp2,s));
		Double diff_chi = sqrt(norm2(chi[i],s));
		Double rel_diff = diff/diff_chi;
		QDPIO::cout << "i= "<<i << " Diff= " << diff << " Rel diff= " << rel_diff << std::endl;
		ASSERT_LT( toDouble(rel_diff), 1.0e-8);
	}
}
#if 1
INSTANTIATE_TEST_CASE_P(MRHSSyssolverTwistedSeoprecProxy,
						MHRSSolverTwistedSeoprecProxyFactoryTest,
                        ::testing::Values(inv_param_multi_rhs_twisted_proxy_seoprec_cg_xml));
#endif

