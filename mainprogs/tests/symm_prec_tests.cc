#include "chromabase.h"

#include "handle.h"

#include "seoprec_linop.h"
#include "eoprec_linop.h"
#include "eoprec_wilstype_fermact_w.h"
#include "seoprec_wilstype_fermact_w.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "util/gauge/reunit.h"
#include "gtest/gtest.h"

#include "io/xml_group_reader.h"
#include "./symm_prec_xml.h"

#ifdef BUILD_QUDA
#include "quda.h"
#include "update/molecdyn/predictor/quda_predictor.h"
#endif

using namespace Chroma;
using namespace QDP;
using namespace SymmPrecTesting;

template<typename TestType>
class SymmFixtureT : public TestType {
public:
	using T = LatticeFermion;
	using Q = multi1d<LatticeColorMatrix>;
	using P = multi1d<LatticeColorMatrix>;

	using S_asymm_T = EvenOddPrecWilsonTypeFermAct<T,P,Q>;
	using S_symm_T =  SymEvenOddPrecWilsonTypeFermAct<T,P,Q>;
	using LinOpSymm_T = SymEvenOddPrecLinearOperator<T,P,Q>;
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

	    M_asymm = S_asymm->linOp(state);
	    M_symm = S_symm->linOp(state);
	}

	void TearDown() {}

	Q u;
	Handle<S_symm_T> S_symm;
	Handle<S_asymm_T> S_asymm;
	Handle<FermState<T,P,Q> > state;
	Handle<LinOpAsymm_T> M_asymm;
	Handle<LinOpSymm_T> M_symm;
};

class SymmFixture : public SymmFixtureT<::testing::Test> {};
class QPropTest : public SymmFixtureT<::testing::TestWithParam<std::string>>{};

TEST_F(SymmFixture, CheckOp)
{
	// Verify A_oo M_symm = M_asymm
	LatticeFermion x=zero;
	gaussian(x, rb[1]);

	LatticeFermion t1 = zero;


	LatticeFermion AM_symm_x = zero;
	LatticeFermion M_asymm_x = zero;

	{
		QDPIO::cout << "Op: " << std::endl;
		// AM_symm_x = A_oo t_1 = A_oo M_symm x
		(*M_symm)(t1,x,PLUS);
		(*M_symm).unprecOddOddLinOp(AM_symm_x,t1,PLUS);


		(*M_asymm)(M_asymm_x, x, PLUS);

		M_asymm_x -= AM_symm_x;
		Double norm_cb0 = sqrt(norm2(M_asymm_x,rb[0]));
		Double norm_cb1 = sqrt(norm2(M_asymm_x,rb[1]));

		QDPIO::cout << "CB=0: || A_oo M_symm_x - M_asymm_x ||= " << norm_cb0 <<std::endl;
		QDPIO::cout << "CB=1: || A_oo M_symm_x - M_asymm_x ||= " << norm_cb1
				<< "     || A_oo M_symm_x - M_asymm_x ||/ ||x||=" << norm_cb1/sqrt(norm2(x,rb[1]))
				<< std::endl;

		ASSERT_LT( toDouble(norm_cb0), 1.0e-14);
		ASSERT_LT( toDouble(norm_cb1), 1.0e-13);
	}
	// Now check herm conj
	{
		QDPIO::cout << "Daggered Op: " << std::endl;
		(*M_symm).unprecOddOddLinOp(t1,x,MINUS);
		(*M_symm)(AM_symm_x,t1,MINUS);

		(*M_asymm)(M_asymm_x, x, MINUS);
		M_asymm_x -= AM_symm_x;
		Double norm_cb0 = sqrt(norm2(M_asymm_x,rb[0]));
		Double norm_cb1 = sqrt(norm2(M_asymm_x,rb[1]));

		QDPIO::cout << "CB=0: || A_oo M_symm_x - M_asymm_x ||= " << norm_cb0 <<std::endl;
		QDPIO::cout << "CB=1: || A_oo M_symm_x - M_asymm_x ||= " << norm_cb1
				<< "     || A_oo M_symm_x - M_asymm_x ||/ ||x||=" << norm_cb1/sqrt(norm2(x,rb[1]))
				<< std::endl;

		ASSERT_LT(  toDouble(norm_cb0), 1.0e-14);
		ASSERT_LT(  toDouble(norm_cb1), 1.0e-13);
	}

}

// Check both symm and asymm linops have same unprec op
TEST_F(SymmFixture, CheckUnprecOp)
{
	T x;
	T unprec_symm_x;
	T unprec_asymm_x;

	gaussian(x);
	{

		QDPIO::cout << "Regular Op:" << std::endl;

		(*M_symm).unprecLinOp( unprec_symm_x, x, PLUS);
		(*M_asymm).unprecLinOp( unprec_asymm_x, x, PLUS);

		unprec_asymm_x -= unprec_symm_x;
		Double norm_cb0 = sqrt(norm2(unprec_asymm_x,rb[0]));
		Double norm_cb1 = sqrt(norm2(unprec_asymm_x,rb[1]));

		QDPIO::cout << "CB=0: || Munprec_symm_x - Munprec_asymm_x ||= " << norm_cb0 <<std::endl;
		QDPIO::cout << "CB=1: || Munprec_symm_x - Munprec_asymm_x ||= " << norm_cb1 << std::endl;

		ASSERT_LT( toDouble(norm_cb0), 1.0e-13);
		ASSERT_LT( toDouble(norm_cb1), 1.0e-13);

	}

	{
		QDPIO::cout << "Daggered Op:" << std::endl;

		(*M_symm).unprecLinOp( unprec_symm_x, x, MINUS);
		(*M_asymm).unprecLinOp( unprec_asymm_x, x, MINUS);
		unprec_asymm_x -= unprec_symm_x;

		Double norm_cb0 = sqrt(norm2(unprec_asymm_x,rb[0]));
		Double norm_cb1 = sqrt(norm2(unprec_asymm_x,rb[1]));

		QDPIO::cout << "CB=0: || Munprec_symm_x - Munprec_asymm_x ||= " << norm_cb0 <<std::endl;
		QDPIO::cout << "CB=1: || Munprec_symm_x - Munprec_asymm_x ||= " << norm_cb1 << std::endl;

		ASSERT_LT( toDouble(norm_cb0), 1.0e-13);
		ASSERT_LT( toDouble(norm_cb1), 1.0e-13);

	}

}

// Check QProp Functionality.
TEST_P(QPropTest, CheckQprop)
{
	LatticeFermion rhs=zero;
	gaussian(rhs);

	std::istringstream inv_param_xml_stream(GetParam());
	XMLReader xml_in(inv_param_xml_stream);

	GroupXML_t inv_param = readXMLGroup(xml_in, "//InvertParam", "invType");
    Handle<SystemSolver<T>>	qprop_solver = S_symm->qprop(state,inv_param);
	LatticeFermion x = zero;

	(*qprop_solver)(x,rhs);

	// Check residuum
	LatticeFermion Ax=zero;
	(*M_symm).unprecLinOp(Ax,x,PLUS);
	Ax -= rhs;

	Double resid_cb0 = sqrt(norm2(Ax,rb[0]));
	Double resid_cb1 = sqrt(norm2(Ax,rb[1]));
	QDPIO::cout << "Qprop: rsd cb0 = " << resid_cb0 << std::endl;
	QDPIO::cout << "Qprop: rsd cb1 = " << resid_cb1 << std::endl;

	Double resid = sqrt(norm2(Ax));
	Double resid_rel = resid/sqrt(norm2(rhs));
	QDPIO::cout << "QProp Check Back: || r || = " << resid << "  || r ||/||b|| = "
			<< resid_rel << std::endl;

	ASSERT_LT(toDouble(resid_rel), 1.0e-8);

}

INSTANTIATE_TEST_CASE_P(PropSyssolver,
                        QPropTest,
                        ::testing::Values(inv_param_syssolver_bicgstab_xml));

#ifdef BUILD_QUDA
INSTANTIATE_TEST_CASE_P(PropQUDASolver,
				        QPropTest,
						::testing::Values(inv_param_quda_bicgstab_xml,
										  inv_param_quda_multigrid_xml));

#endif

class MdagMInvTestSymm : public SymmFixtureT<::testing::TestWithParam<std::string>>{};
class MdagMInvTestAsymm : public SymmFixtureT<::testing::TestWithParam<std::string>>{};

TEST_P(MdagMInvTestSymm, CheckMdagMInvSymm)
{
	std::istringstream inv_param_xml_stream(GetParam());
	XMLReader xml_in(inv_param_xml_stream);

	GroupXML_t inv_param = readXMLGroup(xml_in, "//InvertParam", "invType");
	Handle<MdagMSystemSolver<T>>	MdagM_solver = S_symm->invMdagM(state,inv_param);

	T b; gaussian(b,rb[1]);
	T x; x[rb[1]] = zero;

	(*MdagM_solver)(x,b);

	T tmp, r;

	(*M_symm)(tmp,x,PLUS);
	(*M_symm)(r,tmp,MINUS);

	r[rb[1]] -= b;
	Double resid = sqrt(norm2(r,rb[1]));
	Double resid_rel = resid/sqrt(norm2(b,rb[1]));
	QDPIO::cout << "MdagM check: || r || = " << resid << "   || r || / || b ||=" << resid_rel << std::endl;
	ASSERT_LT( toDouble(resid_rel),1.0e-8);



}

#ifdef BUILD_QUDA
TEST_P(MdagMInvTestSymm, CheckMdagMInvSymmQUDAPredict)
{
	std::istringstream inv_param_xml_stream(GetParam());
	XMLReader xml_in(inv_param_xml_stream);

	GroupXML_t inv_param = readXMLGroup(xml_in, "//InvertParam", "invType");
	Handle<MdagMSystemSolver<T>>	MdagM_solver = S_symm->invMdagM(state,inv_param);

	T b; gaussian(b,rb[1]);
	T x; x[rb[1]] = zero;
	QUDA4DChronoPredictor chrono(5,DEFAULT);
	(*MdagM_solver)(x,b,chrono);

	T tmp, r;

	(*M_symm)(tmp,x,PLUS);
	(*M_symm)(r,tmp,MINUS);

	r[rb[1]] -= b;
	Double resid = sqrt(norm2(r,rb[1]));
	Double resid_rel = resid/sqrt(norm2(b,rb[1]));
	QDPIO::cout << "MdagM check: || r || = " << resid << "   || r || / || b ||=" << resid_rel << std::endl;
	ASSERT_LT( toDouble(resid_rel),1.0e-8);



}
#endif

TEST_P(MdagMInvTestAsymm, CheckMdagMInvAsymm)
{
	std::istringstream inv_param_xml_stream(GetParam());
	XMLReader xml_in(inv_param_xml_stream);

	GroupXML_t inv_param = readXMLGroup(xml_in, "//InvertParam", "invType");
	Handle<MdagMSystemSolver<T>>	MdagM_solver = S_asymm->invMdagM(state,inv_param);

	T b; gaussian(b,rb[1]);
	T x; x[rb[1]] = zero;

	(*MdagM_solver)(x,b);

	T tmp, r;

	(*M_asymm)(tmp,x,PLUS);
	(*M_asymm)(r,tmp,MINUS);

	r[rb[1]] -= b;
	Double resid = sqrt(norm2(r,rb[1]));
	Double resid_rel = resid/sqrt(norm2(b,rb[1]));
	QDPIO::cout << "MdagM check: || r || = " << resid << "   || r || / || b ||=" << resid_rel << std::endl;
	ASSERT_LT( toDouble(resid_rel),1.0e-8);



}

#ifdef BUILD_QUDA
TEST_P(MdagMInvTestAsymm, CheckMdagMInvAsymmQUDAPredict)
{
	std::istringstream inv_param_xml_stream(GetParam());
	XMLReader xml_in(inv_param_xml_stream);

	GroupXML_t inv_param = readXMLGroup(xml_in, "//InvertParam", "invType");
	Handle<MdagMSystemSolver<T>>	MdagM_solver = S_asymm->invMdagM(state,inv_param);

	T b; gaussian(b,rb[1]);
	T x; x[rb[1]] = zero;

	QUDA4DChronoPredictor chrono(5,DEFAULT);

	(*MdagM_solver)(x,b,chrono);

	T tmp, r;

	(*M_asymm)(tmp,x,PLUS);
	(*M_asymm)(r,tmp,MINUS);

	r[rb[1]] -= b;
	Double resid = sqrt(norm2(r,rb[1]));
	Double resid_rel = resid/sqrt(norm2(b,rb[1]));
	QDPIO::cout << "MdagM check: || r || = " << resid << "   || r || / || b ||=" << resid_rel << std::endl;
	ASSERT_LT( toDouble(resid_rel),1.0e-8);

}
#endif

INSTANTIATE_TEST_CASE_P(MdagMInvSysSolver,
                        MdagMInvTestSymm,
                        ::testing::Values(inv_param_syssolver_bicgstab_xml));

#ifdef BUILD_QUDA
INSTANTIATE_TEST_CASE_P(MdagMInvQUDASolver,
                        MdagMInvTestSymm,
						::testing::Values(inv_param_quda_bicgstab_xml,
																  inv_param_quda_multigrid_xml));
#endif

INSTANTIATE_TEST_CASE_P(MdagMInvSysSolver,
                        MdagMInvTestAsymm,
                        ::testing::Values(inv_param_syssolver_bicgstab_xml));

#ifdef BUILD_QUDA
INSTANTIATE_TEST_CASE_P(MdagMInvQUDASolver,
                        MdagMInvTestAsymm,
						::testing::Values(inv_param_quda_bicgstab_asymm_xml,
																  inv_param_quda_multigrid_asymm_xml));
#endif
