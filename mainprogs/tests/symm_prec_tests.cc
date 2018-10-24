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

	    std::istringstream input_asymm_periodic(fermact_xml_asymm_periodic);
	    std::istringstream input_symm_periodic(fermact_xml_symm_periodic);

	    XMLReader xml_in_asymm_periodic(input_asymm_periodic);
	    XMLReader xml_in_symm_periodic(input_symm_periodic);

	    S_asymm = dynamic_cast<S_asymm_T*>(TheFermionActionFactory::Instance().createObject("CLOVER",
	    											   xml_in_asymm,
	    											   "FermionAction"));

	    S_symm = dynamic_cast<S_symm_T*>(TheFermionActionFactory::Instance().createObject("SEOPREC_CLOVER",
	        											   xml_in_symm,
	        											   "FermionAction"));

	    S_asymm_periodic = dynamic_cast<S_asymm_T*>(TheFermionActionFactory::Instance().createObject("CLOVER",
	   	    											   xml_in_asymm_periodic,
	   	    											   "FermionAction"));

	    S_symm_periodic = dynamic_cast<S_symm_T*>(TheFermionActionFactory::Instance().createObject("SEOPREC_CLOVER",
	   	        											   xml_in_symm_periodic,
	   	        											   "FermionAction"));
	    state = S_asymm->createState(u);

	    M_asymm =dynamic_cast<LinOpAsymm_T *>(S_asymm->linOp(state));
	    M_symm = dynamic_cast<LinOpSymm_T *>(S_symm->linOp(state));
	}

	void TearDown() {}

	Q u;
	Handle<S_symm_T> S_symm;
	Handle<S_asymm_T> S_asymm;
	Handle<S_symm_T> S_symm_periodic;
	Handle<S_asymm_T> S_asymm_periodic;
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

class multiMdagMInvTestSymm : public SymmFixtureT<::testing::TestWithParam<std::string>>{};
class multiMdagMInvTestAsymm : public SymmFixtureT<::testing::TestWithParam<std::string>>{};

TEST_P(multiMdagMInvTestSymm, checkMultiShift)
{
	std::istringstream inv_param_xml_stream(GetParam());
	XMLReader xml_in(inv_param_xml_stream);

	GroupXML_t inv_param = readXMLGroup(xml_in, "//InvertParam", "invType");

	Handle<MdagMMultiSystemSolver<T>> multiMdagM = S_symm->mInvMdagM(state,inv_param);

	T rhs;
	gaussian(rhs,rb[1]);

	int n_shift=3;
	multi1d<Real> shifts(n_shift);
	shifts[0] = 0.0001;
	shifts[1] = 0.01;
	shifts[2] = 0.1;

	// Zero the initial guesses
	multi1d<T> solns(n_shift);
	for(int shift=0; shift < n_shift; ++shift) {
		(solns[shift])[rb[1]]=zero;
	}

	//operator() (multi1d< multi1d<T> >& psi, const multi1d<Real>& shifts, const multi1d<T>& chi)
	(*multiMdagM)(solns,shifts,rhs);

	for(int shift = 0; shift < n_shift; ++shift) {
		T r = zero;
		T tmp = zero;
		// r = M^\dag M solns[shift]
		(*M_symm)(tmp,solns[shift],PLUS);
		(*M_symm)(r, tmp, MINUS);

		// r = M^\dag M solns[shift] + shifts[shift] solns[shift]
		//   = (M^\dag M + shifts[shift]) solns[shift]
		r[rb[1]] += shifts[shift]*solns[shift];

		// -residudum
		r[rb[1]] -= rhs;

		Double resid_rel = sqrt(norm2(r,rb[1])/norm2(rhs,rb[1]));
		QDPIO::cout << "shift="<<shift << " || r || / || b ||=" << resid_rel << std::endl;
		ASSERT_LT( toDouble(resid_rel), 1.0e-8);
	}
}

TEST_P(multiMdagMInvTestAsymm, checkMultiShift)
{
	std::istringstream inv_param_xml_stream(GetParam());
	XMLReader xml_in(inv_param_xml_stream);

	GroupXML_t inv_param = readXMLGroup(xml_in, "//InvertParam", "invType");

	Handle<MdagMMultiSystemSolver<T>> multiMdagM = S_asymm->mInvMdagM(state,inv_param);

	T rhs;
	gaussian(rhs,rb[1]);

	int n_shift=3;
	multi1d<Real> shifts(n_shift);
	shifts[0] = 0.0001;
	shifts[1] = 0.01;
	shifts[2] = 0.1;

	// Zero the initial guesses
	multi1d<T> solns(n_shift);
	for(int shift=0; shift < n_shift; ++shift) {
		(solns[shift])[rb[1]]=zero;
	}

	//operator() (multi1d< multi1d<T> >& psi, const multi1d<Real>& shifts, const multi1d<T>& chi)
	(*multiMdagM)(solns,shifts,rhs);

	for(int shift = 0; shift < n_shift; ++shift) {
		T r = zero;
		T tmp = zero;
		// r = M^\dag M solns[shift]
		(*M_asymm)(tmp,solns[shift],PLUS);
		(*M_asymm)(r, tmp, MINUS);

		// r = M^\dag M solns[shift] + shifts[shift] solns[shift]
		//   = (M^\dag M + shifts[shift]) solns[shift]
		r[rb[1]] += shifts[shift]*solns[shift];

		// -residudum
		r[rb[1]] -= rhs;

		Double resid_rel = sqrt(norm2(r,rb[1])/norm2(rhs,rb[1]));
		QDPIO::cout << "shift="<<shift << " || r || / || b ||=" << resid_rel << std::endl;
		ASSERT_LT( toDouble(resid_rel), 1.0e-8);
	}
}

INSTANTIATE_TEST_CASE_P(MultiShiftSysSolver,
                        multiMdagMInvTestAsymm,
                        ::testing::Values(inv_param_multi_cg_xml));

#ifdef BUILD_QUDA
INSTANTIATE_TEST_CASE_P(MultiShiftQUDASolver,
                        multiMdagMInvTestAsymm,
						::testing::Values(inv_param_multi_cg_quda_asymm_xml));
#endif

INSTANTIATE_TEST_CASE_P(MultiShiftSysSolver,
                        multiMdagMInvTestSymm,
                        ::testing::Values(inv_param_multi_cg_xml));

#ifdef BUILD_QUDA
INSTANTIATE_TEST_CASE_P(MultiShiftQUDASolver,
                        multiMdagMInvTestSymm,
						::testing::Values(inv_param_multi_cg_quda_xml));
#endif


// Forces
TEST_F(SymmFixture, TestDeriv)
{
	// M_symm = M_oo^{-1} M_asymm
	// X^\dagger d[ M_symm ] Y
	//  = X^\dagger d[ M_oo^{-1} ] M_asymm Y
	//   +X^\dagger M_oo^{-1} d[ M_asymm ] Y
    //
	// = -X^\dagger M_oo^{-1} d[ M_oo ] M^{-1}_oo M_asymm Y
	//   +X^\dagger M_oo^{-1} d[ M_asymm   ] Y
	//
	// =  -Z^\dagger d[ M_oo ] W
	//    +Z^\dagger d[ M_asymm ] Y
	//
	// with W = M_{oo}^{-1} M_asymm Y
	// and  Z = M_oo^{-dagger} X

	P ds_symm;
	P ds_tmp;
	P rhs;


	T X = zero;
	T Y = zero;

	gaussian(X,rb[1]);
	gaussian(Y,rb[1]);

	T tmp = zero;
	T W = zero;
	T Z = zero;

	// W = M^{-1}_oo M_asymm Y
	(*M_asymm)(tmp,Y,PLUS);
	(*M_symm).unprecOddOddInvLinOp(W,tmp,PLUS);

	(*M_symm).unprecOddOddInvLinOp(Z,X,MINUS);

	// The derivative of M_symm
	(*M_symm).deriv(ds_symm,X,Y,PLUS);

	// rhs = Z^\dagger d[ M_asymm ] Y
	(*M_asymm).deriv(rhs,Z,Y,PLUS);

	// rhs -= Z^\dagger d[ M_oo ] W
	(*M_symm).derivUnprecOddOddLinOp(ds_tmp,Z,W,PLUS);
	rhs -= ds_tmp;

	for(int mu=0; mu < Nd; ++mu) {
		rhs[mu] -= ds_symm[mu];
		Double norm_rhs = sqrt(norm2(rhs[mu]));
		Double norm_rhs_per_number = norm_rhs/Double(3*3*2*Layout::vol());
		QDPIO::cout << "mu=" << mu << " || rhs - ds_symm || = " << norm_rhs
				<< "  || rhs - ds_symm || / number =" << norm_rhs_per_number << std::endl;

		ASSERT_LT(toDouble(norm_rhs_per_number), 1.0e-18 );
	}
}

TEST_F(SymmFixture, TestDerivDagger)
{
	// M_symm^\dagger = M_asymm^\dagger M_oo^{-\dagger}

	// X^\dagger d[ M_symm ] Y
	//  = X^\dagger d[ M_asymm^\dagger ] M_oo^{-\dagger} Y
	//   +X^\dagger M_asymm^\dagger d[ M_oo^{-\dagger} ] Y
    //
	// = +X^\dagger d[M_asymm^\dagger ] M_oo^{-dagger} Y
	//   -X^\dagger M_asymm^\dagger ] M_oo^{-dagger} d[ M_oo^{\dagger} ] M_oo^{-dagger} Y
	//
	// =  +X^\dagger d[ M_asymm^\dagger ] W
	//    -Z^\dagger d[ M_oo^\dagger  ] W
	//
	// with W = M_{oo}^{-\dagger}  Y
	// and  Z = M_oo^{-1} M_asymm  X

	P ds_symm;
	P ds_tmp;
	P rhs;


	T X = zero;
	T Y = zero;

	gaussian(X,rb[1]);
	gaussian(Y,rb[1]);

	T tmp = zero;
	T W = zero;
	T Z = zero;

	// W = M^{-dagger}_oo Y
	(*M_symm).unprecOddOddInvLinOp(W,Y,MINUS);

	//Z = M_oo^{-1} M_asymm  X
	(*M_asymm)(tmp,X,PLUS);
	(*M_symm).unprecOddOddInvLinOp(Z,tmp,PLUS);

	// The derivative of M_symm^dagger
	(*M_symm).deriv(ds_symm,X,Y,MINUS);

	// rhs = X^\dagger d[ M_asymm ] W
	(*M_asymm).deriv(rhs,X,W,MINUS);

	// rhs -= Z^\dagger d[ M_oo ] W
	(*M_symm).derivUnprecOddOddLinOp(ds_tmp,Z,W,MINUS);
	rhs -= ds_tmp;

	for(int mu=0; mu < Nd; ++mu) {
		rhs[mu] -= ds_symm[mu];
		Double norm_rhs = sqrt(norm2(rhs[mu]));
		Double norm_rhs_per_number = norm_rhs/Double(3*3*2*Layout::vol());
		QDPIO::cout << "mu=" << mu << " || rhs - ds_symm || = " << norm_rhs
				<< "  || rhs - ds_symm || / number =" << norm_rhs_per_number << std::endl;

		ASSERT_LT(toDouble(norm_rhs_per_number), 1.0e-18 );
	}
}

TEST_F(SymmFixture,TestLogDetUnitGauge)
{
	Q u_unit;
	u_unit.resize(Nd);
	for(int mu=0; mu < Nd; ++mu) {
		u_unit[mu] = 1;
	}

	Handle< FermState<T,P,Q> > unit_state = S_symm->createState(u_unit);
	QDPIO::cout << "Unit Gauge Symm Linop creation" << std::endl;
	Handle<LinOpSymm_T> lop = dynamic_cast<LinOpSymm_T*>(S_symm->linOp(unit_state));

	Double S_e = lop->logDetEvenEvenLinOp();
	Double S_o = lop->logDetOddOddLinOp();
	QDPIO::cout << "Unit gauge: log Det EvenEven = " << S_e << std::endl;
	QDPIO::cout << "Unit gauge: log Det OddOdd = " << S_o << std::endl;

	Double absdiff = fabs(S_e-S_o);
	QDPIO::cout << "Unit gauge: | S_e - S_o | = " << absdiff << std::endl;
	ASSERT_LT( toDouble(absdiff), 1.0e-14);
}

TEST_F(SymmFixture,TestLogDetShiftedGauge)
{
	Q u_shifted;
	u_shifted.resize(Nd);
	for(int mu=0; mu < Nd; ++mu) {
		u_shifted[mu] = shift(u[mu],FORWARD,3);
	}

	Handle< FermState<T,P,Q> > shifted_state = S_symm->createState(u_shifted);

	QDPIO::cout << "Shifted Gauge Symm Linop creation" << std::endl;

	Handle<LinOpSymm_T> shifted_lop = dynamic_cast<LinOpSymm_T*>(S_symm->linOp(shifted_state));

	Double S_e_asymm = M_asymm->logDetEvenEvenLinOp();
	Double S_e = M_symm->logDetEvenEvenLinOp();
	Double S_o = M_symm->logDetOddOddLinOp();

	Double S_e_shifted = shifted_lop->logDetEvenEvenLinOp();
	Double S_o_shifted = shifted_lop->logDetOddOddLinOp();

	QDPIO::cout << "Asymm op log Det EvenEven = " << S_e_asymm << std::endl;
	QDPIO::cout << "Random gauge: log Det EvenEven = " << S_e << std::endl;
	QDPIO::cout << "Random gauge: log Det OddOdd = " << S_o << std::endl;
	QDPIO::cout << "Shifted gauge: log Det EvenEven = " << S_e_shifted << std::endl;
	QDPIO::cout << "Shifted gauge: log Det OddOdd = " << S_o_shifted << std::endl;

	Double diff_e_symm_asymm = fabs(S_e_asymm - S_e);
	Double diff_eo = fabs(S_e - S_o_shifted)/Double(rb[0].numSiteTable());
	Double diff_oe = fabs(S_o - S_e_shifted)/Double(rb[1].numSiteTable());
	QDPIO::cout << "| logDet_e_asymm - logdet_e_symm | = " << diff_e_symm_asymm << std::endl;
	QDPIO::cout << "| logDet_e - logDet_o_shifted |/site = " << diff_eo << std::endl;
	QDPIO::cout << "| logDet_o - logDet_e_shifted |/site = " << diff_oe << std::endl;

	ASSERT_LT( toDouble(diff_e_symm_asymm), 1.0e-14);
	ASSERT_LT( toDouble(diff_eo), 1.0e-13);
	ASSERT_LT( toDouble(diff_oe), 1.0e-13);

}

class TrLogForceFixture : public SymmFixtureT<::testing::TestWithParam<enum PlusMinus>>{};
TEST_P(TrLogForceFixture,TestShiftedGaugeTrLnForce)
{
	Q u_shifted;
	u_shifted.resize(Nd);
	for(int mu=0; mu < Nd; ++mu) {
		u_shifted[mu] = shift(u[mu],FORWARD,3);
	}

	Handle< FermState<T,P,Q> > periodic_state = S_symm_periodic->createState(u);
	Handle< FermState<T,P,Q> > periodic_shifted_state = S_symm_periodic->createState(u_shifted);

	QDPIO::cout << "Shifted Gauge Symm Linop creation" << std::endl;

	Handle<LinOpAsymm_T> asymm_periodic_op = dynamic_cast<LinOpAsymm_T*>(S_asymm_periodic->linOp(periodic_state));
	Handle<LinOpSymm_T> symm_periodic_op = dynamic_cast<LinOpSymm_T*>(S_symm_periodic->linOp(periodic_state));
	Handle<LinOpSymm_T> shifted_periodic_op = dynamic_cast<LinOpSymm_T*>(S_symm_periodic->linOp(periodic_shifted_state));

	//! Get the force from the EvenEven Trace Log
	//   virtual void derivLogDetEvenEvenLinOp(P& ds_u, enum PlusMinus isign) const = 0;
	P ds_asymm;
	asymm_periodic_op->derivLogDetEvenEvenLinOp(ds_asymm,GetParam());

	P ds_symm;
	symm_periodic_op->derivLogDetEvenEvenLinOp(ds_symm,GetParam());

	P ds_symm_shifted;
	shifted_periodic_op->derivLogDetOddOddLinOp(ds_symm_shifted, GetParam());

	P ds_unshifted; ds_unshifted.resize(Nd);
	P ds_wrong_shifted; ds_wrong_shifted.resize(Nd);
	for(int mu=0; mu < Nd; ++mu ) {
		ds_unshifted[mu] = shift(ds_symm_shifted[mu], BACKWARD,3);
		ds_wrong_shifted[mu] = shift(ds_symm_shifted[mu], FORWARD, 3);
	}
	for(int mu=0; mu < Nd; ++mu ) {
		Double mu_contrib_asymm = norm2(ds_asymm[mu]);

		Double mu_contrib_symm = norm2(ds_symm[mu]);

		Double mu_contrib_shifted = norm2(ds_symm_shifted[mu]);

		QDPIO::cout << "mu=" << mu << "  asymm force norm="<< mu_contrib_asymm << std::endl;
		QDPIO::cout << "mu=" << mu << "  symm force norm ="<< mu_contrib_symm << std::endl;
		QDPIO::cout << "mu=" << mu << "  shifted force norm="<< mu_contrib_shifted << std::endl;


		Double diff_symm_asymm_ee = fabs(mu_contrib_asymm - mu_contrib_symm );
		Double diff_symm_ee_oo = fabs(mu_contrib_symm - mu_contrib_shifted );

		QDPIO::cout << "mu=" << mu << " | F_asymm_ee - F_symm_ee | = " << diff_symm_asymm_ee << std::endl;
		QDPIO::cout << "mu=" << mu << " | F_ee - F_shifted_oo | = " << diff_symm_ee_oo << std::endl;

		ASSERT_LT( toDouble(diff_symm_asymm_ee), 1.0e-15);
		ASSERT_LT( toDouble(diff_symm_ee_oo), 1.0e-15);

		ds_unshifted[mu] -= ds_symm[mu];
		ds_wrong_shifted[mu] -= ds_symm[mu];
		Double diffnorm_unshifted = sqrt(norm2(ds_unshifted[mu]));
		Double diffnorm_per_site = diffnorm_unshifted/Layout::vol();
		QDPIO::cout << "mu=" << mu << " || F_mu - shift(F_shifted, BACKWARD, mu) || =" << diffnorm_unshifted << std::endl;
		QDPIO::cout << "mu=" << mu << " || F_mu - shift(F_shifted, BACKWARD, mu) ||/site =" << diffnorm_per_site << std::endl;

		ASSERT_LT( toDouble(diffnorm_per_site), 1.0e-14 );

		Double diffnorm_wrong_shifted = sqrt(norm2(ds_wrong_shifted[mu]));
		QDPIO::cout << "mu=" << mu << "  DELIBERATELY WRONG: || F_mu - shift(F_shifted, FORWARD, mu) || =" << diffnorm_wrong_shifted << std::endl;
		ASSERT_GT( toDouble(diffnorm_wrong_shifted), 0.5);
		QDPIO::cout << std::endl;
	}

}

INSTANTIATE_TEST_CASE_P(TrLogForces,
		TrLogForceFixture,
		::testing::Values(PLUS,MINUS));
