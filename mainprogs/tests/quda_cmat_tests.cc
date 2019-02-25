#include <iostream>

#include "chromabase.h"

#include "handle.h"

#include "seoprec_linop.h"
#include "eoprec_linop.h"
#include "eoprec_wilstype_fermact_w.h"
#include "seoprec_wilstype_fermact_w.h"
#include "actions/ferm/linop/eoprec_clover_linop_w.h"
#include "actions/ferm/linop/seoprec_clover_linop_w.h"
#include "actions/ferm/linop/shifted_linop_w.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/linop/unprec_clover_plus_igmuA2_linop_w.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "util/gauge/reunit.h"
#include "gtest/gtest.h"

#include "actions/ferm/invert/quda_solvers/syssolver_linop_clover_quda_w.h"
#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_factory.h"
// #include "actions/ferm/linop/unprec_clover_plus_igmuA2_linop_w.h"

#include "io/xml_group_reader.h"
#ifdef BUILD_QUDA
#include "quda.h"
#endif

#include "symm_prec_xml.h"
#include "quda_cmat_xml.h"
#include "actions/ferm/linop/shifted_linop_w.h"
#include "actions/ferm/linop/unprec_clover_plus_igmuA2_linop_w.h"


using namespace Chroma;
using namespace QDP;
using namespace SymmPrecTesting;

template<typename TestType>
class QudaFixtureT : public TestType {
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


		std::istringstream inv_param_xml_stream(inv_param_quda_bicgstab_xml);
		std::istringstream inv_param_twisted_gcr_stream(inv_param_quda_gcr_twisted_xml);
		XMLReader xml_in(inv_param_xml_stream);
		GroupXML_t inv_param = readXMLGroup(xml_in, "//InvertParam", "invType");


	    linop_solver = S_symm->invLinOp(state,inv_param);

	}


	void TearDown() {}

	LinOpSysSolverQUDAClover& getSolver()
	{
		return dynamic_cast<LinOpSysSolverQUDAClover&>(*linop_solver);
	}

	Q u;
	Handle<S_symm_T> S_symm;
	Handle<S_asymm_T> S_asymm;
	Handle<FermState<T,P,Q> > state;

	Handle<LinOpAsymm_T> M_asymm;
	Handle<LinOpSymm_T> M_symm;
	Handle<SystemSolver<T>>	linop_solver;
};

class QudaFixture : public QudaFixtureT<::testing::Test> {};

TEST_F(QudaFixture, TestCloverMat)
{
	auto the_quda_solver = getSolver();
	auto quda_inv_param = the_quda_solver.getQudaInvertParam();

	for(int dagger = 0; dagger < 2; ++dagger) {

		enum PlusMinus isign = ( dagger == 0 ) ? PLUS : MINUS;
		quda_inv_param.dagger = (dagger == 0 ) ? QUDA_DAG_NO : QUDA_DAG_YES;
		std::string op_str = (dagger == 0 ) ? "Op :" : "Dag :" ;
		T src=zero;
		T res=zero;
		T res_quda = zero;

		// Apply symmetric clover operator
		gaussian(src,rb[1]);
		(*M_symm)(res,src,isign);
		auto src_ptr = (void *)&(src.elem(rb[1].start()).elem(0).elem(0).real());
		auto res_quda_ptr = (void *)&(res_quda.elem(rb[1].start()).elem(0).elem(0).real());

		//quda_inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
		// Now I want to apply the QUDA operator.
		MatQuda(res_quda_ptr, src_ptr, &quda_inv_param);

		T diff = zero;
		diff[rb[1]] = res_quda - res;
		Double norm_diff = norm2(diff,rb[1]);
		QDPIO::cout << op_str << " || QUDA - Chroma || = " << sqrt(norm_diff) << std::endl;
		Double sites = Double(Layout::vol()/2);
		Double norm_diff_per_site = sqrt(norm_diff/sites);
		QDPIO::cout << op_str << " || QUDA - Chroma ||/site = " << norm_diff_per_site << std::endl;

		ASSERT_LT( toDouble(norm_diff_per_site), 1.0e-14);
	}
}

TEST_F(QudaFixture, TestUnprecCloverMat)
{
	auto the_quda_solver = getSolver();
	auto quda_inv_param = the_quda_solver.getQudaInvertParam();
	quda_inv_param.solution_type = QUDA_MAT_SOLUTION;
	quda_inv_param.matpc_type = QUDA_MATPC_INVALID;

	for(int dagger = 0; dagger < 2; ++dagger) {

		enum PlusMinus isign = ( dagger == 0 ) ? PLUS : MINUS;
		quda_inv_param.dagger = (dagger == 0 ) ? QUDA_DAG_NO : QUDA_DAG_YES;
		std::string op_str = (dagger == 0 ) ? "Op :" : "Dag :" ;
		T src=zero;
		T res=zero;
		T res_quda = zero;

		// Apply symmetric clover operator
		gaussian(src);
		(*M_symm).unprecLinOp(res,src,isign);
		auto src_ptr = (void *)&(src.elem(rb[0].start()).elem(0).elem(0).real());
		auto res_quda_ptr = (void *)&(res_quda.elem(rb[0].start()).elem(0).elem(0).real());


		//quda_inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
		// Now I want to apply the QUDA operator.
		MatQuda(res_quda_ptr, src_ptr, &quda_inv_param);

		T diff = zero;
		diff = res_quda - res;
		Double norm_diff = norm2(diff);
		QDPIO::cout << op_str << " || QUDA - Chroma || = " << sqrt(norm_diff) << std::endl;
		Double sites = Double(Layout::vol());
		Double norm_diff_per_site = sqrt(norm_diff/sites);
		QDPIO::cout << op_str << " || QUDA - Chroma ||/site = " << norm_diff_per_site << std::endl;

		ASSERT_LT( toDouble(norm_diff_per_site), 1.0e-14);
	}
}

TEST_F(QudaFixture, TestUnprecTwistedMatDispatchZeroTwist)
{
	auto the_quda_solver = getSolver();
	auto quda_inv_param = the_quda_solver.getQudaInvertParam();
	quda_inv_param.solution_type = QUDA_MAT_SOLUTION;
	quda_inv_param.matpc_type = QUDA_MATPC_INVALID;
	quda_inv_param.dslash_type = QUDA_CLOVER_HASENBUSCH_TWIST_DSLASH;
	quda_inv_param.m5 = 0.0;
	QDPIO::cout << "Twist is " << quda_inv_param.m5 << std::endl;
	for(int dagger = 0; dagger < 2; ++dagger) {

		enum PlusMinus isign = ( dagger == 0 ) ? PLUS : MINUS;
		quda_inv_param.dagger = (dagger == 0 ) ? QUDA_DAG_NO : QUDA_DAG_YES;
		std::string op_str = (dagger == 0 ) ? "Op :" : "Dag :" ;
		T src=zero;
		T res=zero;
		T res_quda = zero;

		// Apply symmetric clover operator
		gaussian(src);
		(*M_symm).unprecLinOp(res,src,isign);
		auto src_ptr = (void *)&(src.elem(rb[0].start()).elem(0).elem(0).real());
		auto res_quda_ptr = (void *)&(res_quda.elem(rb[0].start()).elem(0).elem(0).real());


		//quda_inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
		// Now I want to apply the QUDA operator.
		MatQuda(res_quda_ptr, src_ptr, &quda_inv_param);

		T diff = zero;
		diff = res_quda - res;
		Double norm_diff = norm2(diff);
		QDPIO::cout << op_str << " || QUDA - Chroma || = " << sqrt(norm_diff) << std::endl;
		Double sites = Double(Layout::vol());
		Double norm_diff_per_site = sqrt(norm_diff/sites);
		QDPIO::cout << op_str << " || QUDA - Chroma ||/site = " << norm_diff_per_site << std::endl;

		ASSERT_LT( toDouble(norm_diff_per_site), 1.0e-14);
	}
}

TEST_F(QudaFixture, TestUnprecTwistedMatDispatchNonZeroTwist)
{
	auto the_quda_solver = getSolver();
	auto quda_inv_param = the_quda_solver.getQudaInvertParam();
	quda_inv_param.solution_type = QUDA_MAT_SOLUTION;
	quda_inv_param.matpc_type = QUDA_MATPC_INVALID;
	quda_inv_param.dslash_type = QUDA_CLOVER_HASENBUSCH_TWIST_DSLASH;


	CloverFermActParams p;
	p.Mass = Real(0.1);
	p.clovCoeffR=Real(1.0);
	p.clovCoeffT=Real(1.0);
	p.u0 = Real(1);
	p.anisoParam.anisoP=false;
	p.anisoParam.t_dir=3;
	p.anisoParam.xi_0 =Real(1);
	p.anisoParam.nu=Real(1);
	p.twisted_m_usedP = true;
	p.twisted_m = Real(0.2345);

	// There is a normalization difference with a sign of -4
	quda_inv_param.m5 =-4*toDouble(p.twisted_m);


	UnprecCloverPlusIG5MuA2Linop M_u_tw(state,p);

	QDPIO::cout << "Twist is " << quda_inv_param.m5 << std::endl;
	for(int dagger = 0; dagger < 2; ++dagger) {

		enum PlusMinus isign = ( dagger == 0 ) ? PLUS : MINUS;
		quda_inv_param.dagger = (dagger == 0 ) ? QUDA_DAG_NO : QUDA_DAG_YES;
		std::string op_str = (dagger == 0 ) ? "Op :" : "Dag :" ;
		T src=zero;
		T res=zero;
		T res_quda = zero;

		// Apply symmetric clover operator
		gaussian(src);
		M_u_tw(res,src,isign);
		auto src_ptr = (void *)&(src.elem(rb[0].start()).elem(0).elem(0).real());
		auto res_quda_ptr = (void *)&(res_quda.elem(rb[0].start()).elem(0).elem(0).real());


		//quda_inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
		// Now I want to apply the QUDA operator.
		MatQuda(res_quda_ptr, src_ptr, &quda_inv_param);

		T diff = zero;
		diff = res_quda - res;
		Double norm_diff = norm2(diff);
		QDPIO::cout << op_str << " || QUDA - Chroma || = " << sqrt(norm_diff) << std::endl;
		Double sites = Double(Layout::vol());
		Double norm_diff_per_site = sqrt(norm_diff/sites);
		QDPIO::cout << op_str << " || QUDA - Chroma ||/site = " << norm_diff_per_site << std::endl;

		ASSERT_LT( toDouble(norm_diff_per_site), 1.0e-14);
	}
}



TEST_F(QudaFixture, TestPrecTwistedAsymmMatZeroTwist)
{
	auto the_quda_solver = getSolver();
	auto quda_inv_param = the_quda_solver.getQudaInvertParam();
	quda_inv_param.solution_type = QUDA_MATPC_SOLUTION;
	quda_inv_param.matpc_type = QUDA_MATPC_ODD_ODD_ASYMMETRIC;
	quda_inv_param.dslash_type = QUDA_CLOVER_HASENBUSCH_TWIST_DSLASH;

	// There is a normalization difference with a sign of -4
	quda_inv_param.m5 = 0;

	QDPIO::cout << "Twist is " << quda_inv_param.m5 << std::endl;
	for(int dagger = 0; dagger < 2; ++dagger) {

		enum PlusMinus isign = ( dagger == 0 ) ? PLUS : MINUS;
		quda_inv_param.dagger = (dagger == 0 ) ? QUDA_DAG_NO : QUDA_DAG_YES;
		std::string op_str = (dagger == 0 ) ? "Op :" : "Dag :" ;
		T src=zero;
		T res=zero;
		T res_quda = zero;

		// Apply symmetric clover operator
		gaussian(src, rb[1]);
		(*M_asymm)(res,src,isign);
		auto src_ptr = (void *)&(src.elem(rb[1].start()).elem(0).elem(0).real());
		auto res_quda_ptr = (void *)&(res_quda.elem(rb[1].start()).elem(0).elem(0).real());


		//quda_inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
		// Now I want to apply the QUDA operator.
		MatQuda(res_quda_ptr, src_ptr, &quda_inv_param);

		T diff = zero;
		diff[rb[1]] = res_quda - res;
		Double norm_diff = norm2(diff,rb[1]);
		QDPIO::cout << op_str << " || QUDA - Chroma || = " << sqrt(norm_diff) << std::endl;
		Double sites = Double(Layout::vol())/Double(2);
		Double norm_diff_per_site = sqrt(norm_diff/sites);
		QDPIO::cout << op_str << " || QUDA - Chroma ||/site = " << norm_diff_per_site << std::endl;

		ASSERT_LT( toDouble(norm_diff_per_site), 1.0e-14);
	}
}

TEST_F(QudaFixture, TestPrecTwistedAsymmMatNonZeroTwist)
{
	auto the_quda_solver = getSolver();
	auto quda_inv_param = the_quda_solver.getQudaInvertParam();
	quda_inv_param.solution_type = QUDA_MATPC_SOLUTION;
	quda_inv_param.matpc_type = QUDA_MATPC_ODD_ODD_ASYMMETRIC;
	quda_inv_param.dslash_type = QUDA_CLOVER_HASENBUSCH_TWIST_DSLASH;

	Real twist = Real(0.2345);

	// There is a normalization difference with a sign of -4
	quda_inv_param.m5 = -toDouble(twist);

	QDPIO::cout << "Twist is " << quda_inv_param.m5 << std::endl;

	TwistedShiftedLinOp<T, P, Q, EvenOddPrecLinearOperator> M_shifted(*M_asymm, twist);

	for(int dagger = 0; dagger < 2; ++dagger) {

		enum PlusMinus isign = ( dagger == 0 ) ? PLUS : MINUS;
		quda_inv_param.dagger = (dagger == 0 ) ? QUDA_DAG_NO : QUDA_DAG_YES;
		std::string op_str = (dagger == 0 ) ? "Op :" : "Dag :" ;
		T src=zero;
		T res=zero;
		T res_quda = zero;

		// Apply symmetric clover operator
		gaussian(src, rb[1]);
		M_shifted(res,src,isign);
		auto src_ptr = (void *)&(src.elem(rb[1].start()).elem(0).elem(0).real());
		auto res_quda_ptr = (void *)&(res_quda.elem(rb[1].start()).elem(0).elem(0).real());


		//quda_inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
		// Now I want to apply the QUDA operator.
		MatQuda(res_quda_ptr, src_ptr, &quda_inv_param);

		T diff = zero;
		diff[rb[1]] = res_quda - res;
		Double norm_diff = norm2(diff,rb[1]);
		QDPIO::cout << op_str << " || QUDA - Chroma || = " << sqrt(norm_diff) << std::endl;
		Double sites = Double(Layout::vol())/Double(2);
		Double norm_diff_per_site = sqrt(norm_diff/sites);
		QDPIO::cout << op_str << " || QUDA - Chroma ||/site = " << norm_diff_per_site << std::endl;

		ASSERT_LT( toDouble(norm_diff_per_site), 1.0e-14);
	}
}



TEST_F(QudaFixture, TestPrecTwistedMatZeroTwist)
{
	auto the_quda_solver = getSolver();
	auto quda_inv_param = the_quda_solver.getQudaInvertParam();
	quda_inv_param.solution_type = QUDA_MATPC_SOLUTION;
	quda_inv_param.matpc_type = QUDA_MATPC_ODD_ODD;
	quda_inv_param.dslash_type = QUDA_CLOVER_HASENBUSCH_TWIST_DSLASH;

	// There is a normalization difference with a sign of -4
	quda_inv_param.m5 = 0;

	QDPIO::cout << "Twist is " << quda_inv_param.m5 << std::endl;
	for(int dagger = 0; dagger < 2; ++dagger) {

		enum PlusMinus isign = ( dagger == 0 ) ? PLUS : MINUS;
		quda_inv_param.dagger = (dagger == 0 ) ? QUDA_DAG_NO : QUDA_DAG_YES;
		std::string op_str = (dagger == 0 ) ? "Op :" : "Dag :" ;
		T src=zero;
		T res=zero;
		T res_quda = zero;

		// Apply symmetric clover operator
		gaussian(src, rb[1]);
		(*M_symm)(res,src,isign);
		auto src_ptr = (void *)&(src.elem(rb[1].start()).elem(0).elem(0).real());
		auto res_quda_ptr = (void *)&(res_quda.elem(rb[1].start()).elem(0).elem(0).real());


		//quda_inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
		// Now I want to apply the QUDA operator.
		MatQuda(res_quda_ptr, src_ptr, &quda_inv_param);

		T diff = zero;
		diff[rb[1]] = res_quda - res;
		Double norm_diff = norm2(diff,rb[1]);
		QDPIO::cout << op_str << " || QUDA - Chroma || = " << sqrt(norm_diff) << std::endl;
		Double sites = Double(Layout::vol())/Double(2);
		Double norm_diff_per_site = sqrt(norm_diff/sites);
		QDPIO::cout << op_str << " || QUDA - Chroma ||/site = " << norm_diff_per_site << std::endl;

		ASSERT_LT( toDouble(norm_diff_per_site), 1.0e-14);
	}
}


TEST_F(QudaFixture, TestPrecTwistedMatNonZeroTwist)
{
	auto the_quda_solver = getSolver();
	auto quda_inv_param = the_quda_solver.getQudaInvertParam();
	quda_inv_param.solution_type = QUDA_MATPC_SOLUTION;
	quda_inv_param.matpc_type = QUDA_MATPC_ODD_ODD;
	quda_inv_param.dslash_type = QUDA_CLOVER_HASENBUSCH_TWIST_DSLASH;

	Real twist = Real(0.2345);
	// There is a normalization difference with a sign of -4
	quda_inv_param.m5 = -2.0*toDouble(twist);

	QDPIO::cout << "Twist is " << quda_inv_param.m5 << std::endl;

	TwistedShiftedLinOp<T, P, Q, SymEvenOddPrecLogDetLinearOperator> M_shifted(*M_symm, twist);

	// Let us not do dagger initially, get the op working.
	for(int dagger = 0; dagger < 2; ++dagger) {

		enum PlusMinus isign = ( dagger == 0 ) ? PLUS : MINUS;
		quda_inv_param.dagger = (dagger == 0 ) ? QUDA_DAG_NO : QUDA_DAG_YES;
		std::string op_str = (dagger == 0 ) ? "Op :" : "Dag :" ;
		T src=zero;
		T res=zero;
		T res_quda = zero;

		// Apply symmetric clover operator
		gaussian(src, rb[1]);
		M_shifted(res,src,isign);
		auto src_ptr = (void *)&(src.elem(rb[1].start()).elem(0).elem(0).real());
		auto res_quda_ptr = (void *)&(res_quda.elem(rb[1].start()).elem(0).elem(0).real());


		//quda_inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
		// Now I want to apply the QUDA operator.
		MatQuda(res_quda_ptr, src_ptr, &quda_inv_param);

		T diff = zero;
		diff[rb[1]] = res_quda - res;
		Double norm_diff = norm2(diff,rb[1]);
		QDPIO::cout << op_str << " || QUDA - Chroma || = " << sqrt(norm_diff) << std::endl;
		Double sites = Double(Layout::vol())/Double(2);
		Double norm_diff_per_site = sqrt(norm_diff/sites);
		QDPIO::cout << op_str << " || QUDA - Chroma ||/site = " << norm_diff_per_site << std::endl;

		ASSERT_LT( toDouble(norm_diff_per_site), 1.0e-14);
	}
}

TEST_F(QudaFixture, TestSymmPrecTwistedNonZeroTwistMdagM)
{
	auto the_quda_solver = getSolver();
	auto quda_inv_param = the_quda_solver.getQudaInvertParam();
	quda_inv_param.solution_type = QUDA_MATPC_SOLUTION;
	quda_inv_param.matpc_type = QUDA_MATPC_ODD_ODD;
	quda_inv_param.dslash_type = QUDA_CLOVER_HASENBUSCH_TWIST_DSLASH;

	Real twist = Real(0.2345);
	// There is a normalization difference with a sign of -4
	quda_inv_param.m5 = -2.0*toDouble(twist);

	QDPIO::cout << "Twist is " << quda_inv_param.m5 << std::endl;

	TwistedShiftedLinOp<T, P, Q, SymEvenOddPrecLogDetLinearOperator> M_shifted(*M_symm, twist);

	quda_inv_param.dagger = QUDA_DAG_NO;

	T src=zero;
	T res=zero;
	T res2=zero;
	T res_quda = zero;

	// Apply symmetric clover operator
		gaussian(src, rb[1]);
		M_shifted(res2,src,PLUS);
		M_shifted(res,res2,MINUS);

		auto src_ptr = (void *)&(src.elem(rb[1].start()).elem(0).elem(0).real());
		auto res_quda_ptr = (void *)&(res_quda.elem(rb[1].start()).elem(0).elem(0).real());


		//quda_inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
		// Now I want to apply the QUDA operator.
		MatDagMatQuda(res_quda_ptr, src_ptr, &quda_inv_param);

		T diff = zero;
		diff[rb[1]] = res_quda - res;
		Double norm_diff = norm2(diff,rb[1]);
		QDPIO::cout << "MdagM: || QUDA - Chroma || = " << sqrt(norm_diff) << std::endl;
		Double sites = Double(Layout::vol())/Double(2);
		Double norm_diff_per_site = sqrt(norm_diff/sites);
		QDPIO::cout << "MdagM || QUDA - Chroma ||/site = " << norm_diff_per_site << std::endl;

		ASSERT_LT( toDouble(norm_diff_per_site), 1.0e-14);

}

TEST_F(QudaFixture, TestAsymmPrecTwistedNonZeroTwistMdagM)
{
	auto the_quda_solver = getSolver();
	auto quda_inv_param = the_quda_solver.getQudaInvertParam();
	quda_inv_param.solution_type = QUDA_MATPC_SOLUTION;
	quda_inv_param.matpc_type = QUDA_MATPC_ODD_ODD_ASYMMETRIC;
	quda_inv_param.dslash_type = QUDA_CLOVER_HASENBUSCH_TWIST_DSLASH;

	Real twist = Real(0.2345);
	// There is a normalization difference with a sign of -4
	quda_inv_param.m5 = -toDouble(twist);

	QDPIO::cout << "Twist is " << quda_inv_param.m5 << std::endl;

	TwistedShiftedLinOp<T, P, Q, EvenOddPrecLinearOperator> M_shifted(*M_asymm, twist);

	quda_inv_param.dagger = QUDA_DAG_NO;

	T src=zero;
	T res=zero;
	T res2=zero;
	T res_quda = zero;

	// Apply symmetric clover operator
		gaussian(src, rb[1]);
		M_shifted(res2,src,PLUS);
		M_shifted(res,res2,MINUS);

		auto src_ptr = (void *)&(src.elem(rb[1].start()).elem(0).elem(0).real());
		auto res_quda_ptr = (void *)&(res_quda.elem(rb[1].start()).elem(0).elem(0).real());


		//quda_inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
		// Now I want to apply the QUDA operator.
		MatDagMatQuda(res_quda_ptr, src_ptr, &quda_inv_param);

		T diff = zero;
		diff[rb[1]] = res_quda - res;
		Double norm_diff = norm2(diff,rb[1]);
		QDPIO::cout << "MdagM: || QUDA - Chroma || = " << sqrt(norm_diff) << std::endl;
		Double sites = Double(Layout::vol())/Double(2);
		Double norm_diff_per_site = sqrt(norm_diff/sites);
		QDPIO::cout << "MdagM || QUDA - Chroma ||/site = " << norm_diff_per_site << std::endl;

		ASSERT_LT( toDouble(norm_diff_per_site), 5.0e-14);

}

TEST_F(QudaFixture, TestShiftedSymmGCRSolver )
{
	Real twist = Real(0.05);
	std::istringstream inv_param_twisted_gcr_stream(inv_param_quda_gcr_twisted_xml);
	XMLReader xml_in_gcr_twisted(inv_param_twisted_gcr_stream);
	GroupXML_t inv_param_gcr_twisted = readXMLGroup(xml_in_gcr_twisted, "//InvertParam", "invType");
	Handle<SystemSolver<T>>	twisted_gcr_solver=S_symm->invLinOp(state,inv_param_gcr_twisted);
	TwistedShiftedLinOp<T, P, Q, SymEvenOddPrecLogDetLinearOperator> M_shifted(*M_symm, twist);

	auto downcast_solver=dynamic_cast<LinOpSysSolverQUDAClover&>(*twisted_gcr_solver);
	auto quda_inv_param = downcast_solver.getQudaInvertParam();
	// Let us not do dagger initially, get the op working.
	for(int dagger = 0; dagger < 2; ++dagger) {

		enum PlusMinus isign = ( dagger == 0 ) ? PLUS : MINUS;
		quda_inv_param.dagger = (dagger == 0 ) ? QUDA_DAG_NO : QUDA_DAG_YES;
		std::string op_str = (dagger == 0 ) ? "Op :" : "Dag :" ;
		T src=zero;
		T res=zero;
		T res_quda = zero;

		// Apply symmetric clover operator
		gaussian(src, rb[1]);
		M_shifted(res,src,isign);
		auto src_ptr = (void *)&(src.elem(rb[1].start()).elem(0).elem(0).real());
		auto res_quda_ptr = (void *)&(res_quda.elem(rb[1].start()).elem(0).elem(0).real());


		//quda_inv_param.matpc_type = QUDA_MATPC_EVEN_EVEN;
		// Now I want to apply the QUDA operator.
		MatQuda(res_quda_ptr, src_ptr, &quda_inv_param);

		T diff = zero;
		diff[rb[1]] = res_quda - res;
		Double norm_diff = norm2(diff,rb[1]);
		QDPIO::cout << op_str << " || QUDA - Chroma || = " << sqrt(norm_diff) << std::endl;
		Double sites = Double(Layout::vol())/Double(2);
		Double norm_diff_per_site = sqrt(norm_diff/sites);
		QDPIO::cout << op_str << " || QUDA - Chroma ||/site = " << norm_diff_per_site << std::endl;

		ASSERT_LT( toDouble(norm_diff_per_site), 1.0e-14);
	}


}

TEST_F(QudaFixture, TestShiftedSymmGCRSolverSolve )
{
	// Create Twisted Op
	Real twist = Real(0.05);
	Handle<LinearOperator<T>> M_shifted(
	  new TwistedShiftedLinOp<T, P, Q, SymEvenOddPrecLogDetLinearOperator>(*M_symm, twist));

	// Create a Group XML from the solver input XML
	std::istringstream inv_param_twisted_gcr_stream(inv_param_quda_gcr_twisted_xml);
	XMLReader xml_in_gcr_twisted(inv_param_twisted_gcr_stream);
	GroupXML_t inv_param_gcr_twisted = readXMLGroup(xml_in_gcr_twisted, "//InvertParam", "invType");

	// Create a soplver from the factory using the twisted op.
	std::istringstream  xml(inv_param_gcr_twisted.xml);
	XMLReader  paramtop(xml);

	Handle<SystemSolver<T>>	twisted_gcr_solver=
			TheLinOpFermSystemSolverFactory::Instance().createObject(inv_param_gcr_twisted.id,
					paramtop,
					inv_param_gcr_twisted.path,
					state,
					M_shifted);



	T src = zero;
	T res = zero;
			// Apply symmetric clover operator
	gaussian(src, rb[1]);

	(*twisted_gcr_solver)(res, src);

	T M_res=zero;
	(*M_shifted)(M_res,res,PLUS);
	M_res[rb[1]] -= src;
	Double norm_diff = norm2(M_res,rb[1]);
	QDPIO::cout << " || b - Mx || = " << sqrt(norm_diff) << std::endl;
	Double rel_norm_diff = sqrt(norm_diff/norm2(src,rb[1]));
	QDPIO::cout <<  " || b - Mx ||/ b= " << rel_norm_diff << std::endl;

	ASSERT_LT( toDouble(rel_norm_diff), 1.0e-8);


}


TEST_F(QudaFixture, TestShiftedAsymmGCRSolverSolve )
{
	// Create Twisted Op
	Real twist = Real(0.05);
	Handle<LinearOperator<T>> M_shifted(
			new TwistedShiftedLinOp<T, P, Q, EvenOddPrecLinearOperator>(*M_asymm, twist));

	// Create a Group XML from the solver input XML
	std::istringstream inv_param_twisted_gcr_stream(inv_param_quda_gcr_twisted_asymm_xml);
	XMLReader xml_in_gcr_twisted(inv_param_twisted_gcr_stream);
	GroupXML_t inv_param_gcr_twisted = readXMLGroup(xml_in_gcr_twisted, "//InvertParam", "invType");

	// Create a soplver from the factory using the twisted op.
	std::istringstream  xml(inv_param_gcr_twisted.xml);
	XMLReader  paramtop(xml);

	Handle<SystemSolver<T>>	twisted_gcr_solver=
			TheLinOpFermSystemSolverFactory::Instance().createObject(inv_param_gcr_twisted.id,
					paramtop,
					inv_param_gcr_twisted.path,
					state,
					M_shifted);



	T src = zero;
	T res = zero;
			// Apply symmetric clover operator
	gaussian(src, rb[1]);

	(*twisted_gcr_solver)(res, src);

	T M_res=zero;
	(*M_shifted)(M_res,res,PLUS);
	M_res[rb[1]] -= src;
	Double norm_diff = norm2(M_res,rb[1]);
	QDPIO::cout << " || b - Mx || = " << sqrt(norm_diff) << std::endl;
	Double rel_norm_diff = sqrt(norm_diff/norm2(src,rb[1]));
	QDPIO::cout <<  " || b - Mx ||/ b= " << rel_norm_diff << std::endl;

	ASSERT_LT( toDouble(rel_norm_diff), 1.0e-8);


}

TEST_F(QudaFixture, TestShiftedSymmGCRMdagMSolverSolve )
{
	// Create Twisted Op
	Real twist = Real(0.05);
	Handle<LinearOperator<T>> M_shifted(
	  new TwistedShiftedLinOp<T, P, Q, SymEvenOddPrecLogDetLinearOperator>(*M_symm, twist));

	// Create a Group XML from the solver input XML
	std::istringstream inv_param_twisted_gcr_stream(inv_param_quda_gcr_twisted_xml);
	XMLReader xml_in_gcr_twisted(inv_param_twisted_gcr_stream);
	GroupXML_t inv_param_gcr_twisted = readXMLGroup(xml_in_gcr_twisted, "//InvertParam", "invType");

	// Create a soplver from the factory using the twisted op.
	std::istringstream  xml(inv_param_gcr_twisted.xml);
	XMLReader  paramtop(xml);

	Handle<SystemSolver<T>>	twisted_gcr_solver=
			TheMdagMFermSystemSolverFactory::Instance().createObject(inv_param_gcr_twisted.id,
					paramtop,
					inv_param_gcr_twisted.path,
					state,
					M_shifted);



	T src = zero;
	T res = zero;
			// Apply symmetric clover operator
	gaussian(src, rb[1]);

	(*twisted_gcr_solver)(res, src);

	T tmp=zero;
	T M_res=zero;

	(*M_shifted)(tmp,res,PLUS);
	(*M_shifted)(M_res,tmp,MINUS);

	M_res[rb[1]] -= src;
	Double norm_diff = norm2(M_res,rb[1]);
	QDPIO::cout << " || b - M^dag M x || = " << sqrt(norm_diff) << std::endl;
	Double rel_norm_diff = sqrt(norm_diff/norm2(src,rb[1]));
	QDPIO::cout <<  " || b - M^dag M x ||/ b= " << rel_norm_diff << std::endl;

	ASSERT_LT( toDouble(rel_norm_diff), 1.0e-8);


}

TEST_F(QudaFixture, TestShiftedAsymmGCRMdagMSolverSolve )
{
	// Create Twisted Op
	Real twist = Real(0.05);
	Handle<LinearOperator<T>> M_shifted(
	  new TwistedShiftedLinOp<T, P, Q, EvenOddPrecLinearOperator>(*M_asymm, twist));

	// Create a Group XML from the solver input XML
	std::istringstream inv_param_twisted_gcr_stream(inv_param_quda_gcr_twisted_asymm_xml);
	XMLReader xml_in_gcr_twisted(inv_param_twisted_gcr_stream);
	GroupXML_t inv_param_gcr_twisted = readXMLGroup(xml_in_gcr_twisted, "//InvertParam", "invType");

	// Create a soplver from the factory using the twisted op.
	std::istringstream  xml(inv_param_gcr_twisted.xml);
	XMLReader  paramtop(xml);

	Handle<SystemSolver<T>>	twisted_gcr_solver=
			TheMdagMFermSystemSolverFactory::Instance().createObject(inv_param_gcr_twisted.id,
					paramtop,
					inv_param_gcr_twisted.path,
					state,
					M_shifted);



	T src = zero;
	T res = zero;
			// Apply symmetric clover operator
	gaussian(src, rb[1]);

	(*twisted_gcr_solver)(res, src);

	T tmp=zero;
	T M_res=zero;

	(*M_shifted)(tmp,res,PLUS);
	(*M_shifted)(M_res,tmp,MINUS);

	M_res[rb[1]] -= src;
	Double norm_diff = norm2(M_res,rb[1]);
	QDPIO::cout << " || b - M^dag M x || = " << sqrt(norm_diff) << std::endl;
	Double rel_norm_diff = sqrt(norm_diff/norm2(src,rb[1]));
	QDPIO::cout <<  " || b - M^dag M x ||/ b= " << rel_norm_diff << std::endl;

	ASSERT_LT( toDouble(rel_norm_diff), 1.0e-8);


}

// Check QProp Functionality.
TEST_F(QudaFixture, CheckQprop)
{

	std::istringstream twisted_fermact(fermact_xml_symm_twisted);
	XMLReader xml_in_twisted(twisted_fermact);
	Handle<S_symm_T> S_symm_twisted = dynamic_cast<S_symm_T*>(TheFermionActionFactory::Instance().createObject("SEOPREC_CLOVER",
				    											   xml_in_twisted,
				    											   "FermionAction"));

	std::istringstream inv_param_xml_stream(inv_param_quda_gcr_twisted_xml);
	XMLReader xml_in(inv_param_xml_stream);

	GroupXML_t inv_param = readXMLGroup(xml_in, "//InvertParam", "invType");
	Handle<SystemSolver<T>>	qprop_solver = S_symm_twisted->qprop(state,inv_param);

    Handle<SymEvenOddPrecLinearOperator<T,P,Q>> M = S_symm_twisted->linOp(state);

    LatticeFermion rhs=zero;
    LatticeFermion x = zero;
    gaussian(rhs);
	(*qprop_solver)(x,rhs);

	// Check residuum
	LatticeFermion Ax=zero;
	//M_u(Ax,x,PLUS);
	M->unprecLinOp(Ax,x,PLUS);
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

