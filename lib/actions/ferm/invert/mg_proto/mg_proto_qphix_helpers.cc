/*
 * mg_proto_helpers.cpp
 *
 *  Created on: Mar 24, 2017
 *      Author: bjoo
 */
#include "chromabase.h"
#include "meas/inline/io/named_objmap.h"
#include "actions/ferm/invert/mg_proto/mg_proto_qphix_helpers.h"
#include "lattice/fine_qdpxx/mg_params_qdpxx.h"
#include "lattice/qphix/mg_level_qphix.h"
#include "lattice/qphix/vcycle_recursive_qphix.h"
#include "lattice/qphix/qphix_clover_linear_operator.h"

#include <vector>
#include <memory>

using std::shared_ptr;
using std::make_shared;
using std::vector;

using namespace QDP;


namespace Chroma {
namespace  MGProtoHelpersQPhiX{

template<typename QPhiXLinOpT>
shared_ptr<QPhiXLinOpT>
createFineLinOpT( const MGProtoSolverParams& params, const multi1d<LatticeColorMatrix>& u,
    const MG::LatticeInfo& info)
{
	shared_ptr<QPhiXLinOpT> M_fine=nullptr;

	bool anisoP = params.CloverParams.anisoParam.anisoP;
	double m_q=toDouble(params.CloverParams.Mass);
	double u0 = toDouble(params.CloverParams.u0);

#if 0
	int t_bc = params.AntiPeriodicT ? -1 : 1; - ignore as links already have BC on.
#else
	int t_bc = 1;
#endif


	if( !anisoP ) {
		double c_sw=toDouble(params.CloverParams.clovCoeffR);
		M_fine = make_shared<QPhiXLinOpT>(info,m_q,c_sw,t_bc,u);
	}
	else {
		QDPIO::cout << "Using aniso interface" << std::endl;

		double xi0 = toDouble(params.CloverParams.anisoParam.xi_0);
		double nu=toDouble(params.CloverParams.anisoParam.nu);
		double c_sw_r = toDouble(params.CloverParams.clovCoeffR);
		double c_sw_t = toDouble(params.CloverParams.clovCoeffT);
		M_fine = make_shared<QPhiXLinOpT>(info,m_q,u0,xi0,nu,c_sw_r,c_sw_t,t_bc,u);
	}

	return M_fine;
}

std::shared_ptr<MG::QPhiXWilsonCloverLinearOperatorF>
createFineLinOpF( const MGProtoSolverParams& params, const multi1d<LatticeColorMatrix>& u,
    const MG::LatticeInfo& info)
{
  return createFineLinOpT<MG::QPhiXWilsonCloverLinearOperatorF>(params,u,info);
}

std::shared_ptr<MG::QPhiXWilsonCloverLinearOperator>
createFineLinOp( const MGProtoSolverParams& params, const multi1d<LatticeColorMatrix>& u,
    const MG::LatticeInfo& info)
{
  return createFineLinOpT<MG::QPhiXWilsonCloverLinearOperator>(params,u,info);
}

void
createMGPreconditioner( const MGProtoSolverParams& params, const multi1d<LatticeColorMatrix>& u)
{
	START_CODE();
	StopWatch swatch;
	swatch.reset();
	swatch.start();

	const std::string& subspaceId = params.SubspaceId;
	QDPIO::cout << "Creating MG_Proto preconditioner with subspaceID=" << subspaceId << std::endl;
	// Look up in named object map
	QDPIO::cout << "Looking Up " << subspaceId << " in the Named Onject Map" << std::endl;
	if( TheNamedObjMap::Instance().check(subspaceId) ) {
		// Found it... Delete it.
		QDPIO::cout << "  ... Subspace ID found... Deleting" <<std::endl;
		deleteMGPreconditioner(subspaceId);
	}


	// Now make a new one.
	shared_ptr<MG::QPhiXMultigridLevels> mg_levels = make_shared<MG::QPhiXMultigridLevels>();

	// First make an info from the lattice parameters:
  IndexArray latdims = {{ QDP::Layout::subgridLattSize()[0],
                QDP::Layout::subgridLattSize()[1],
                QDP::Layout::subgridLattSize()[2],
                QDP::Layout::subgridLattSize()[3] }};

  (mg_levels->fine_level).info = std::make_shared<LatticeInfo>( latdims, 4,3,NodeInfo());

	// First make M
	QDPIO::cout << "Creating M..." ;
	shared_ptr<MG::QPhiXWilsonCloverLinearOperator> M=createFineLinOp( params, u, *((mg_levels->fine_level).info) );
	QDPIO::cout << "Done" << std::endl;

	QDPIO::cout << "Creating single prec M...";
	shared_ptr<MG::QPhiXWilsonCloverLinearOperatorF> M_f=createFineLinOpF( params, u, *((mg_levels->fine_level).info) );
	QDPIO::cout << "Done" << std::endl;

	QDPIO::cout << "Creating MG Levels ... " ;
	MG::SetupParams level_params;
	int n_levels = params.MGLevels;
	level_params.n_levels = n_levels;
	level_params.n_vecs.resize(n_levels-1);
	level_params.null_solver_max_iter.resize(n_levels-1);
	level_params.null_solver_rsd_target.resize(n_levels -1);
	level_params.null_solver_verboseP.resize(n_levels -1);
	for(int l=0; l < n_levels-1;++l) {
		QDPIO::cout << "Level L=" << l << " Null Vecs=" << params.NullVecs[l] << std::endl;

		level_params.n_vecs[l] = params.NullVecs[l];
		level_params.null_solver_max_iter[l]=params.NullSolverMaxIters[l];
		level_params.null_solver_rsd_target[l]=toDouble(params.NullSolverRsdTarget[l]);
		level_params.null_solver_verboseP[l]=toDouble(params.NullSolverVerboseP[l]);
	}
	level_params.block_sizes.resize(n_levels-1);
	for(int l=0; l < n_levels-1;++l) {
		for(int mu=0; mu < 4; ++mu) {
			(level_params.block_sizes[l])[mu] = (params.Blocking[l])[mu];
		}
	}

	MG::SetupQPhiXMGLevels(level_params, *mg_levels, M_f );
	QDPIO::cout << "... Done " << std::endl;

	if( (mg_levels->fine_level).M == nullptr ) {
	  QDPIO::cout << "Error... Barfaroni. Fine Level M is null after subspace setup" << std::endl;
	  QDP_abort(1);
	}

	QDPIO::cout << "Creating VCycle Parameters..." << std::endl;
	vector<MG::VCycleParams> v_params(n_levels-1);
	for(int l=0; l < n_levels-1;++l) {
		QDPIO::cout << "   Level " << l << std::endl;
		v_params[l].pre_smoother_params.MaxIter=params.VCyclePreSmootherMaxIters[l];
		v_params[l].pre_smoother_params.RsdTarget=toDouble(params.VCyclePreSmootherRsdTarget[l]);
		v_params[l].pre_smoother_params.VerboseP =params.VCyclePreSmootherVerboseP[l];
		v_params[l].pre_smoother_params.Omega =toDouble(params.VCyclePreSmootherRelaxOmega[l]);

		v_params[l].post_smoother_params.MaxIter=params.VCyclePostSmootherMaxIters[l];
		v_params[l].post_smoother_params.RsdTarget=toDouble(params.VCyclePostSmootherRsdTarget[l]);
		v_params[l].post_smoother_params.VerboseP =params.VCyclePostSmootherVerboseP[l];
		v_params[l].post_smoother_params.Omega =toDouble(params.VCyclePostSmootherRelaxOmega[l]);

		v_params[l].bottom_solver_params.MaxIter=params.VCycleBottomSolverMaxIters[l];
		v_params[l].bottom_solver_params.NKrylov = params.VCycleBottomSolverNKrylov[l];
		v_params[l].bottom_solver_params.RsdTarget= toDouble(params.VCycleBottomSolverRsdTarget[l]);
		v_params[l].bottom_solver_params.VerboseP = params.VCycleBottomSolverVerboseP[l];

		v_params[l].cycle_params.MaxIter=params.VCycleMaxIters[l];
		v_params[l].cycle_params.RsdTarget=toDouble(params.VCycleRsdTarget[l]);
		v_params[l].cycle_params.VerboseP = params.VCycleVerboseP[l];
	}

	QDPIO::cout << "Creating VCycle Preconditioner...";

	shared_ptr<MG::VCycleRecursiveQPhiX> v_cycle=make_shared<MG::VCycleRecursiveQPhiX>(v_params, *mg_levels);

	QDPIO::cout << "Done";

	QDPIO::cout << "Saving in Map" << std::endl;
	QDPIO::cout << "Creating Named Object Map Entry for subspace" << std::endl;
	XMLBufferWriter file_xml;
	push(file_xml, "FileXML");
	pop(file_xml);

	XMLBufferWriter record_xml;
	push(record_xml, "RecordXML");
	write(record_xml, "InvertParam", params);
	pop(record_xml);

	TheNamedObjMap::Instance().create<shared_ptr<MGPreconditioner>>(subspaceId);
	TheNamedObjMap::Instance().get(subspaceId).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(subspaceId).setRecordXML(record_xml);
	TheNamedObjMap::Instance().getData<shared_ptr<MGPreconditioner>>(subspaceId)=make_shared<MGPreconditioner>();
	TheNamedObjMap::Instance().getData<shared_ptr<MGPreconditioner>>(subspaceId)->mg_levels = mg_levels;
	TheNamedObjMap::Instance().getData<shared_ptr<MGPreconditioner>>(subspaceId)->v_cycle = v_cycle;
	TheNamedObjMap::Instance().getData<shared_ptr<MGPreconditioner>>(subspaceId)->M = M;

	swatch.stop();
	QDPIO::cout << "MG_PROTO_QPHIX_SETUP: Subspace Creation Took : " << swatch.getTimeInSeconds() << " sec" << std::endl;

	END_CODE();
}

void deleteMGPreconditioner(const std::string& subspaceId)
{

	QDPIO::cout << "Deleting MG_Proto preconditioner with subspaceID=" << subspaceId << std::endl;
	// Look up in named object map
	QDPIO::cout << "Looking Up " << subspaceId << " in the Named Onject Map" << std::endl;
	if( TheNamedObjMap::Instance().check(subspaceId) ) {
		// Found it... Delete it.
		QDPIO::cout << "  ... Subspace ID found... Deleting" <<std::endl;

		// This will erase the MGPreconditioner, which has in it shared pointers.
		// If the shared pointers are destroyed, and there are no references remaining
		// then MG Levels will be destroyed as weill v_cycle
		// These only hold allocated data by shared pointer, so they too should be cleaned
		TheNamedObjMap::Instance().erase(subspaceId);
	}
}

shared_ptr<MGPreconditioner> getMGPreconditioner(const std::string& subspaceId)
{
	shared_ptr<MGPreconditioner> ret_val = nullptr;
	if( TheNamedObjMap::Instance().check(subspaceId) ) {
		// Found it... Delete it.
		QDPIO::cout << "  ... Subspace ID found... returning" <<std::endl;

		// This will erase the MGPreconditioner, which has in it shared pointers.
		// If the shared pointers are destroyed, and there are no references remaining
		// then MG Levels will be destroyed as weill v_cycle
		// These only hold allocated data by shared pointer so they too should be cleaned
		ret_val = TheNamedObjMap::Instance().getData<shared_ptr<MGPreconditioner>>(subspaceId);

	}
	else {
		QDPIO::cout << "Object Not Found... Returning Null" << std::endl;
		ret_val = nullptr;

	}
	return ret_val;

}




}

}
