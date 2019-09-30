/*
 * mg_proto_helpers.h
 *
 *  Created on: Mar 24, 2017
 *      Author: bjoo
 */

#ifndef LIB_ACTIONS_FERM_INVERT_MG_PROTO_MG_PROTO_QPHIX_HELPERS_H_
#define LIB_ACTIONS_FERM_INVERT_MG_PROTO_MG_PROTO_QPHIX_HELPERS_H_

#include <memory>
#include <string>
#include <lattice/lattice_info.h>
#include <lattice/qphix/mg_level_qphix.h>
#include <lattice/qphix/vcycle_recursive_qphix.h>
#include <actions/ferm/invert/mg_proto/mgproto_solver_params.h>
#include <lattice/lattice_info.h>
#include <lattice/qphix/qphix_clover_linear_operator.h>
#include <lattice/qphix/qphix_eo_clover_linear_operator.h>

namespace Chroma {

namespace MGProtoHelpersQPhiX {

template<typename MyLevelT, typename MyVCycleT, typename MyLinOpT, typename MyLinOpFT>
struct MGPreconditionerT
{
	using LevelT = MyLevelT;
	using VCycleT = MyVCycleT;
	using LinOpT = MyLinOpT;
	using LinOpFT = MyLinOpFT;

	std::shared_ptr<LevelT> mg_levels;
	std::shared_ptr<VCycleT> v_cycle;
	std::shared_ptr<LinOpT> M;
};

using MGPreconditioner = MGPreconditionerT<MG::QPhiXMultigridLevels, MG::VCycleRecursiveQPhiX,MG::QPhiXWilsonCloverLinearOperator, MG::QPhiXWilsonCloverLinearOperatorF>;
using MGPreconditionerEO = MGPreconditionerT<MG::QPhiXMultigridLevelsEO, MG::VCycleRecursiveQPhiXEO2,MG::QPhiXWilsonCloverEOLinearOperator, MG::QPhiXWilsonCloverEOLinearOperatorF>;


// for testing
template<typename QPhiXLinOpT>
std::shared_ptr<QPhiXLinOpT>
createFineLinOpT( const MGProtoSolverParams& params, const multi1d<LatticeColorMatrix>& u,
    const MG::LatticeInfo& info);

std::shared_ptr<MG::QPhiXWilsonCloverLinearOperator>
createFineLinOp( const MGProtoSolverParams& params, const multi1d<LatticeColorMatrix>& u,
    const MG::LatticeInfo& info);

std::shared_ptr<MG::QPhiXWilsonCloverLinearOperatorF>
createFineLinOpF( const MGProtoSolverParams& params, const multi1d<LatticeColorMatrix>& u,
    const MG::LatticeInfo& info);

void createMGPreconditioner(const MGProtoSolverParams& params, const multi1d<LatticeColorMatrix>& u);
void deleteMGPreconditioner(const std::string& subspaceID);
std::shared_ptr<MGPreconditioner> getMGPreconditioner(const std::string& subspaceId);


// EO versions
std::shared_ptr<MG::QPhiXWilsonCloverEOLinearOperator>
createFineEOLinOp( const MGProtoSolverParams& params, const multi1d<LatticeColorMatrix>& u,
    const MG::LatticeInfo& info);

std::shared_ptr<MG::QPhiXWilsonCloverEOLinearOperatorF>
createFineEOLinOpF( const MGProtoSolverParams& params, const multi1d<LatticeColorMatrix>& u,
    const MG::LatticeInfo& info);

void createMGPreconditionerEO(const MGProtoSolverParams& params, const multi1d<LatticeColorMatrix>& u);
void deleteMGPreconditionerEO(const std::string& subspaceID);
std::shared_ptr<MGPreconditionerEO> getMGPreconditionerEO(const std::string& subspaceId);

}  // Namespace MGProtoHelpers


} // Namespace Chroma






#endif /* LIB_ACTIONS_FERM_INVERT_MG_PROTO_MG_PROTO_HELPERS_H_ */
