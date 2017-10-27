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

namespace Chroma {

namespace MGProtoHelpersQPhiX {

struct MGPreconditioner
{
	std::shared_ptr<MG::QPhiXMultigridLevels> mg_levels;
	std::shared_ptr<MG::VCycleRecursiveQPhiX> v_cycle;
	std::shared_ptr<MG::QPhiXWilsonCloverLinearOperator> M;

};

// for testing
std::shared_ptr<MG::QPhiXWilsonCloverLinearOperator>
createFineLinOp( const MGProtoSolverParams& params, const multi1d<LatticeColorMatrix>& u,
    const MG::LatticeInfo& info);

std::shared_ptr<MG::QPhiXWilsonCloverLinearOperatorF>
createFineLinOpF( const MGProtoSolverParams& params, const multi1d<LatticeColorMatrix>& u,
    const MG::LatticeInfo& info);

void createMGPreconditioner(const MGProtoSolverParams& params, const multi1d<LatticeColorMatrix>& u);
void deleteMGPreconditioner(const std::string& subspaceID);
std::shared_ptr<MGPreconditioner> getMGPreconditioner(const std::string& subspaceId);

}  // Namespace MGProtoHelpers


} // Namespace Chroma






#endif /* LIB_ACTIONS_FERM_INVERT_MG_PROTO_MG_PROTO_HELPERS_H_ */
