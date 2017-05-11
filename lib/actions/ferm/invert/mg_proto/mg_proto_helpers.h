/*
 * mg_proto_helpers.h
 *
 *  Created on: Mar 24, 2017
 *      Author: bjoo
 */

#ifndef LIB_ACTIONS_FERM_INVERT_MG_PROTO_MG_PROTO_HELPERS_H_
#define LIB_ACTIONS_FERM_INVERT_MG_PROTO_MG_PROTO_HELPERS_H_

#include <memory>
#include <string>
#include "lattice/fine_qdpxx/mg_level_qdpxx.h"
#include "lattice/fine_qdpxx/vcycle_recursive_qdpxx.h"
#include "actions/ferm/invert/mg_proto/mgproto_solver_params.h"
#include "lattice/fine_qdpxx/wilson_clover_linear_operator.h"

namespace Chroma {

namespace MGProtoHelpers {

struct MGPreconditioner
{
	std::shared_ptr<MG::MultigridLevels> mg_levels;
	std::shared_ptr<MG::VCycleRecursiveQDPXX> v_cycle;
};

// for testing
std::shared_ptr<const MG::QDPWilsonCloverLinearOperator>
createFineLinOp( const MGProtoSolverParams& params, const multi1d<LatticeColorMatrix>& u);

void createMGPreconditioner(const MGProtoSolverParams& params, const multi1d<LatticeColorMatrix>& u);
void deleteMGPreconditioner(const std::string& subspaceID);
std::shared_ptr<MGPreconditioner> getMGPreconditioner(const std::string& subspaceId);

}  // Namespace MGProtoHelpers


} // Namespace Chroma






#endif /* LIB_ACTIONS_FERM_INVERT_MG_PROTO_MG_PROTO_HELPERS_H_ */
