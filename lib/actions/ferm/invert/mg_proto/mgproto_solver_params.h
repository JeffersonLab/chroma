/*
 * mgproto_solver_params.h
 *
 *  Created on: Mar 23, 2017
 *      Author: bjoo
 */

#ifndef LIB_ACTIONS_FERM_INVERT_MG_PROTO_MGPROTO_SOLVER_PARAMS_H
#define LIB_ACTIONS_FERM_INVERT_MG_PROTO_MGPROTO_SOLVER_PARAMS_H

#include "chromabase.h"
#include "io/xml_group_reader.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"

#include <string>

using namespace QDP;

namespace Chroma
{

struct MGProtoSolverParams {
public:
	MGProtoSolverParams(){ }

	MGProtoSolverParams(XMLReader& xml, const std::string& path);

	// Clover Params
	Chroma::CloverFermActParams CloverParams;
	bool AntiPeriodicT;
	// Details of the NullSpace Construction
	int MGLevels;
	multi1d< multi1d<int> > Blocking;
	multi1d< int > NullVecs;
	multi1d< int > NullSolverMaxIters;
	multi1d< Double > NullSolverRsdTarget;
	multi1d< bool > NullSolverVerboseP;

	// Details of Outer Solver
	int OuterSolverNKrylov;
	Double OuterSolverRsdTarget;
    Double RsdToleranceFactor;
	int OuterSolverMaxIters;
	bool OuterSolverVerboseP;
    int ThresholdCount;

	// VCycle Details (MGLevels - 1 of these)
	// Presmoother
	multi1d<int> VCyclePreSmootherMaxIters;
	multi1d<Double> VCyclePreSmootherRsdTarget;
	multi1d<bool> VCyclePreSmootherVerboseP;
	multi1d<Real> VCyclePreSmootherRelaxOmega;

	// Post Smoother
	multi1d<int> VCyclePostSmootherMaxIters;
	multi1d<Double> VCyclePostSmootherRsdTarget;
	multi1d<bool> VCyclePostSmootherVerboseP;
	multi1d<Real> VCyclePostSmootherRelaxOmega;

	// Bottom Solver
	multi1d<int> VCycleBottomSolverNKrylov;
	multi1d<int> VCycleBottomSolverMaxIters;
	multi1d<Double> VCycleBottomSolverRsdTarget;
	multi1d<bool> VCycleBottomSolverVerboseP;

	// VCycle Iterations/Target Residuum
	multi1d<int> VCycleMaxIters;
	multi1d<Double> VCycleRsdTarget;
	multi1d<bool> VCycleVerboseP;

	std::string SubspaceId;

};

void read(XMLReader& xml, const std::string& path, MGProtoSolverParams& p);
void write(XMLWriter& xml, const std::string& path, const MGProtoSolverParams& p);

};


#endif


