/*
 * mgproto_solver_params.cc
 *
 *  Created on: Mar 23, 2017
 *      Author: bjoo
 */

#include "chromabase.h"
#include "actions/ferm/invert/mg_proto/mgproto_solver_params.h"
#include "io/xml_group_reader.h"
#include <string>

using std::string;

using namespace QDP;

namespace Chroma
{

MGProtoSolverParams::MGProtoSolverParams(XMLReader& xml, const std::string& path)
{
	XMLReader paramtop(xml, path);
	read( paramtop, "CloverParams", CloverParams);
	read( paramtop, "MGLevels", MGLevels);
	read( paramtop, "Blocking", Blocking);
	read( paramtop, "NullVecs", NullVecs);
	read( paramtop, "NullSolverMaxIters", NullSolverMaxIters);
	read( paramtop, "NullSolverRsdTarget", NullSolverRsdTarget);

	read( paramtop, "OuterSolverNKrylov", OuterSolverNKrylov);
	read( paramtop, "OuterSolverRsdTarget", OuterSolverRsdTarget);
	read( paramtop, "OuterSolverMaxIters", OuterSolverMaxIters);

	read( paramtop, "VCyclePreSmootherMaxIters", VCyclePreSmootherMaxIters);
	read( paramtop, "VCyclePreSmootehrRsdTarget", VCyclePreSmootherRsdTarget);
	read( paramtop, "VCyclePreSmootehrRelaxOmega", VCyclePreSmootherRelaxOmega);

	read( paramtop, "VCyclePostSmootherMaxIters", VCyclePostSmootherMaxIters);
	read( paramtop, "VCyclePostSmootehrRsdTarget", VCyclePostSmootherRsdTarget);
	read( paramtop, "VCyclePostSmootehrRelaxOmega", VCyclePostSmootherRelaxOmega);

	read( paramtop, "VCycleBottomSolverNKrylov", VCycleBottomSolverNKrylov);
	read( paramtop, "VCycleBottomSolverMaxIters", VCycleBottomSolverMaxIters);
	read( paramtop, "VCycleBottomSolverRsdTarget", VCycleBottomSolverRsdTarget);

	read( paramtop, "VCycleMaxIters", VCycleMaxIters);
	read( paramtop, "VCycleRsdTarget", VCycleRsdTarget);
	read( paramtop, "SubspaceID", SubspaceId);
}


void read(XMLReader& xml, const std::string& path, MGProtoSolverParams& p)
{
	MGProtoSolverParams ret_val(xml,path);
	p = ret_val;
}

void write(XMLWriter& xml, const std::string& path, const MGProtoSolverParams& p)
{
	push(xml, path);
	write(xml, "CloverParams", p.CloverParams);
	write(xml, "MGLevels", p.MGLevels);
	write(xml, "Blocking", p.Blocking);
	write(xml, "NullVecs", p.NullVecs);
	write(xml, "NullSolverMaxIters", p.NullSolverMaxIters);
	write(xml, "NullSolverRsdTarget", p.NullSolverRsdTarget);
	write(xml, "OuterSolverNKrylov", p.OuterSolverNKrylov);
	write(xml, "OuterSolverRsdTarget", p.OuterSolverRsdTarget);
	write(xml, "OuterSolverMaxIters", p.OuterSolverMaxIters);

	write(xml, "VCyclePreSmootherMaxIters", p.VCyclePreSmootherMaxIters);
	write(xml, "VCyclePreSmootherRsdTarget", p.VCyclePreSmootherRsdTarget);
	write(xml, "VCyclePreSmootherRelaxOmega", p.VCyclePreSmootherRelaxOmega);

	write(xml, "VCyclePostSmootherMaxIters", p.VCyclePostSmootherMaxIters);
	write(xml, "VCyclePostSmootherRsdTarget", p.VCyclePostSmootherRsdTarget);
	write(xml, "VCyclePostSmootherRelaxOmega", p.VCyclePostSmootherRelaxOmega);

	write(xml, "VCycleBottomSolverNKrylov", p.VCycleBottomSolverNKrylov);
	write(xml, "VCycleBottomSolverMaxIters", p.VCycleBottomSolverMaxIters);
	write(xml, "VCycleBottomSolverRsdTarget", p.VCycleBottomSolverRsdTarget);

	write(xml, "VCycleMaxIters", p.VCycleMaxIters);
	write(xml, "VCycleRsdTarget", p.VCycleMaxIters);

	write(xml, "SubspaceID", p.SubspaceId);

}


}
