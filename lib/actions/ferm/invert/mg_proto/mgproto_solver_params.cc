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

// Anonymous namespace
namespace  {
 template<typename T>
 void read(XMLReader& xml, const std::string& path, multi1d<T>& result, int N)
 {
	 multi1d<T> read_result;
	 read(xml, path, read_result);
	 result.resize(N);

	 if( read_result.size() == 1 ) {
		 // Broadcast
		 for(int i=0; i < N; ++i) {
			 result[i] = read_result[0];
		 }
	 }
	 else {
		 if( read_result.size() == N ) {
			 // Copy
			 for(int i=0; i < N; ++i) {
				 result[i] = read_result[i];
			 }
		 }
		 else {
			 QDPIO::cerr << "Wrong number of elements reading" << path << " should be either 1 or " << N << std::endl;
			 QDP_abort(1);
		 }
	 }
 }

} // end anonymous namespace
MGProtoSolverParams::MGProtoSolverParams(XMLReader& xml, const std::string& path)
{
	try {
	XMLReader paramtop(xml, path);
	read( paramtop, "CloverParams", CloverParams);

	read( paramtop, "AntiPeriodicT", AntiPeriodicT );
	read( paramtop, "MGLevels", MGLevels);

	read( paramtop, "Blocking", Blocking, MGLevels-1);

	read( paramtop, "NullVecs", NullVecs, MGLevels-1);
	read( paramtop, "NullSolverMaxIters", NullSolverMaxIters, MGLevels-1);
	read( paramtop, "NullSolverRsdTarget", NullSolverRsdTarget, MGLevels-1);
	read( paramtop, "NullSolverVerboseP", NullSolverVerboseP, MGLevels-1);

	read( paramtop, "OuterSolverNKrylov", OuterSolverNKrylov);
	read( paramtop, "OuterSolverRsdTarget", OuterSolverRsdTarget);
	read( paramtop, "OuterSolverMaxIters", OuterSolverMaxIters);
	read( paramtop, "OuterSolverVerboseP", OuterSolverVerboseP);
    if( paramtop.count("RsdToleranceFactor") > 0 ) {
        read(paramtop, "RsdToleranceFactor", RsdToleranceFactor);
    } else {
        RsdToleranceFactor = 10;
    }

    if( paramtop.count("ThresholdCount") > 0 ) {
        read(paramtop, "ThresholdCount", ThresholdCount);
    } else {
        ThresholdCount = 50; // current setup for test
    }

	read( paramtop, "VCyclePreSmootherMaxIters", VCyclePreSmootherMaxIters, MGLevels-1);
	read( paramtop, "VCyclePreSmootherRsdTarget", VCyclePreSmootherRsdTarget, MGLevels-1);
	read( paramtop, "VCyclePreSmootherRelaxOmega", VCyclePreSmootherRelaxOmega, MGLevels-1);
	read( paramtop, "VCyclePreSmootherVerboseP", VCyclePreSmootherVerboseP, MGLevels-1);

	read( paramtop, "VCyclePostSmootherMaxIters", VCyclePostSmootherMaxIters, MGLevels-1);
	read( paramtop, "VCyclePostSmootherRsdTarget", VCyclePostSmootherRsdTarget, MGLevels-1);
	read( paramtop, "VCyclePostSmootherRelaxOmega", VCyclePostSmootherRelaxOmega, MGLevels-1);
	read( paramtop, "VCyclePostSmootherVerboseP", VCyclePostSmootherVerboseP, MGLevels-1);

	read( paramtop, "VCycleBottomSolverNKrylov", VCycleBottomSolverNKrylov, MGLevels-1);
	read( paramtop, "VCycleBottomSolverMaxIters", VCycleBottomSolverMaxIters, MGLevels-1);
	read( paramtop, "VCycleBottomSolverRsdTarget", VCycleBottomSolverRsdTarget,MGLevels-1);
	read( paramtop, "VCycleBottomSolverVerboseP", VCycleBottomSolverVerboseP, MGLevels-1);

	read( paramtop, "VCycleMaxIters", VCycleMaxIters, MGLevels-1);
	read( paramtop, "VCycleRsdTarget", VCycleRsdTarget, MGLevels-1);
	read( paramtop, "VCycleVerboseP", VCycleVerboseP, MGLevels-1);

	read( paramtop, "SubspaceId", SubspaceId);

	}
	catch( const std::string& e) {
		QDPIO::cout << "Caught exception " << e << std::endl;
		QDP_abort(1);
	}
	catch(...) {
		QDPIO::cout << "Caught unknown exception " << std::endl;
		QDP_abort(1);
	}
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
	write(xml, "AntiPeriodicT", p.AntiPeriodicT);

	write(xml, "MGLevels", p.MGLevels);
	write(xml, "Blocking", p.Blocking);
	write(xml, "NullVecs", p.NullVecs);
	write(xml, "NullSolverMaxIters", p.NullSolverMaxIters);
	write(xml, "NullSolverRsdTarget", p.NullSolverRsdTarget);
	write(xml, "NullSolverVerboseP", p.NullSolverVerboseP);

	write(xml, "OuterSolverNKrylov", p.OuterSolverNKrylov);
	write(xml, "OuterSolverRsdTarget", p.OuterSolverRsdTarget);
	write(xml, "OuterSolverMaxIters", p.OuterSolverMaxIters);
	write(xml, "OuterSolverVerboseP", p.OuterSolverVerboseP);

	write(xml, "VCyclePreSmootherMaxIters", p.VCyclePreSmootherMaxIters);
	write(xml, "VCyclePreSmootherRsdTarget", p.VCyclePreSmootherRsdTarget);
	write(xml, "VCyclePreSmootherRelaxOmega", p.VCyclePreSmootherRelaxOmega);
	write(xml, "VCyclePreSmootherVerboseP", p.VCyclePreSmootherVerboseP);

	write(xml, "VCyclePostSmootherMaxIters", p.VCyclePostSmootherMaxIters);
	write(xml, "VCyclePostSmootherRsdTarget", p.VCyclePostSmootherRsdTarget);
	write(xml, "VCyclePostSmootherRelaxOmega", p.VCyclePostSmootherRelaxOmega);
	write(xml, "VCyclePostSmootherVerboseP", p.VCyclePostSmootherVerboseP);

	write(xml, "VCycleBottomSolverNKrylov", p.VCycleBottomSolverNKrylov);
	write(xml, "VCycleBottomSolverMaxIters", p.VCycleBottomSolverMaxIters);
	write(xml, "VCycleBottomSolverRsdTarget", p.VCycleBottomSolverRsdTarget);
	write(xml, "VCycleBottomSolverVerboseP", p.VCycleBottomSolverVerboseP);

	write(xml, "VCycleMaxIters", p.VCycleMaxIters);
	write(xml, "VCycleRsdTarget", p.VCycleMaxIters);
	write(xml, "VCycleVerboseP", p.VCycleVerboseP);
	write(xml, "SubspaceId", p.SubspaceId);

}


}
