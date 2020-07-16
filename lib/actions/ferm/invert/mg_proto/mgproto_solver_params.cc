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

MGProtoMGDeflationParams::MGProtoMGDeflationParams(XMLReader& xml, const std::string& path)
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

	read( paramtop, "EigenSolverMaxRestartSize", EigenSolverMaxRestartSize);
	read( paramtop, "EigenSolverMaxRank", EigenSolverMaxRank);
	read( paramtop, "EigenSolverRsdTarget", EigenSolverRsdTarget);
	read( paramtop, "EigenSolverMaxIters", EigenSolverMaxIters);
	read( paramtop, "EigenSolverVerboseP", EigenSolverVerboseP);

	read( paramtop, "BottomSolverNKrylov", BottomSolverNKrylov);
	read( paramtop, "BottomSolverMaxIters", BottomSolverMaxIters);
	read( paramtop, "BottomSolverRsdTarget", BottomSolverRsdTarget);
	read( paramtop, "BottomSolverVerboseP", BottomSolverVerboseP);

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


void read(XMLReader& xml, const std::string& path, MGProtoMGDeflationParams& p)
{
	MGProtoMGDeflationParams ret_val(xml,path);
	p = ret_val;
}

void write(XMLWriter& xml, const std::string& path, const MGProtoMGDeflationParams& p)
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

	write(xml, "EigenSolverMaxRestartSize", p.EigenSolverMaxRestartSize);
	write(xml, "EigenSolverMaxRank", p.EigenSolverMaxRank);
	write(xml, "EigenSolverRsdTarget", p.EigenSolverRsdTarget);
	write(xml, "EigenSolverMaxIters", p.EigenSolverMaxIters);
	write(xml, "EigenSolverVerboseP", p.EigenSolverVerboseP);

	write(xml, "BottomSolverNKrylov", p.BottomSolverNKrylov);
	write(xml, "BottomSolverMaxIters", p.BottomSolverMaxIters);
	write(xml, "BottomSolverRsdTarget", p.BottomSolverRsdTarget);
	write(xml, "BottomSolverVerboseP", p.BottomSolverVerboseP);

	write(xml, "SubspaceId", p.SubspaceId);

}

MGProtoALIPrecParams::MGProtoALIPrecParams(XMLReader& xml, const std::string& path)
{
	try {
	XMLReader paramtop(xml, path);
	read( paramtop, "CloverParams", CloverParams);
	read( paramtop, "AntiPeriodicT", AntiPeriodicT );

        XMLReader defl(paramtop, "MGDeflation");
	read( defl, "MGLevels", Deflation.MGLevels);

	read( defl, "Blocking", Deflation.Blocking, Deflation.MGLevels-1);

	read( defl, "NullVecs", Deflation.NullVecs, Deflation.MGLevels-1);
	read( defl, "NullSolverMaxIters", Deflation.NullSolverMaxIters, Deflation.MGLevels-1);
	read( defl, "NullSolverRsdTarget", Deflation.NullSolverRsdTarget, Deflation.MGLevels-1);
	read( defl, "NullSolverVerboseP", Deflation.NullSolverVerboseP, Deflation.MGLevels-1);

	read( defl, "EigenSolverMaxRestartSize", Deflation.EigenSolverMaxRestartSize);
	read( defl, "EigenSolverMaxRank", Deflation.EigenSolverMaxRank);
	read( defl, "EigenSolverRsdTarget", Deflation.EigenSolverRsdTarget);
	read( defl, "EigenSolverMaxIters", Deflation.EigenSolverMaxIters);
	read( defl, "EigenSolverVerboseP", Deflation.EigenSolverVerboseP);

	read( defl, "BottomSolverNKrylov", Deflation.BottomSolverNKrylov);
	read( defl, "BottomSolverMaxIters", Deflation.BottomSolverMaxIters);
	read( defl, "BottomSolverRsdTarget", Deflation.BottomSolverRsdTarget);
	read( defl, "BottomSolverVerboseP", Deflation.BottomSolverVerboseP);

        XMLReader recon(paramtop, "Reconstruction");
	read( recon, "ali_distance", Reconstruction.ali_distance);
	read( recon, "probing_distance", Reconstruction.probing_distance);

	read( recon, "MGLevels", Reconstruction.MGLevels);

	read( recon, "Blocking", Reconstruction.Blocking, Reconstruction.MGLevels-1);

	read( recon, "NullVecs", Reconstruction.NullVecs, Reconstruction.MGLevels-1);
	read( recon, "NullSolverMaxIters", Reconstruction.NullSolverMaxIters, Reconstruction.MGLevels-1);
	read( recon, "NullSolverRsdTarget", Reconstruction.NullSolverRsdTarget, Reconstruction.MGLevels-1);
	read( recon, "NullSolverVerboseP", Reconstruction.NullSolverVerboseP, Reconstruction.MGLevels-1);

	read( recon, "OuterSolverNKrylov", Reconstruction.OuterSolverNKrylov);
	read( recon, "OuterSolverRsdTarget", Reconstruction.OuterSolverRsdTarget);
	read( recon, "OuterSolverMaxIters", Reconstruction.OuterSolverMaxIters);
	read( recon, "OuterSolverVerboseP", Reconstruction.OuterSolverVerboseP);

	read( recon, "VCyclePreSmootherMaxIters", Reconstruction.VCyclePreSmootherMaxIters, Reconstruction.MGLevels-1);
	read( recon, "VCyclePreSmootherRsdTarget", Reconstruction.VCyclePreSmootherRsdTarget, Reconstruction.MGLevels-1);
	read( recon, "VCyclePreSmootherRelaxOmega", Reconstruction.VCyclePreSmootherRelaxOmega, Reconstruction.MGLevels-1);
	read( recon, "VCyclePreSmootherVerboseP", Reconstruction.VCyclePreSmootherVerboseP, Reconstruction.MGLevels-1);

	read( recon, "VCyclePostSmootherMaxIters", Reconstruction.VCyclePostSmootherMaxIters, Reconstruction.MGLevels-1);
	read( recon, "VCyclePostSmootherRsdTarget", Reconstruction.VCyclePostSmootherRsdTarget, Reconstruction.MGLevels-1);
	read( recon, "VCyclePostSmootherRelaxOmega", Reconstruction.VCyclePostSmootherRelaxOmega, Reconstruction.MGLevels-1);
	read( recon, "VCyclePostSmootherVerboseP", Reconstruction.VCyclePostSmootherVerboseP, Reconstruction.MGLevels-1);

	read( recon, "VCycleBottomSolverNKrylov", Reconstruction.VCycleBottomSolverNKrylov, Reconstruction.MGLevels-1);
	read( recon, "VCycleBottomSolverMaxIters", Reconstruction.VCycleBottomSolverMaxIters, Reconstruction.MGLevels-1);
	read( recon, "VCycleBottomSolverRsdTarget", Reconstruction.VCycleBottomSolverRsdTarget,Reconstruction.MGLevels-1);
	read( recon, "VCycleBottomSolverVerboseP", Reconstruction.VCycleBottomSolverVerboseP, Reconstruction.MGLevels-1);

	read( recon, "VCycleMaxIters", Reconstruction.VCycleMaxIters, Reconstruction.MGLevels-1);
	read( recon, "VCycleRsdTarget", Reconstruction.VCycleRsdTarget, Reconstruction.MGLevels-1);
	read( recon, "VCycleVerboseP", Reconstruction.VCycleVerboseP, Reconstruction.MGLevels-1);


	read( paramtop, "OuterSolverNKrylov", OuterSolverNKrylov);
	read( paramtop, "OuterSolverRsdTarget", OuterSolverRsdTarget);
	read( paramtop, "OuterSolverMaxIters", OuterSolverMaxIters);
	read( paramtop, "OuterSolverVerboseP", OuterSolverVerboseP);

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


void read(XMLReader& xml, const std::string& path, MGProtoALIPrecParams& p)
{
	MGProtoALIPrecParams ret_val(xml,path);
	p = ret_val;
}

void write(XMLWriter& xml, const std::string& path, const MGProtoALIPrecParams& p)
{
	push(xml, path);
	write(xml, "CloverParams", p.CloverParams);
	write(xml, "AntiPeriodicT", p.AntiPeriodicT);

        push(xml, "MGDeflation");
	write( xml, "MGLevels", p.Deflation.MGLevels);

	write( xml, "Blocking", p.Deflation.Blocking);

	write( xml, "NullVecs", p.Deflation.NullVecs);
	write( xml, "NullSolverMaxIters", p.Deflation.NullSolverMaxIters);
	write( xml, "NullSolverRsdTarget", p.Deflation.NullSolverRsdTarget);
	write( xml, "NullSolverVerboseP", p.Deflation.NullSolverVerboseP);

	write( xml, "EigenSolverMaxRestartSize", p.Deflation.EigenSolverMaxRestartSize);
	write( xml, "EigenSolverMaxRank", p.Deflation.EigenSolverMaxRank);
	write( xml, "EigenSolverRsdTarget", p.Deflation.EigenSolverRsdTarget);
	write( xml, "EigenSolverMaxIters", p.Deflation.EigenSolverMaxIters);
	write( xml, "EigenSolverVerboseP", p.Deflation.EigenSolverVerboseP);

	write( xml, "BottomSolverNKrylov", p.Deflation.BottomSolverNKrylov);
	write( xml, "BottomSolverMaxIters", p.Deflation.BottomSolverMaxIters);
	write( xml, "BottomSolverRsdTarget", p.Deflation.BottomSolverRsdTarget);
	write( xml, "BottomSolverVerboseP", p.Deflation.BottomSolverVerboseP);
	pop(xml);

        push(xml, "Reconstruction");
	write( xml, "ali_distance", p.Reconstruction.ali_distance);
	write( xml, "probing_distance", p.Reconstruction.probing_distance);

	write( xml, "MGLevels", p.Reconstruction.MGLevels);

	write( xml, "Blocking", p.Reconstruction.Blocking);

	write( xml, "NullVecs", p.Reconstruction.NullVecs);
	write( xml, "NullSolverMaxIters", p.Reconstruction.NullSolverMaxIters);
	write( xml, "NullSolverRsdTarget", p.Reconstruction.NullSolverRsdTarget);
	write( xml, "NullSolverVerboseP", p.Reconstruction.NullSolverVerboseP);

	write( xml, "OuterSolverNKrylov", p.Reconstruction.OuterSolverNKrylov);
	write( xml, "OuterSolverRsdTarget", p.Reconstruction.OuterSolverRsdTarget);
	write( xml, "OuterSolverMaxIters", p.Reconstruction.OuterSolverMaxIters);
	write( xml, "OuterSolverVerboseP", p.Reconstruction.OuterSolverVerboseP);

	write( xml, "VCyclePreSmootherMaxIters", p.Reconstruction.VCyclePreSmootherMaxIters);
	write( xml, "VCyclePreSmootherRsdTarget", p.Reconstruction.VCyclePreSmootherRsdTarget);
	write( xml, "VCyclePreSmootherRelaxOmega", p.Reconstruction.VCyclePreSmootherRelaxOmega);
	write( xml, "VCyclePreSmootherVerboseP", p.Reconstruction.VCyclePreSmootherVerboseP);

	write( xml, "VCyclePostSmootherMaxIters", p.Reconstruction.VCyclePostSmootherMaxIters);
	write( xml, "VCyclePostSmootherRsdTarget", p.Reconstruction.VCyclePostSmootherRsdTarget);
	write( xml, "VCyclePostSmootherRelaxOmega", p.Reconstruction.VCyclePostSmootherRelaxOmega);
	write( xml, "VCyclePostSmootherVerboseP", p.Reconstruction.VCyclePostSmootherVerboseP);

	write( xml, "VCycleBottomSolverNKrylov", p.Reconstruction.VCycleBottomSolverNKrylov);
	write( xml, "VCycleBottomSolverMaxIters", p.Reconstruction.VCycleBottomSolverMaxIters);
	write( xml, "VCycleBottomSolverRsdTarget", p.Reconstruction.VCycleBottomSolverRsdTarget);
	write( xml, "VCycleBottomSolverVerboseP", p.Reconstruction.VCycleBottomSolverVerboseP);

	write( xml, "VCycleMaxIters", p.Reconstruction.VCycleMaxIters);
	write( xml, "VCycleRsdTarget", p.Reconstruction.VCycleRsdTarget);
	write( xml, "VCycleVerboseP", p.Reconstruction.VCycleVerboseP);
	pop(xml);

	write( xml, "OuterSolverNKrylov", p.OuterSolverNKrylov);
	write( xml, "OuterSolverRsdTarget", p.OuterSolverRsdTarget);
	write( xml, "OuterSolverMaxIters", p.OuterSolverMaxIters);
	write( xml, "OuterSolverVerboseP", p.OuterSolverVerboseP);

	write(xml, "SubspaceId", p.SubspaceId);
}

}
