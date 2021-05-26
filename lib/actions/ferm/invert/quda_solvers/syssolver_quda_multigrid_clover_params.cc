#include "actions/ferm/invert/quda_solvers/syssolver_quda_multigrid_clover_params.h"
#include "chromabase.h"
#include "io/xml_group_reader.h"
#include "chroma_config.h"



using namespace QDP;

namespace Chroma {

  template <typename T>
  void read(XMLReader& xml, const std::string& path, T& result, T default_value)
  {
    if (xml.count(path) > 0)
    {
      read(xml, path, result);
    }
    else
    {
      result = default_value;
    }
  }

  SysSolverQUDAMULTIGRIDCloverParams::SysSolverQUDAMULTIGRIDCloverParams(XMLReader& xml,
									 const std::string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "MaxIter", MaxIter);
    read(paramtop, "RsdTarget", RsdTarget);
    read(paramtop, "CloverParams", CloverParams);
    read(paramtop, "AntiPeriodicT", AntiPeriodicT);

    read(paramtop, "Delta", Delta);

    read(paramtop, "SolverType", solverType);

    read(paramtop, "Verbose", verboseP, false);
    read(paramtop, "AsymmetricLinop", asymmetricP, false);
    read(paramtop, "CudaPrecision", cudaPrecision, DEFAULT);
    read(paramtop, "CudaSloppyPrecision", cudaSloppyPrecision, DEFAULT);
    read(paramtop, "CudaReconstruct", cudaReconstruct, RECONS_12);
    read(paramtop, "CudaSloppyReconstruct", cudaSloppyReconstruct, RECONS_12);
    read(paramtop, "AxialGaugeFix", axialGaugeP, false);
    read(paramtop, "SilentFail", SilentFailP, false);
    read(paramtop, "RsdToleranceFactor", RsdToleranceFactor, Real(10));
    read(paramtop, "AutotuneDslash", tuneDslashP, false);

    read(paramtop, "SubspaceID", SaveSubspaceID);

    read(paramtop, "ThresholdCount", ThresholdCount, 2 * MaxIter + 1);

    read(paramtop, "Pipeline", Pipeline, 1);

    if (paramtop.count("MULTIGRIDParams") > 0)
    {
      MULTIGRIDParams = new MULTIGRIDSolverParams(paramtop, "./MULTIGRIDParams");
      MULTIGRIDParamsP = true;
    }
    else
    {
      MULTIGRIDParamsP = false;
    }

    if (paramtop.count("BackupSolverParam") > 0)
    {
      // If user specified a backup solver, let's read it
      backup_invP = true;
      backup_inv_param = readXMLGroup(paramtop, "./BackupSolverParam", "invType");
    }
    else
    {
      // No backup
      backup_invP = false;
      backup_inv_param.xml = "";
    }

    read(paramtop, "DumpOnFail", dump_on_failP, false);

    read(paramtop, "SolutionCheckP", SolutionCheckP, true);
  }

  void read(XMLReader& xml, const std::string& path, SysSolverQUDAMULTIGRIDCloverParams& p)
  {
    SysSolverQUDAMULTIGRIDCloverParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, const std::string& path, const SysSolverQUDAMULTIGRIDCloverParams& p)
  {
    push(xml, path);
    write(xml, "MaxIter", p.MaxIter);
    write(xml, "RsdTarget", p.RsdTarget);
    write(xml, "CloverParams", p.CloverParams);
    write(xml, "AntiPeriodicT", p.AntiPeriodicT);
    write(xml, "Delta", p.Delta);
    write(xml, "SolverType", p.solverType);
    write(xml, "Verbose", p.verboseP);
    write(xml, "AsymmetricLinop", p.asymmetricP);
    write(xml, "CudaPrecision", p.cudaPrecision);
    write(xml, "CudaReconstruct", p.cudaReconstruct);
    write(xml, "CudaSloppyPrecision", p.cudaSloppyPrecision);
    write(xml, "CudaSloppyReconstruct", p.cudaSloppyReconstruct);
    write(xml, "AxialGaugeFix", p.axialGaugeP);
    write(xml, "SilentFail", p.SilentFailP);
    write(xml, "RsdToleranceFactor", p.RsdToleranceFactor);

    write(xml, "AutotuneDslash", p.tuneDslashP);

    //Write the MG persistence params.
    write(xml, "SubspaceID", p.SaveSubspaceID);
    write(xml, "ThresholdCount", p.ThresholdCount);
    write(xml, "Pipeline", p.Pipeline);
 
    if( p.MULTIGRIDParamsP ) { 
      write(xml, "MULTIGRIDParams", *(p.MULTIGRIDParams));
    }

    write(xml, "DumpOnFail", p.dump_on_failP);
    write(xml, "SolutionCheckP", p.SolutionCheckP);

    if( p.backup_invP ) { 
      // Need to dump out the XML for the back up solver here...
      // Turn XML into an std::istringstream
      std::istringstream is( p.backup_inv_param.xml);
      XMLReader tmp_reader(is);
      // I wonder if this method works....
      write(xml, "BackupSolverParam", tmp_reader);
    }
    pop(xml);
  }

  MugiqMGDeflationCloverParams::MugiqMGDeflationCloverParams(XMLReader& xml,
							     const std::string& path)
    : SysSolverQUDAMULTIGRIDCloverParams(xml, path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "EigenSolverMaxRestartSize", EigenSolverMaxRestartSize);
    read(paramtop, "EigenSolverMaxRank", EigenSolverMaxRank);
  }

  void read(XMLReader& xml, const std::string& path, MugiqMGDeflationCloverParams& p)
  {
    MugiqMGDeflationCloverParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, const std::string& path, const MugiqMGDeflationCloverParams& p)
  {
    push(xml, path);
    write(xml, "MaxIter", p.MaxIter);
    write(xml, "RsdTarget", p.RsdTarget);
    write(xml, "CloverParams", p.CloverParams);
    write(xml, "AntiPeriodicT", p.AntiPeriodicT);
    write(xml, "Delta", p.Delta);
    write(xml, "SolverType", p.solverType);
    write(xml, "Verbose", p.verboseP);

    write(xml, "EigenSolverMaxRestartSize", p.EigenSolverMaxRestartSize);
    write(xml, "EigenSolverMaxRank", p.EigenSolverMaxRank);

    write(xml, "AsymmetricLinop", p.asymmetricP);
    write(xml, "CudaPrecision", p.cudaPrecision);
    write(xml, "CudaReconstruct", p.cudaReconstruct);
    write(xml, "CudaSloppyPrecision", p.cudaSloppyPrecision);
    write(xml, "CudaSloppyReconstruct", p.cudaSloppyReconstruct);
    write(xml, "AxialGaugeFix", p.axialGaugeP);
    write(xml, "SilentFail", p.SilentFailP);
    write(xml, "RsdToleranceFactor", p.RsdToleranceFactor);

    write(xml, "AutotuneDslash", p.tuneDslashP);

    //Write the MG persistence params.
    write(xml, "SubspaceID", p.SaveSubspaceID);
    write(xml, "ThresholdCount", p.ThresholdCount);
    write(xml, "Pipeline", p.Pipeline);
 
    if( p.MULTIGRIDParamsP ) { 
      write(xml, "MULTIGRIDParams", *(p.MULTIGRIDParams));
    }

    write(xml, "DumpOnFail", p.dump_on_failP);
    write(xml, "SolutionCheckP", p.SolutionCheckP);

    if( p.backup_invP ) { 
      // Need to dump out the XML for the back up solver here...
      // Turn XML into an std::istringstream
      std::istringstream is( p.backup_inv_param.xml);
      XMLReader tmp_reader(is);
      // I wonder if this method works....
      write(xml, "BackupSolverParam", tmp_reader);
    }
    pop(xml);
  }
}
