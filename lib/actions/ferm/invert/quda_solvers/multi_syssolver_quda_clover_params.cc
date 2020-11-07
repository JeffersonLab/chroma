#include "actions/ferm/invert/quda_solvers/multi_syssolver_quda_clover_params.h"
#include "chromabase.h"
#include "io/xml_group_reader.h"
#include "chroma_config.h"



using namespace QDP;

namespace Chroma {

  MultiSysSolverQUDACloverParams::MultiSysSolverQUDACloverParams(XMLReader& xml, 
						       const std::string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "MaxIter", MaxIter);
    read(paramtop, "RsdTarget", RsdTarget);
    read(paramtop, "CloverParams", CloverParams);
    read(paramtop, "AntiPeriodicT", AntiPeriodicT);

    read(paramtop, "Delta", Delta);
  

    read(paramtop, "SolverType", solverType);

    if ( paramtop.count("Verbose") > 0 ) { 
      read(paramtop, "Verbose", verboseP);
    }
    else { 
      verboseP = false;
    }
    if ( paramtop.count("AsymmetricLinop") > 0 ) { 
      read(paramtop, "AsymmetricLinop", asymmetricP);
    }
    else { 
      asymmetricP = false; // Symmetric is default 
    }

    if( paramtop.count("CudaPrecision") > 0 ) {
      read(paramtop, "CudaPrecision", cudaPrecision);
    }
    else { 
      cudaPrecision = DEFAULT;
    }

    if( paramtop.count("CudaSloppyPrecision") > 0 ) {
      read(paramtop, "CudaSloppyPrecision", cudaSloppyPrecision);
    }
    else { 
      cudaSloppyPrecision = DEFAULT;
    }

    if( paramtop.count("CudaRefinementPrecision") > 0 ) {
    	read(paramtop, "CudaRefinementPrecision", cudaRefinementPrecision);
    }
    else {
    	cudaRefinementPrecision = DEFAULT;
    }

    if( paramtop.count("CudaReconstruct") > 0 ) {
      read(paramtop, "CudaReconstruct", cudaReconstruct);
    }
    else { 
      cudaReconstruct = RECONS_12;
    }

    if( paramtop.count("CudaSloppyReconstruct") > 0 ) {
      read(paramtop, "CudaSloppyReconstruct", cudaSloppyReconstruct);
    }
    else { 
      cudaSloppyReconstruct = RECONS_12;
    }

    if( paramtop.count("CudaRefinementReconstruct") > 0 ) {
    	read(paramtop, "CudaRefinementReconstruct", cudaRefinementReconstruct);
    }
    else {
    	cudaRefinementReconstruct = RECONS_12;
    }

    if( paramtop.count("AxialGaugeFix") > 0 ) {
      read(paramtop, "AxialGaugeFix", axialGaugeP);
    }
    else { 
      axialGaugeP = false;
    }

    if( paramtop.count("SilentFail") > 0) {
      read(paramtop, "SilentFail", SilentFailP);
    }
    else { 
      SilentFailP = false;
    }

    if( paramtop.count("RsdToleranceFactor") > 0 ) { 
       read(paramtop, "RsdToleranceFactor", RsdToleranceFactor);
    }
    else { 
       RsdToleranceFactor = Real(10); // Tolerate an order of magnitude difference by default.
    }

    if( paramtop.count("AutotuneDslash") > 0 ) { 
      read(paramtop, "AutotuneDslash", tuneDslashP);
    }
    else { 
      tuneDslashP = false;
    }
    QDPIO::cout << "tuneDslasP = " << tuneDslashP << std::endl;


    if( paramtop.count("GCRInnerParams") > 0 ) {
      innerParams = new GCRInnerSolverParams(paramtop, "./GCRInnerParams");
      innerParamsP = true;
    }
    else { 
      innerParamsP = false;
    }

    if ( paramtop.count("CheckShifts") == 1 ) {
      read(paramtop, "CheckShifts", checkShiftsP );
    }
    // Backward compatibility
    if( paramtop.count("check_shifts") == 1) {
      read(paramtop, "check_shifts", checkShiftsP) ;
    }

    if( paramtop.count("Pipeline") == 1 ) { 
      read(paramtop, "Pipeline", Pipeline);
    }
    else {
      Pipeline = 1;
    }

  }

  void read(XMLReader& xml, const std::string& path, 
	    MultiSysSolverQUDACloverParams& p)
  {
    MultiSysSolverQUDACloverParams tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, const std::string& path, 
	     const MultiSysSolverQUDACloverParams& p) {
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
    write(xml, "CheckShifts", p.checkShiftsP);
    write(xml, "AutotuneDslash", p.tuneDslashP);
    write(xml, "Pipeline", p.Pipeline);

    if( p.innerParamsP ) { 
      write(xml, "GCRInnerParams", *(p.innerParams));
    }
    pop(xml);

  }



}
