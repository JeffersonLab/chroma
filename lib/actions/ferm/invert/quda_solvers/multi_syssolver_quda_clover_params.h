#ifndef __MULTI_SYSSOLVER_QUDA_CLOVER_PARAMS_H__
#define __MULTI_SYSSOLVER_QUDA_CLOVER_PARAMS_H__

#include "chromabase.h"
#include "io/xml_group_reader.h"

#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/invert/quda_solvers/enum_quda_io.h"
#include "actions/ferm/invert/quda_solvers/quda_gcr_params.h"
#include <string>
#include "handle.h"

namespace Chroma 
{
  struct MultiSysSolverQUDACloverParams { 
    MultiSysSolverQUDACloverParams(XMLReader& xml, const std::string& path);
    MultiSysSolverQUDACloverParams() {
      solverType=CG;
      cudaPrecision=DEFAULT;
      cudaReconstruct=RECONS_12;
      cudaSloppyPrecision=DEFAULT;
      cudaSloppyReconstruct=RECONS_12;
      cudaRefinementPrecision=DEFAULT;
      cudaRefinementReconstruct=RECONS_12;
      asymmetricP = false; //< Use asymmetric version of the linear operator
      axialGaugeP = false; //< Fix Axial Gauge?
      SilentFailP = false; //< If set to true ignore lack of convergence. Default is 'loud' 
      RsdToleranceFactor = Real(10); //< Tolerate if the solution achived is better (less) than rsdToleranceFactor*RsdTarget
      tuneDslashP = false ; //< v0.3 autotune feature
      verboseP = false;
      innerParamsP = false;
      checkShiftsP = true;
      Pipeline=0;
    };
    MultiSysSolverQUDACloverParams( const MultiSysSolverQUDACloverParams& p) {
      CloverParams = p.CloverParams;
      AntiPeriodicT = p.AntiPeriodicT;
      MaxIter = p.MaxIter;
      RsdTarget = p.RsdTarget;
      Delta = p.Delta;
      solverType = p.solverType;
      verboseP = p.verboseP;
      asymmetricP = p.asymmetricP;
      cudaPrecision = p.cudaPrecision;
      cudaReconstruct = p.cudaReconstruct;
      cudaSloppyPrecision = p.cudaSloppyPrecision;
      cudaSloppyReconstruct = p.cudaSloppyReconstruct;
      cudaRefinementPrecision = p.cudaRefinementPrecision;
      cudaRefinementReconstruct = p.cudaRefinementReconstruct;
      axialGaugeP = p.axialGaugeP;
      SilentFailP = p.SilentFailP;
      RsdToleranceFactor = p.RsdToleranceFactor;
      tuneDslashP = p.tuneDslashP;
      innerParamsP = p.innerParamsP;
      innerParams = p.innerParams;
      checkShiftsP = p.checkShiftsP;
      Pipeline = p.Pipeline;
    }

   
    CloverFermActParams CloverParams;
    bool AntiPeriodicT;
    int MaxIter;
    multi1d<Real> RsdTarget;
    Real Delta;
    QudaSolverType solverType;
    bool verboseP;
    bool asymmetricP;
    QudaPrecisionType cudaPrecision;
    QudaReconsType cudaReconstruct;
    QudaPrecisionType cudaSloppyPrecision;
    QudaReconsType cudaSloppyReconstruct;
    QudaPrecisionType cudaRefinementPrecision;
    QudaReconsType cudaRefinementReconstruct;
    bool axialGaugeP;
    bool SilentFailP;
    Real RsdToleranceFactor;
    bool tuneDslashP;
    bool innerParamsP;
    bool checkShiftsP;
    int Pipeline;

    // GCR Specific params
    // Params for the preconditioner
    Handle<GCRInnerSolverParams> innerParams;


  };

  void read(XMLReader& xml, const std::string& path, MultiSysSolverQUDACloverParams& p);

  void write(XMLWriter& xml, const std::string& path, 
	     const MultiSysSolverQUDACloverParams& param);



}

#endif


