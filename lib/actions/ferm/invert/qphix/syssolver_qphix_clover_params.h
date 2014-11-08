#ifndef __SYSSOLVER_QPHIX_CLOVER_PARAMS_H__
#define __SYSSOLVER_QPHIX_CLOVER_PARAMS_H__

#include "chromabase.h"
#include "io/xml_group_reader.h"
#include "io/enum_io/enum_type_map.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include <string>
#include "handle.h"

namespace Chroma 
{

  enum QPhiXSolverType 
  { 
    CG, 
    BICGSTAB 
  };

  namespace QPhiXSolverTypeEnv { 
    extern const std::string typeIDString;
    extern bool registered;
    bool registerAll(void);
  }

  // A singleton to hold the typemap
  typedef SingletonHolder<EnumTypeMap<QPhiXSolverType> > theQPhiXSolverTypeMap;

  //! Read an WaveStateType enum
  void read(XMLReader& r, const std::string& path, QPhiXSolverType& t);

  //! Write an WaveStateType enum
  void write(XMLWriter& w, const std::string& path, const QPhiXSolverType& t);



  struct SysSolverQPhiXCloverParams { 
    SysSolverQPhiXCloverParams(XMLReader& xml, const std::string& path);
    SysSolverQPhiXCloverParams() {
       RsdToleranceFactor = Real(10); //< Tolerate if the solution achived is better (less) than rsdToleranceFactor*RsdTarget
      TuneP = false ; //< v0.3 autotune feature
      VerboseP = false;
      MinCt = 1;
      Delta = Real(-1);
      SolverType = CG;
    };

    SysSolverQPhiXCloverParams( const SysSolverQPhiXCloverParams& p) {
      CloverParams = p.CloverParams;
      AntiPeriodicT = p.AntiPeriodicT;
      MaxIter = p.MaxIter;
      Delta = p.Delta;
      RsdTarget = p.RsdTarget;
      VerboseP = p.VerboseP;
      NCores = p.NCores;
      By=p.By;
      Bz=p.Bz;
      Sy=p.Sy;
      Sz=p.Sz;
      PadXY = p.PadXY;
      PadXYZ = p.PadXYZ;
      RsdToleranceFactor = p.RsdToleranceFactor;
      TuneP = p.TuneP;
      MinCt = p.MinCt;
      SolverType = p.SolverType;
    }

   
    CloverFermActParams CloverParams;
    QPhiXSolverType SolverType;
    bool AntiPeriodicT;
    int MaxIter;
    Real RsdTarget;
    Real Delta;
    bool VerboseP;
    int NCores;
    int By;
    int Bz;
    int Sy;
    int Sz;
    int PadXY;
    int PadXYZ;
    int MinCt;
    Real RsdToleranceFactor;
    bool TuneP;
  };

  void read(XMLReader& xml, const std::string& path, SysSolverQPhiXCloverParams& p);

  void write(XMLWriter& xml, const std::string& path, 
	     const SysSolverQPhiXCloverParams& param);



}

#endif


