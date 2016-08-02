#include "actions/ferm/invert/qphix/syssolver_qphix_clover_params.h"
#include "chromabase.h"
#include "io/xml_group_reader.h"
#include "chroma_config.h"



using namespace QDP;

namespace Chroma {

  namespace QPhiXSolverTypeEnv { 

    bool registerAll(void) 
    {
      bool success; 
      success = theQPhiXSolverTypeMap::Instance().registerPair(std::string("CG"),
							   CG);

      success &=theQPhiXSolverTypeMap::Instance().registerPair(std::string("BICGSTAB" ), 
							       BICGSTAB);

      return success;
    }

    const std::string typeIDString = "QPhiXSolverType";
    bool registered = registerAll();
  };
  
  using namespace QPhiXSolverTypeEnv;
  //! Read an WaveType enum
  void read(XMLReader& xml_in,  const std::string& path, QPhiXSolverType& t) {
    theQPhiXSolverTypeMap::Instance().read(typeIDString, xml_in, path,t);
  }
  
  //! Write an WaveType enum
  void write(XMLWriter& xml_out, const std::string& path, const QPhiXSolverType& t) {
    theQPhiXSolverTypeMap::Instance().write(typeIDString, xml_out, path, t);
  }



  SysSolverQPhiXCloverParams::SysSolverQPhiXCloverParams(XMLReader& xml, 
							 const std::string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "MaxIter", MaxIter);
    read(paramtop, "RsdTarget", RsdTarget);
    read(paramtop, "CloverParams", CloverParams);
    read(paramtop, "SolverType", SolverType);
    read(paramtop, "AntiPeriodicT", AntiPeriodicT);
    
    if ( paramtop.count("Verbose") > 0 ) { 
      read(paramtop, "Verbose", VerboseP);
    }
    else { 
      VerboseP = false;
    }

    if (paramtop.count("Delta") > 0 ) { 
      read(paramtop,"Delta", Delta);
    }
    else { 
      Delta = Real(-1);
    }
    
    if( paramtop.count("RsdToleranceFactor") > 0 ) { 
      read(paramtop, "RsdToleranceFactor", RsdToleranceFactor);
    }
    else { 
      RsdToleranceFactor = Real(10); // Tolerate an order of magnitude difference by default.
    }

    if( paramtop.count("Tune") > 0 ) { 
      read(paramtop, "Tune", TuneP);
    }
    else { 
      TuneP = false;
    }
    
  }
  
  void read(XMLReader& xml, const std::string& path, 
	    SysSolverQPhiXCloverParams& p)
  {
    SysSolverQPhiXCloverParams tmp(xml, path);
    p = tmp;
  }
  
  void write(XMLWriter& xml, const std::string& path, 
	     const SysSolverQPhiXCloverParams& p) {
    push(xml, path);
    write(xml, "MaxIter", p.MaxIter);

    // Hack. Write delta only if it is greater than 0. 
    // -ve sign indicates user did not supply
    if( toBool(p.Delta > 0 ) ) {
      write(xml, "Delta", p.Delta);
    }

    write(xml, "RsdTarget", p.RsdTarget);
    write(xml, "Verbose", p.VerboseP);
    write(xml, "CloverParams", p.CloverParams);
    write(xml, "SolverType", p.SolverType);
    write(xml, "AntiPeriodicT", p.AntiPeriodicT);
    write(xml, "RsdToleranceFactor", p.RsdToleranceFactor);
    write(xml, "Tune", p.TuneP);
    

  }



}
