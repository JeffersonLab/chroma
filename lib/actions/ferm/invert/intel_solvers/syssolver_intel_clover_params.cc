#include "actions/ferm/invert/intel_solvers/syssolver_intel_clover_params.h"
#include "chromabase.h"
#include "io/xml_group_reader.h"
#include "chroma_config.h"



using namespace QDP;

namespace Chroma {

  SysSolverIntelCloverParams::SysSolverIntelCloverParams(XMLReader& xml, 
							 const std::string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "MaxIter", MaxIter);
    read(paramtop, "RsdTarget", RsdTarget);
    read(paramtop, "CloverParams", CloverParams);
    read(paramtop, "AntiPeriodicT", AntiPeriodicT);
    
    if ( paramtop.count("Verbose") > 0 ) { 
      read(paramtop, "Verbose", VerboseP);
    }
    else { 
      VerboseP = false;
    }
    
    read(paramtop, "NCores", NCores);
    read(paramtop, "ThreadsPerCore", ThreadsPerCore);
    read(paramtop, "By", By);
    read(paramtop, "Bz", Bz);
    read(paramtop, "Sy", Sy);
    read(paramtop, "Sz", Sz);
    
    if( paramtop.count("PadXY") > 0 ) {
      read(paramtop, "PadXY", PadXY);
    }
    else { 
      PadXY = 0;
    }
    
    if( paramtop.count("PadXYZ") > 0 ) {
      read(paramtop, "PadXYZ", PadXYZ);
    }
    else { 
      PadXYZ = 0;
    }
#if 0    
    read(paramtop, "SOALEN", Soalen);
    read(paramtop, "VECLEN", Veclen);
#endif

    if( paramtop.count("MinCt") > 0 ) {
      read(paramtop, "MinCt", MinCt);
    }
    
    if( paramtop.count("RsdToleranceFactor") > 0 ) { 
      read(paramtop, "RsdToleranceFactor", RsdToleranceFactor);
    }
    else { 
      RsdToleranceFactor = Real(10); // Tolerate an order of magnitude difference by default.
    }
    
    if( paramtop.count("Compress") > 0 ) { 
      read(paramtop, "Compress", CompressP);
    }
    else {
      CompressP = false;
    }

    if( paramtop.count("Tune") > 0 ) { 
      read(paramtop, "Tune", TuneP);
    }
    else { 
      TuneP = false;
    }
    
  }
  
  void read(XMLReader& xml, const std::string& path, 
	    SysSolverIntelCloverParams& p)
  {
    SysSolverIntelCloverParams tmp(xml, path);
    p = tmp;
  }
  
  void write(XMLWriter& xml, const std::string& path, 
	     const SysSolverIntelCloverParams& p) {
    push(xml, path);
    write(xml, "MaxIter", p.MaxIter);
    write(xml, "RsdTarget", p.RsdTarget);
    write(xml, "Verbose", p.VerboseP);
    write(xml, "CloverParams", p.CloverParams);
    write(xml, "AntiPeriodicT", p.AntiPeriodicT);
    write(xml, "RsdToleranceFactor", p.RsdToleranceFactor);
    write(xml, "NCores", p.NCores);
    write(xml, "ThreadsPerCore", p.ThreadsPerCore);
    write(xml, "By", p.By);
    write(xml, "Bz", p.Bz);
    write(xml, "Sy", p.Sy);
    write(xml, "Sz", p.Sz);
    write(xml, "PadXY", p.PadXY);
    write(xml, "PadXYZ", p.PadXYZ);
#if 0
    write(xml, "Soalen", p.Soalen);
    write(xml, "Veclen", p.Veclen);
#endif

    write(xml, "MinCt", p.MinCt);
    write(xml, "Compress", p.CompressP);
    write(xml, "Tune", p.TuneP);
    

  }



}
