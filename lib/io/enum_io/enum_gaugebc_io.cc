#include "enum_gaugebc_io.h"

#include <string>
using namespace std;
using namespace Chroma;

namespace Chroma { 


  /*****************  SOURCES *****************************/
  namespace GaugeBCTypeEnv { 

    bool registerAll(void) 
    {
      bool success; 
      success = theGaugeBCTypeMap::Instance().registerPair(string( "ALL_PERIODIC" ),
							   GAUGEBC_ALL_PERIODIC);

      success &=theGaugeBCTypeMap::Instance().registerPair(string("SCHROEDINRGER_1LINK" ), 
							   GAUGEBC_SCHROEDINGER_1LINK );
      
      success &=theGaugeBCTypeMap::Instance().registerPair(string("SCHROEDINGER_2LINK" ), 
							   GAUGEBC_SCHROEDINGER_2LINK);

      success &=theGaugeBCTypeMap::Instance().registerPair(string( "SIMPLE" ), 
							   GAUGEBC_SIMPLE);

      return success;
    }
    const string typeIDString = "GaugeBCType";
    bool registered = registerAll();
  };
  
  //! Read an GaugeBCType enum
  void read(XMLReader& xml_in,  const string& path, GaugeBCType& t) {
    theGaugeBCTypeMap::Instance().read(GaugeBCTypeEnv::typeIDString,
				       xml_in, path,t);
  }
  
  //! Write an GaugeBCType enum
  void write(XMLWriter& xml_out, const string& path, const GaugeBCType& t) {
    theGaugeBCTypeMap::Instance().write(GaugeBCTypeEnv::typeIDString,
					xml_out, path, t);
  }

  /*****************  SINKS *****************************/
  namespace SchrFunTypeEnv { 

    bool registerAll(void) 
    {
      bool success; 
      success = theSchrFunTypeMap::Instance().registerPair(string(  "NONE" ),
							SF_NONE);

      success &=theSchrFunTypeMap::Instance().registerPair(string("TRIVIAL" ), 
							   SF_TRIVIAL);
      
      success &=theSchrFunTypeMap::Instance().registerPair(string("NONPERTURBATIVE" ), 
							   SF_NONPERT);

      success &=theSchrFunTypeMap::Instance().registerPair(string( "COUPLING" ), 
							   SF_COUPLING);
      
      success &=theSchrFunTypeMap::Instance().registerPair(string("CHROMOMAGNETIC" ), 
							   SF_CHROMOMAG);

      success &=theSchrFunTypeMap::Instance().registerPair(string( "DIRICHLET" ), 
							   SF_DIRICHLET);
      
      return success;
    }

    const string typeIDString = "SchrFunType";
    bool registered = registerAll();
  };

  //! Read an SchrFunType enum
  void read(XMLReader& xml_in,  const string& path, SchrFunType& t) {
    theSchrFunTypeMap::Instance().read(SchrFunTypeEnv::typeIDString,
				       xml_in, path,t);
  }
  
  //! Write an SchrFunType enum
  void write(XMLWriter& xml_out, const string& path, const SchrFunType& t) {
    theSchrFunTypeMap::Instance().write(SchrFunTypeEnv::typeIDString,
					xml_out, path, t);
  }


};
