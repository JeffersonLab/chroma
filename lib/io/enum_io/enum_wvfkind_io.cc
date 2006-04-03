// -*- C++ -*-
// $Id: enum_wvfkind_io.cc,v 3.0 2006-04-03 04:58:57 edwards Exp $
/*! \file
 * \brief Wavekind enum
 */
#include "enum_wvfkind_io.h"

#include <string>

namespace Chroma { 

  namespace WvfKindEnv { 

    bool registerAll(void) 
    {
      bool success; 
      success = theWvfKindMap::Instance().registerPair(string("GAUSSIAN"),
						       WVF_KIND_GAUSSIAN   );

      success &=theWvfKindMap::Instance().registerPair(string("EXPONENTIAL"), 
						       WVF_KIND_EXPONENTIAL );

      success &=theWvfKindMap::Instance().registerPair(string("GAUGE_INV_GAUSSIAN"), 
						       WVF_KIND_GAUGE_INV_GAUSSIAN);

      success &=theWvfKindMap::Instance().registerPair(string( "WUPPERTAL"), 
						       WVF_KIND_WUPPERTAL);

      success &=theWvfKindMap::Instance().registerPair(string("JACOBI"), 
						       WVF_KIND_JACOBI);

      
      return success;
    }

    const string typeIDString = "WvfKind";
    bool registered = registerAll();
  };

  using namespace WvfKindEnv;

  //! Read an Smearing Wavefunction kind enum
  void read(XMLReader& xml_in,  const string& path, WvfKind& t) {
    theWvfKindMap::Instance().read(typeIDString, xml_in, path,t);
  }
  
  //! Write an Smearing Wavefunction kind enum
  void write(XMLWriter& xml_out, const string& path, const WvfKind& t) {
    theWvfKindMap::Instance().write(typeIDString, xml_out, path, t);
  }
};
