// $Id: enum_quarkspintype_io.cc,v 3.0 2006-04-03 04:58:56 edwards Exp $
/*! \file
 * \brief Enum for what spin components of a quark prop to compute
 */

#include "enum_quarkspintype_io.h"

namespace Chroma 
{ 
  namespace QuarkSpinTypeEnv 
  { 
    bool registerAll(void) 
    {
      bool success = true; 
      success &= theQuarkSpinTypeMap::Instance().registerPair(string("FULL"), QUARK_SPIN_TYPE_FULL );
      success &= theQuarkSpinTypeMap::Instance().registerPair(string("UPPER"), QUARK_SPIN_TYPE_UPPER);
      success &= theQuarkSpinTypeMap::Instance().registerPair(string("LOWER"), QUARK_SPIN_TYPE_LOWER);
      return success;
    }

    bool registered = registerAll();
    const string typeIDString = "QuarkSpinType";
  };
  using namespace QuarkSpinTypeEnv;

  //! Read a quark spin type enum
  void read(XMLReader& xml_in,  const string& path, QuarkSpinType& t) 
  {
    theQuarkSpinTypeMap::Instance().read(typeIDString, xml_in, path,t);
  }
  
  //! Write a quark spin type enum
  void write(XMLWriter& xml_out, const string& path, const QuarkSpinType& t) 
  {
    theQuarkSpinTypeMap::Instance().write(typeIDString, xml_out, path, t);
  }
}
