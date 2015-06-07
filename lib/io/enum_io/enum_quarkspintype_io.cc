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
      success &= theQuarkSpinTypeMap::Instance().registerPair(std::string("FULL"), QUARK_SPIN_TYPE_FULL );
      success &= theQuarkSpinTypeMap::Instance().registerPair(std::string("UPPER"), QUARK_SPIN_TYPE_UPPER);
      success &= theQuarkSpinTypeMap::Instance().registerPair(std::string("LOWER"), QUARK_SPIN_TYPE_LOWER);
      return success;
    }

    bool registered = registerAll();
    const std::string typeIDString = "QuarkSpinType";
  };
  using namespace QuarkSpinTypeEnv;

  //! Read a quark spin type enum
  void read(XMLReader& xml_in,  const std::string& path, QuarkSpinType& t) 
  {
    theQuarkSpinTypeMap::Instance().read(typeIDString, xml_in, path,t);
  }
  
  //! Write a quark spin type enum
  void write(XMLWriter& xml_out, const std::string& path, const QuarkSpinType& t) 
  {
    theQuarkSpinTypeMap::Instance().write(typeIDString, xml_out, path, t);
  }
}
