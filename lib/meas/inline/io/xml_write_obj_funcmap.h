// -*- C++ -*-
// $Id: xml_write_obj_funcmap.h,v 3.1 2006-09-20 20:28:03 edwards Exp $
/*! \file
 *  \brief Write object function std::map
 */

#ifndef __xml_write_obj_funcmap_h__
#define __xml_write_obj_funcmap_h__

#include "singleton.h"
#include "funcmap.h"
#include "chromabase.h"

namespace Chroma
{

  //! Write object function std::map
  /*! \ingroup inlineio */
  namespace XMLWriteObjCallMapEnv
  { 
    struct DumbDisambiguator {};

    //! Write object function std::map
    /*! \ingroup inlineio */
    typedef SingletonHolder< 
      FunctionMap<DumbDisambiguator,
		  void,
		  std::string,
		  TYPELIST_2(const std::string&,
			     const std::string&),
		  void (*)(const std::string& buffer_id,
			   const std::string& filename),
		  StringFunctionMapError> >
    TheXMLWriteObjFuncMap;

    bool registerAll();
  }

} // end namespace Chroma


#endif
