// -*- C++ -*-
// $Id: szin_read_obj_funcmap.h,v 3.1 2006-09-20 20:28:03 edwards Exp $
/*! \file
 *  \brief Read object function std::map
 */

#ifndef __szin_read_obj_funcmap_h__
#define __szin_read_obj_funcmap_h__

#include "singleton.h"
#include "funcmap.h"
#include "chromabase.h"

namespace Chroma
{

  //! Read object function std::map
  /*! \ingroup inlineio */
  namespace SZINReadObjCallMapEnv
  { 
    struct DumbDisambiguator {};

    //! Read object function std::map
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
    TheSZINReadObjFuncMap;

    bool registerAll();
  }

} // end namespace Chroma


#endif
