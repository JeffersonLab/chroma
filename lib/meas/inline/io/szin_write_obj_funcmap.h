// -*- C++ -*-
/*! \file
 *  \brief Write object function std::map
 */

#ifndef __szin_write_obj_funcmap_h__
#define __szin_write_obj_funcmap_h__

#include "singleton.h"
#include "funcmap.h"
#include "chromabase.h"

namespace Chroma
{

  //! Write object function std::map
  /*! \ingroup inlineio */
  namespace SZINWriteObjCallMapEnv
  { 
    struct DumbDisambiguator {};

    //! Write object function std::map
    /*! \ingroup inlineio */
    typedef Chroma::SingletonHolder< 
      FunctionMap<DumbDisambiguator,
		  void,
		  std::string,
		  TYPELIST_5(const std::string&,
			     const std::string&, 
			     int, int, int),
		  void (*)(const std::string& buffer_id,
			   const std::string& filename, 
			   int j_decay, int t_start, int t_end),
		  StringFunctionMapError> >
    TheSZINWriteObjFuncMap;

    bool registerAll();
  }

} // end namespace Chroma


#endif
