// -*- C++ -*-
// $Id: szin_write_obj_funcmap.h,v 3.1 2006-09-20 20:28:03 edwards Exp $
/*! \file
 *  \brief Write object function map
 */

#ifndef __szin_write_obj_funcmap_h__
#define __szin_write_obj_funcmap_h__

#include "singleton.h"
#include "funcmap.h"
#include "chromabase.h"

namespace Chroma
{

  //! Write object function map
  /*! \ingroup inlineio */
  namespace SZINWriteObjCallMapEnv
  { 
    struct DumbDisambiguator {};

    //! Write object function map
    /*! \ingroup inlineio */
    typedef SingletonHolder< 
      FunctionMap<DumbDisambiguator,
		  void,
		  std::string,
		  TYPELIST_5(const string&,
			     const string&, 
			     int, int, int),
		  void (*)(const string& buffer_id,
			   const string& filename, 
			   int j_decay, int t_start, int t_end),
		  StringFunctionMapError> >
    TheSZINWriteObjFuncMap;

    bool registerAll();
  }

} // end namespace Chroma


#endif
