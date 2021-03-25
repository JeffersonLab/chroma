// -*- C++ -*-
/*! \file
 *  \brief Write object function std::map
 */

#ifndef __qio_write_obj_funcmap_h__
#define __qio_write_obj_funcmap_h__

#include "singleton.h"
#include "funcmap.h"
#include "chromabase.h"

namespace Chroma
{

  //! Write object function std::map
  /*! \ingroup inlineio */
  namespace QIOWriteObjCallMapEnv
  { 
    struct DumbDisambiguator {};

    //! Write object function std::map
    /*! \ingroup inlineio */
    typedef Chroma::SingletonHolder< 
      FunctionMap<DumbDisambiguator,
		  void,
		  std::string,
		  TYPELIST_4(const std::string&,
			     const std::string&, 
			     QDP_volfmt_t, QDP_serialparallel_t),
		  void (*)(const std::string& buffer_id,
			   const std::string& filename, 
			   QDP_volfmt_t volfmt, QDP_serialparallel_t serpar),
		  StringFunctionMapError> >
    TheQIOWriteObjFuncMap;

    bool registerAll();
  }

} // end namespace Chroma


#endif
