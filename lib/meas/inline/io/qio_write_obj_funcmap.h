// -*- C++ -*-
// $Id: qio_write_obj_funcmap.h,v 2.1 2006-03-24 22:16:40 edwards Exp $
/*! \file
 *  \brief Write object function map
 */

#ifndef __qio_write_obj_funcmap_h__
#define __qio_write_obj_funcmap_h__

#include "singleton.h"
#include "funcmap.h"
#include "chromabase.h"

namespace Chroma
{

  //! Write object function map
  /*! \ingroup inlineio */
  namespace QIOWriteObjCallMapEnv
  { 
    struct DumbDisambiguator {};

    //! Write object function map
    /*! \ingroup inlineio */
    typedef SingletonHolder< 
      FunctionMap<DumbDisambiguator,
		  void,
		  std::string,
		  TYPELIST_4(const string&,
			     const string&, 
			     QDP_volfmt_t, QDP_serialparallel_t),
		  void (*)(const string& buffer_id,
			   const string& filename, 
			   QDP_volfmt_t volfmt, QDP_serialparallel_t serpar),
		  StringFunctionMapError> >
    TheQIOWriteObjFuncMap;

    extern bool registered;   // forward decl
  };

} // end namespace Chroma


#endif
