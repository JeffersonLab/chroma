// -*- C++ -*-
// $Id: qio_read_obj_funcmap.h,v 3.1 2006-09-20 20:28:03 edwards Exp $
/*! \file
 *  \brief Read object function map
 */

#ifndef __qio_read_obj_funcmap_h__
#define __qio_read_obj_funcmap_h__

#include "singleton.h"
#include "funcmap.h"

namespace Chroma
{

  //! Read object function map
  /*! \ingroup inlineio */
  namespace QIOReadObjCallMapEnv
  { 
    struct DumbDisambiguator {};

    //! Read object function map
    /*! \ingroup inlineio */
    typedef SingletonHolder< 
      FunctionMap<DumbDisambiguator,
		  void, 
		  std::string,
		  TYPELIST_3(const string&,
			     const string&, 
			     QDP_serialparallel_t),
		  void (*)(const string& buffer_id,
			   const string& filename, 
			   QDP_serialparallel_t serpar),
		  StringFunctionMapError> >
    TheQIOReadObjFuncMap;

    //! MapObj factory (foundry)
    /*! \ingroup inlineio */
    typedef SingletonHolder< 
      ObjectFactory<MapObject<int,EVPair<LatticeColorVector> >, 
		    std::string,
		    TYPELIST_2(XMLReader&, const std::string&),
		    MapObject<int,EVPair<LatticeColorVector> >* (*)(XMLReader&,
								    const std::string&), 
		    StringFactoryError> >
    TheMapObjIntKeyColorEigenVecFactory;

    bool registerAll();
  }

} // end namespace Chroma


#endif
