// -*- C++ -*-
/*! 
 *  \brief Factory for creating various kinds of map objects
 */

#ifndef __map_obj_factory_w_h__
#define __map_obj_factory_w_h__

#include "chromabase.h"
#include "handle.h"
#include "state.h"
#include "singleton.h"
#include "typelist.h"
#include "objfactory.h"
#include "util/ferm/map_obj.h"
#include "util/ferm/key_block_prop.h"
#include "util/ferm/key_grid_prop.h"
#include "util/ferm/key_prop_colorvec.h"

namespace Chroma
{
  namespace { 
  //! MapObj factory (foundry)
  /*! @ingroup invert */
    typedef SingletonHolder< 
      ObjectFactory<MapObject<KeyBlockProp_t,LatticeFermion>, 
		    std::string,
		    TYPELIST_2(XMLReader&, const std::string&),
		    MapObject<KeyBlockProp_t,LatticeFermion>* (*)(XMLReader&,
								  const std::string&), 
		    StringFactoryError> >
    TheMapObjKeyBlockPropFactory;

    //! MapObj factory (foundry)
    /*! @ingroup invert */
    typedef SingletonHolder< 
      ObjectFactory<MapObject<KeyGridProp_t,LatticeFermion>, 
		    std::string,
		    TYPELIST_2(XMLReader&, const std::string&),
		    MapObject<KeyGridProp_t,LatticeFermion>* (*)(XMLReader&,
								 const std::string&), 
		    StringFactoryError> >
    TheMapObjKeyGridPropFactory;

    //! MapObj factory (foundry)
    /*! @ingroup invert */
    typedef SingletonHolder< 
      ObjectFactory<MapObject<KeyPropColorVec_t,LatticeFermion>, 
		    std::string,
		    TYPELIST_2(XMLReader&, const std::string&),
		    MapObject<KeyPropColorVec_t,LatticeFermion>* (*)(XMLReader&,
								     const std::string&), 
		    StringFactoryError> >
    TheMapObjKeyPropColorVecFactory;


  }

}


#endif
