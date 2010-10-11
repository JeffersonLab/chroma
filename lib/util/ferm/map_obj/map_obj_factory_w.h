// -*- C++ -*-
/*! 
 *  \brief Factory for creating various kinds of map objects
 */

#ifndef __map_obj_factory_w_h__
#define __map_obj_factory_w_h__

#include "chromabase.h"
#include "handle.h"
#include "singleton.h"
#include "typelist.h"
#include "objfactory.h"
#include "qdp_map_obj.h"
#include "util/ferm/subset_vectors.h"
#include "util/ferm/key_prop_colorvec.h"

namespace Chroma
{
  //! MapObj factory (foundry)
  /*! @ingroup ferm */
  typedef SingletonHolder< 
    ObjectFactory<QDP::MapObject<int,EVPair<LatticeColorVector> >, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  QDP::MapObject<int,EVPair<LatticeColorVector> >* (*)(XMLReader&,
								       const std::string&), 
		  StringFactoryError> >
  TheMapObjIntKeyColorEigenVecFactory;

  //! MapObj factory (foundry)
  /*! @ingroup ferm */
  typedef SingletonHolder< 
    ObjectFactory<QDP::MapObject<KeyPropColorVec_t,LatticeFermion>, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  QDP::MapObject<KeyPropColorVec_t,LatticeFermion>* (*)(XMLReader&,
									const std::string&), 
		  StringFactoryError> >
  TheMapObjKeyPropColorVecFactory;

}


#endif
