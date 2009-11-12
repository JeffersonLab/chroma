// -*- C++ -*-
/*
 * \file map_obj_aggregate_w.h
 * \brief Header file for map obj aggregate registrations 
 */

#ifndef __map_obj_aggregate_w_h__
#define __map_obj_aggregage_w_h__

#include "chromabase.h"

namespace Chroma {

  //! Registration aggregator
  namespace MapObjectWilson4DEnv {

    //! aggregate for  MapObject<KeyBlockProp_t, LatticeFermion> factory
    bool registerKeyBlockPropLFAll();

    //! aggregate for  MapObject<KeyGridProp_t, LatticeFermion> factory
    bool registerKeyGridPropLFAll();

    //! aggregate for MapObject<KeyPropColorVec_t, LatticeFermion> factory
    bool registerKeyPropColorVecLFAll();

    //! aggregate everything
    bool registerAll();

  }

}



#endif
