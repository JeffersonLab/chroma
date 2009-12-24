// -*- C++ -*-
#ifndef __map_obj_disk_traits_h__
#define __map_obj_disk_traits_h__

#include "chromabase.h"
#include "util/ferm/key_prop_colorvec.h"
#include "util/ferm/subset_ev_pair.h"
#include <string>

namespace Chroma {

  /* Type strings are initialized in map_obj_disk.cc */

  //! Type trait to take Key/Value pair to unsigned int
  template<typename K, typename V> 
  struct MapObjTraitsNum {};
  

  /*---------  KeyPropColorVec_t/LatticeFermion -----------------*/
  template<>
  struct MapObjTraitsNum<KeyPropColorVec_t, LatticeFermion> {
    static const std::string type_string;
  };

  /* -------- int / LatticeColorVec ---------------------------*/
  template<>
  struct MapObjTraitsNum<int, EVPair<LatticeColorVector> >{
    static const std::string type_string;
  };

  /* -------- char / float (test ) ---------------------------*/
  template<>
  struct MapObjTraitsNum<char, float> {
    static const std::string type_string;
  };


}



#endif
