// -*- C++ -*-
#ifndef __map_obj_disk_traits_h__
#define __map_obj_disk_traits_h__

#include "chromabase.h"
#include "util/ferm/key_prop_colorvec.h"
#include "util/ferm/key_block_prop.h"
#include "util/ferm/key_grid_prop.h"
#include "util/ferm/subset_ev_pair.h"
#include <string>

namespace Chroma {

  //! Type trait to take Key/Value pair to unsigned int
  template<typename K, typename V> 
  struct MapObjTraitsNum {};
  
  //! Type trait to take int to Key/Value pair 
  template<unsigned int N> 
  struct MapObjTraitsKeyVal {};
  

  /*---------  KeyPropColorVec_t/LatticeFermion -----------------*/


  template<>
  struct MapObjTraitsNum<KeyPropColorVec_t, LatticeFermion> {
    static const unsigned int filenum=1;
    //    static std::string type_string("KeyTKeyPropColorVec_tValTLatticeFermion");

  };

  template<>
  struct MapObjTraitsKeyVal<1> {
    typedef KeyPropColorVec_t KeyType_t;
    typedef LatticeFermion    ValType_t;
  };

  /*--------- KeyBlockProp_t/LatticeFermion ---------------------*/


  template<>
  struct MapObjTraitsNum<KeyBlockProp_t, LatticeFermion> {
    static const unsigned int filenum=2;
    // static std::string type_string("KeyTKeyBlockProp_tValTLatticeFermion");

  };

  template<>
  struct MapObjTraitsKeyVal<2> {
    typedef KeyBlockProp_t KeyType_t;
    typedef LatticeFermion ValType_t;
  };

  /* -------- KeyGridProp_t/LatticeFermion ----------------------*/

  
  template<>
  struct MapObjTraitsNum<KeyGridProp_t, LatticeFermion> { 
    static const unsigned int filenum=3;
    // static std::string type_string("KeyTKeyGridProp_tValTLatticeFermion");
  };

  template<>
  struct MapObjTraitsKeyVal<3> {
    typedef KeyGridProp_t KeyType_t;
    typedef LatticeFermion ValType_t;
  };


  /* -------- int / LatticeColorVec ---------------------------*/
  template<>
  struct MapObjTraitsNum<int, EVPair<LatticeColorVector> >{
    static const unsigned int filenum=4;
    // static std::string type_string("KeyTintValTEVPairLatticeColorVector");
  };

  template<>
  struct MapObjTraitsKeyVal<4> {
    typedef int                KeyType_t;
    typedef EVPair<LatticeColorVector> ValType_t;
  };


  /* -------- char / float (test ) ---------------------------*/
  template<>
  struct MapObjTraitsNum<char, float> {
    static const unsigned int filenum=5;
    // static std::string type_string("KeyTcharValTfloat");
  };

  template<>
  struct MapObjTraitsKeyVal<5> {
    typedef char                KeyType_t;
    typedef float               ValType_t;
  };


}



#endif
