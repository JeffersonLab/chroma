// -*- C++ -*-
#ifndef __map_obj_memory_w_h__
#define __map_obj_memory_w_h__

namespace Chroma {
  //! Private Namespace 
  namespace MapObjectMemoryEnv { 
    

    //! MapObjectMemory<KeyBlockProp_t,LatticeFermion> registrations
    bool registerKeyBlockPropLF();

    //! MapObjectMemory<KeyGridProp_t,LatticeFermion> registrations
    bool registerKeyGridPropLF();

    //! MapObjectMemory<KeyPropColorVec_t,LatticeFermion> registrations
    bool registerKeyPropColorVecLF();

  }
}

#endif
