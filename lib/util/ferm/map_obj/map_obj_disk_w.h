// -*- C++ -*-
#ifndef __map_obj_disk_w_h__
#define __map_obj_disk_w_h__

namespace Chroma {
  //! Private Namespace 
  namespace MapObjectDiskEnv { 
    

    //! MapObjectDisk<KeyBlockProp_t,LatticeFermion> registrations
    bool registerKeyBlockPropLF();

    //! MapObjectDisk<KeyGridProp_t,LatticeFermion> registrations
    bool registerKeyGridPropLF();

    //! MapObjectDisk<KeyPropColorVec_t,LatticeFermion> registrations
    bool registerKeyPropColorVecLF();

  }
}

#endif
