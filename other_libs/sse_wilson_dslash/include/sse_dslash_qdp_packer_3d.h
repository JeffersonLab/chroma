#ifndef SSE_DSLASH_QDP_PACKER_3D_H
#define SSE_DSLASH_QDP_PACKER_3D_H

#include "sse_config.h"

#ifndef QDP_INCLUDE
#include "qdp.h"
#endif 
using namespace QDP;

namespace SSEDslash3D { 

  // This is a little hacky, but switching this on for general Nc
  typedef PColorMatrix<RComplex<REAL>, Nc> PrimitiveSU3Matrix;

  void qdp_pack_gauge_3d(const multi1d<LatticeColorMatrix>&_u, multi1d<PrimitiveSU3Matrix>& u_tmp);

};

#endif
