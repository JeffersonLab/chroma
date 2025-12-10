#ifndef SSE_DSLASH_QDP_PACKER_H
#define SSE_DSLASH_QDP_PACKER_H

#include "sse_config.h"

#ifndef QDP_INCLUDE
#include "qdp.h"
#endif 

using namespace QDP;

namespace SSEDslash { 

  // This is a little hacky, but switching this on for general Nc
  typedef PColorMatrix<RComplex<REAL>, Nc> PrimitiveSU3Matrix;


  void qdp_pack_gauge(const multi1d<LatticeColorMatrix>&_u, multi1d<PrimitiveSU3Matrix>& u_tmp);
  void qdp_pack_gauge3d(const multi1d<LatticeColorMatrix>&_u, multi1d<PrimitiveSU3Matrix>& u_tmp);

};

#endif
