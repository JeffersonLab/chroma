#ifndef __clov_triang_qdp_w_h__
#define __clov_triang_qdp_w_h__

#include <chromabase.h>

namespace Chroma
{

  //! Special structure used for triangular objects
  template <typename R>
  struct PrimitiveClovTriang {
    RScalar<R> diag[2][2 * Nc];
    RComplex<R> offd[2][2 * Nc * Nc - Nc];
  };

  template <typename R>
  struct QUDAPackedClovSite {
    R diag1[6];
    R offDiag1[15][2];
    R diag2[6];
    R offDiag2[15][2];
  };

}

#endif
