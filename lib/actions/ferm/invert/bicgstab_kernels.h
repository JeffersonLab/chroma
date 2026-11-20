#ifndef BICGSTAB_KERNELS_H
#define BICGSTAB_KERNELS_H

#include "chroma_config.h"
// Switchbox file for bicgstab kernel specialisations
// Always include the naive ones

#include "actions/ferm/invert/bicgstab_kernels_naive.h"

#if defined(BUILD_SCALARSITE_BICGSTAB) && !defined(QDP_IS_QDPJIT)
#warning "Using Scalarsite BiCGStab Kernels"
#include "actions/ferm/invert/bicgstab_kernels_scalarsite.h"
#endif

namespace Chroma {

  namespace BiCGStabKernels { 

    inline
    void initKernels()
    {

#if defined(BUILD_SCALARSITE_BICGSTAB) && !defined(QDP_IS_QDPJIT)
      initScalarSiteKernels();
#endif


    }

    inline
    void finishKernels()
    {
#if defined(BUILD_SCALARSITE_BICGSTAB) && !defined(QDP_IS_QDPJIT)
      finishScalarSiteKernels();
#endif
    }
  }


}
#endif
