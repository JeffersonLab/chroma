#ifndef BICGSTAB_KERNELS_H
#define BICGSTAB_KERNELS_H

#include "chroma_config.h"
// Switchbox file for bicgstab kernel specialisations
// Always include the naive ones

#include "actions/ferm/invert/bicgstab_kernels_naive.h"

#ifdef BUILD_SCALARSITE_BICGSTAB
#warning "Using Scalarsite BiCGStab Kernels"
#include "actions/ferm/invert/bicgstab_kernels_scalarsite.h"
#endif

namespace Chroma {

  namespace BiCGStabKernels { 

    inline
    void initKernels()
    {

#ifdef BUILD_SCALARSITE_BICGSTAB
      initScalarSiteKernels();
#endif


    }

    inline
    void finishKernels()
    {
#ifdef BUILD_SCALARSITE_BICGSTAB
      finishScalarSiteKernels();
#endif
    }
  }


}
#endif
