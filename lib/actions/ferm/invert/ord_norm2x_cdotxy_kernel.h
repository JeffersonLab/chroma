#ifndef  ORD_NORM2X_CDOTXY_KERNEL_H
#define  ORD_NORM2X_CDOTXY_KERNEL_H

#ifdef BUILD_SSE_SCALARSITE_BICGSTAB
#include "ord_norm2x_cdotxy_kernel_sse.h"
#warning "Using SSE for BiCGStab Kernels"
#else
#include "ord_norm2x_cdotxy_kernel_generic.h"
#endif

#endif
