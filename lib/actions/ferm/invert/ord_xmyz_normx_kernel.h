#ifndef  ORD_XMYZ_NORMX_KERNEL_H
#define  ORD_XMYZ_NORMX_KERNEL_H

#ifdef BUILD_SSE_SCALARSITE_BICGSTAB
#include "ord_xmyz_normx_kernel_sse.h"
#warning "Using SSE for BiCGStab Kernels"
#else
#include "ord_xmyz_normx_kernel_generic.h"
#endif

#endif
