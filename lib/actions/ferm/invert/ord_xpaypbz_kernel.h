#ifndef  ORD_XPAYPBZ_KERNEL_H
#define  ORD_XPAYPBZ_KERNEL_H

#ifdef BUILD_SSE_SCALARSITE_BICGSTAB
#include "ord_xpaypbz_kernel_sse.h"
#warning "Using SSE for BiCGStab Kernels"
#else
#include "ord_xpaypbz_kernel_generic.h"
#endif

#endif
