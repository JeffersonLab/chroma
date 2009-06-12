#ifndef  ORD_CXMAYF_KERNEL_H
#define  ORD_CXMAYF_KERNEL_H

#ifdef BUILD_SSE_SCALARSITE_BICGSTAB
#include "ord_cxmayf_kernel_sse.h"
#warning "Using SSE for BiCGStab Kernels"
#else
#include "ord_cxmayf_kernel_generic.h"
#endif

#endif
