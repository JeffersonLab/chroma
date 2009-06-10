#ifndef  ORD_YXPAYMABZ_KERNEL_H
#define  ORD_YXPAYMABZ_KERNEL_H

#ifdef BUILD_SSE_SCALARSITE_BICGSTAB
#include "ord_yxpaymabz_kernel_sse.h"
#warning "Using SSE for BiCGStab Kernels"
#else
#include "ord_yxpaymabz_kernel_generic.h"
#endif

#endif
