#ifndef  ORD_XMAY_NORMX_CDOTZX_KERNEL_H
#define  ORD_XMAY_NORMX_CDOTZX_KERNEL_H

#include "qdp_config.h"

#ifdef BUILD_SSE_SCALARSITE_BICGSTAB
#include "ord_xmay_normx_cdotzx_kernel_sse.h"
#warning "Using SSE for BiCGStab Kernels"
#else
#include "ord_xmay_normx_cdotzx_kernel_generic.h"
#endif

#endif
