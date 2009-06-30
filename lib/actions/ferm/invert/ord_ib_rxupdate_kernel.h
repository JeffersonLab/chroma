#ifndef  ORD_IB_RXUPDATE_KERNEL_H
#define  ORD_IB_RXUPDATE_KERNEL_H

#ifdef BUILD_SSE_SCALARSITE_BICGSTAB
#include "ord_ib_rxupdate_kernel_sse.h"
#warning "Using SSE for BiCGStab Kernels"
#else
#include "ord_ib_rxupdate_kernel_generic.h"
#endif

#endif
