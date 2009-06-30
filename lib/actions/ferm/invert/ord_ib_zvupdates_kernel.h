#ifndef  ORD_IB_ZVUPDATES_KERNEL_H
#define  ORD_IB_ZVUPDATES_KERNEL_H

#ifdef BUILD_SSE_SCALARSITE_BICGSTAB
#include "ord_ib_zvupdates_kernel_sse.h"
#warning "Using SSE for BiCGStab Kernels"
#else
#include "ord_ib_zvupdates_kernel_generic.h"
#endif

#endif
