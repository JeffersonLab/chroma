#ifndef SSE_CONFIG_H
#define SSE_CONFIG_H

/* Undef the unwanted from the environment -- eg the compiler command line */
#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION

#include "sse_config_internal.h"

/* Prefix everything with SSE */
#define SSE_PACKAGE PACKAGE
#define SSE_PACKAGE_BUGREPORT PACKAGE_BUGREPORT
#define SSE_PACKAGE_NAME PACKAGE_NAME
#define SSE_PACKAGE_STRING PACKAGE_STRING
#define SSE_PACKAGE_TARNAME PACKAGE_TARNAME
#define SSE_PACKAGE_VERSION PACKAGE_VERSION
#define SSE_VERSION VERSION
                                                                                
                                                                                
/* Undef the unwanted */
#undef PACKAGE
#undef PACKAGE_BUGREPORT
#undef PACKAGE_NAME
#undef PACKAGE_STRING
#undef PACKAGE_TARNAME
#undef PACKAGE_VERSION
#undef VERSION

/* Define the real type */
#if SSE_PRECISION == 32
typedef float SSEREAL;
#elif SSE_PRECISION == 64
typedef double SSEREAL;
#else
#error "Precision Not supported, define SSE_PRECISION"
#endif

#endif
