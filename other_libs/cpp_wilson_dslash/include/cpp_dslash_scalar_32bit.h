#ifndef CPP_DSLASH_SCALAR_32BIT_H
#define CPP_DSLASH_SCALAR_32BIT_H

#include "dslash_config.h"
#ifdef DSLASH_USE_SSE2
#include "cpp_dslash_scalar_32bit_sse.h"
#else
#include "cpp_dslash_scalar_32bit_c.h"
#endif

#endif
