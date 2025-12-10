#ifndef CPP_DSLASH_SCALAR_64BIT_H
#define CPP_DSLASH_SCALAR_64BIT_H

#include "dslash_config.h"
#ifdef DSLASH_USE_SSE2
#include "cpp_dslash_scalar_64bit_sse.h"
#else
#include "cpp_dslash_scalar_64bit_c.h"
#endif

#endif
