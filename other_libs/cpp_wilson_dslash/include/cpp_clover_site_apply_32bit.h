#ifndef CPP_CLOVER_SITE_APPLY_32_BIT_H
#define CPP_CLOVER_SITE_APPLY_32_BIT_H


#include <dslash_config.h>
#if (DSLASH_USE_SSE2 || DSLASH_USE_SSE3)
// No single prec SSE yet.
#include "cpp_clover_site_apply_32bit_sse.h"
#else
#include "cpp_clover_site_apply_32bit_c.h"
#endif

#endif
