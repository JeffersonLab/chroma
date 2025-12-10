#ifndef CPP_CLOVER_SITE_APPLY_64_BIT_H
#define CPP_CLOVER_SITE_APPLY_64_BIT_H

#include <dslash_config.h>
#if (DSLASH_USE_SSE2 || DSLASH_USE_SSE3)
#include "cpp_clover_site_apply_64bit_sse.h"
#else
#include "cpp_clover_site_apply_64bit_c.h"
#endif

#endif
