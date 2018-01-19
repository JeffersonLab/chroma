#ifndef QDPJIT_MEMORY_WRAPPER_H
#define QDPJIT_MEMORY_WRAPPER_H

#include "chromabase.h"
#include "qdp_config.h"

#ifdef QDP_IS_QDPJIT
#ifdef QDPJIT_IS_QDPJITNVVM
// NVVM Wrapper
#define  GetMemoryPtr(a)  QDP_get_global_cache().getDevicePtr( a )
#else
// OLD JIT/PTX wrapper
#define  GetMemoryPtr(a)  QDPCache::Instance().getDevicePtr( a )
#endif

#endif /* QDP_IS_QDPJIT */

#endif /* QDPJIT_MEMORY_WRAPPER */ 
