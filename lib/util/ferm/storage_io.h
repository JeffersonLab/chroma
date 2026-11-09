#pragma once
#include "chromabase.h"
#include <datalib/storage/record_store.h>
#ifdef BUILD_SB
// Match the CORR reader profile: tensor-storage class layouts depend on this.
# ifndef USE_SUPERBBLAS
#  define USE_SUPERBBLAS
# endif
# ifndef USE_SUPERBBLAS_WITH_GPU_SUPPORT
#  define USE_SUPERBBLAS_WITH_GPU_SUPPORT
# endif
#endif
namespace Chroma {
template<class T> std::string encodeStorageRecord(const T& value) {
  QDP::LocalBinaryBufferWriter out;
  write(out, value);
  return out.strPrimaryNode();
}
}
