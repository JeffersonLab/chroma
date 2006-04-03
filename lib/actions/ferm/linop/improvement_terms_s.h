#ifndef IMPROVEMENT_TERMS_S
#define IMPROVEMENT_TERMS_S


#include "chromabase.h"

namespace Chroma 
{ 
  void Fat7_Links(multi1d<LatticeColorMatrix>& u, multi1d<LatticeColorMatrix>& u_fat, Real u0);

  void Triple_Links(multi1d<LatticeColorMatrix>& u, multi1d<LatticeColorMatrix>& u_triple, Real u0);

} // End Namespace Chroma


#endif

