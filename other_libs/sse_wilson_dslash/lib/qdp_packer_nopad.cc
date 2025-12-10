#include "sse_dslash_qdp_packer.h"

using namespace QDP;

namespace SSEDslash {

void qdp_pack_gauge(const multi1d<LatticeColorMatrix>&_u, multi1d<PrimitiveSU3Matrix>& u_tmp)
{
  int volume = Layout::sitesOnNode();
  
  for(int ix = 0; ix < volume; ix++) 
  {
    for(int mu = 0; mu < 4; mu++) 
    { 
      u_tmp[ mu + 4*(ix) ] =
	transpose( _u[mu].elem(ix).elem() );
    }
  }
}

}
