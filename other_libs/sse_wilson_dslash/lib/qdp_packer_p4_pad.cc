#include "sse_dslash_qdp_packer.h"


using namespace QDP;

namespace SSEDslash { 
void qdp_pack_gauge(const multi1d<LatticeColorMatrix>&_u, multi1d<PrimitiveSU3Matrix>& u_tmp)
{
  multi1d<PrimitiveSU3Matrix> v(8);
  int volume = Layout::sitesOnNode();
  
  for(int ix = 0; ix < volume; ix+=2) 
  {
    for(int mu = 0; mu < 4; mu++) 
    { 
      v[2*mu] = transpose(_u[mu].elem(ix).elem());
      v[2*mu+1] = transpose(_u[mu].elem(ix + 1).elem());
    }

    for(int mu = 0; mu < 4; mu++) 
    {
      u_tmp[ mu + Nd*(ix) ] = v[mu];
      u_tmp[ mu + Nd*(ix + 1) ] = v[Nd + mu]; 
    }
  }
}

}
