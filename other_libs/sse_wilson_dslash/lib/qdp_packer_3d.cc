#include "sse_dslash_qdp_packer_3d.h"

using namespace QDP;

namespace SSEDslash3D {

  /* Straightforward packing for now. No interleaving - just pack as
     u[ x ][ mu ] mu=0,1,2
  */

void qdp_pack_gauge_3d(const multi1d<LatticeColorMatrix>&_u, multi1d<PrimitiveSU3Matrix>& u_tmp)
{
  int Nd3=4;

  int volume = Layout::sitesOnNode();
  
  for(int ix = 0; ix < volume; ix++) 
  {
    for(int mu = 0; mu < Nd3; mu++) 
    { 
      u_tmp[ mu + Nd3*(ix) ] =
	transpose( _u[mu].elem(ix).elem() );
    }
  }
}

}
