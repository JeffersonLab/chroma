#include "cpp_dslash_qdp_packer.h"

using namespace QDP;

namespace CPlusPlusWilsonDslash {

  /* Straightforward packing for now. No interleaving - just pack as
     u[ x ][ mu ] mu=0,1,2
  */
  template<typename T1, typename T2>
  void qdpPackGauge3D(const multi1d<T1>&_u, multi1d<T2>& u_tmp)
  {
    int Nd3=4;
    
    int volume = Layout::sitesOnNode();
    
    for(int ix = 0; ix < volume; ix++)  {
      for(int mu = 0; mu < Nd3; mu++) {
	u_tmp[ mu + Nd3*(ix) ] =
	  transpose( _u[mu].elem(ix).elem() );
      }
    }
  }

  void qdp_pack_gauge_3d(const multi1d<LatticeColorMatrixF>&_u, multi1d<PrimitiveSU3MatrixF>& u_tmp)
  {
    qdpPackGauge3D<>(_u, u_tmp);
  }

  void qdp_pack_gauge_3d(const multi1d<LatticeColorMatrixD>&_u, multi1d<PrimitiveSU3MatrixD>& u_tmp)
  {
    qdpPackGauge3D<>(_u, u_tmp);
  }
}
