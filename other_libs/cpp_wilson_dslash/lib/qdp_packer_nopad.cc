#include "cpp_dslash_qdp_packer.h"

using namespace QDP;

namespace CPlusPlusWilsonDslash {

  template<typename T1, typename T2>
void qdpPackGauge(const multi1d<T1>&_u,  T2* u_tmp)
{
	int volume = Layout::sitesOnNode();

#pragma omp parallel for
	for(int ix = 0; ix < volume; ix++) 
	{
		for(int mu = 0; mu < 4; mu++) 
		{ 
			u_tmp[ mu + 4*(ix) ] =
				transpose( _u[mu].elem(ix).elem() );
		}
	}
}

 
void qdp_pack_gauge(const multi1d<LatticeColorMatrixF>&_u, PrimitiveSU3MatrixF* u_tmp)
  {
    qdpPackGauge<>(_u, u_tmp);
  }

void qdp_pack_gauge(const multi1d<LatticeColorMatrixD>&_u, PrimitiveSU3MatrixD* u_tmp)
  {
    qdpPackGauge<>(_u, u_tmp);
  }

}
