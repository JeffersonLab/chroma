/*
 *  P-wave fermion source
 */

#ifndef __d_src_h__
#define __d_src_h__

namespace Chroma {

void d_src(const multi1d<LatticeColorMatrix>& u, 
	   LatticeColorVector& chi,
	   int direction);

void d_src(const multi1d<LatticeColorMatrix>& u, 
	   LatticePropagator& chi,
	   int direction);

void d_src(const multi1d<LatticeColorMatrix>& u, 
	   LatticeFermion& chi,
	   int direction);

}  // end namespace Chroma
#endif
