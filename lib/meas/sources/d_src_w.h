/*
 *  P-wave fermion source
 */

#ifndef __d_src_h__
#define __d_src_h__

void d_src(const multi1d<LatticeColorMatrix>& u, 
	   LatticeFermion& chi,
	   int direction);

#endif
