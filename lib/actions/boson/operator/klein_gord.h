// -*- C++ -*-
// $Id: klein_gord.h,v 1.4 2005-01-14 20:13:04 edwards Exp $

#ifndef KLEIN_GORD_INCLUDE
#define KLEIN_GORD_INCLUDE

namespace Chroma
{
  void klein_gord(const multi1d<LatticeColorMatrix>& u, 
		  const LatticeColorVector& psi, 
		  LatticeColorVector& chi, 
		  const Real& mass_sq, int j_decay);

  void klein_gord(const multi1d<LatticeColorMatrix>& u, 
		  const LatticePropagator& psi, 
		  LatticePropagator& chi, 
		  const Real& mass_sq, int j_decay);
}


#endif
