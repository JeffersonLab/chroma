// -*- C++ -*-
// $Id: klein_gord.h,v 1.2 2003-10-10 03:46:46 edwards Exp $

#ifndef KLEIN_GORD_INCLUDE
#define KLEIN_GORD_INCLUDE

void klein_gord(const multi1d<LatticeColorMatrix>& u, 
		const LatticeColorVector& psi, 
		LatticeColorVector& chi, 
		const Real& mass_sq, int j_decay);

void klein_gord(const multi1d<LatticeColorMatrix>& u, 
		const LatticePropagator& psi, 
		LatticePropagator& chi, 
		const Real& mass_sq, int j_decay);

#endif
