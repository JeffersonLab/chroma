// -*- C++ -*-
// $Id: klein_gord.h,v 1.1 2003-04-09 05:57:56 edwards Exp $

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
