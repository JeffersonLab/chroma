// -*- C++ -*-
// $Id: gaus_smear.h,v 1.2 2003-03-07 05:22:01 edwards Exp $

#ifndef GAUS_SMEAR_INCLUDE
#define GAUS_SMEAR_INCLUDE

void gausSmear(const multi1d<LatticeColorMatrix>& u, 
	       LatticeColorVector& chi, 
	       const Real& width, int ItrGaus, int j_decay);

void gausSmear(const multi1d<LatticeColorMatrix>& u, 
	       LatticePropagator& chi, 
	       const Real& width, int ItrGaus, int j_decay);

#endif
