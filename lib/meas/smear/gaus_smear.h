// -*- C++ -*-
// $Id: gaus_smear.h,v 1.4 2003-10-10 03:46:47 edwards Exp $

#ifndef __gaus_smear_h__
#define __gaus_smear_h__

void gausSmear(const multi1d<LatticeColorMatrix>& u, 
	       LatticeColorVector& chi, 
	       const Real& width, int ItrGaus, int j_decay);

void gausSmear(const multi1d<LatticeColorMatrix>& u, 
	       LatticePropagator& chi, 
	       const Real& width, int ItrGaus, int j_decay);

#endif
