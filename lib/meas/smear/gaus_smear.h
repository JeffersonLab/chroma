// -*- C++ -*-
// $Id: gaus_smear.h,v 1.3 2003-04-01 01:56:19 edwards Exp $

#ifndef __gaus_smear_h__
#define __gaus_smear_h__

void gausSmear(const multi1d<LatticeColorMatrix>& u, 
	       LatticeColorVector& chi, 
	       const Real& width, int ItrGaus, int j_decay);

void gausSmear(const multi1d<LatticeColorMatrix>& u, 
	       LatticePropagator& chi, 
	       const Real& width, int ItrGaus, int j_decay);

#endif
