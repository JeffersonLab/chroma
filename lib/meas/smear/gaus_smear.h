// -*- C++ -*-
// $Id: gaus_smear.h,v 1.1 2003-02-15 05:54:26 edwards Exp $

#ifndef GAUS_SMEAR_INCLUDE
#define GAUS_SMEAR_INCLUDE

void gausSmear(const multi1d<LatticeColorMatrix>& u, 
	       LatticeColorVector& chi, 
	       const Real& width, int ItrGaus, int j_decay);

#endif
