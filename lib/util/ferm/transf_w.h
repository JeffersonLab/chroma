// -*- C++ -*-
// $Id: transf_w.h,v 1.1 2003-02-15 05:54:26 edwards Exp $

#ifndef TRANSF_INCLUDE
#define TRANSF_INCLUDE

void CvToFerm(const LatticeColorVector& a, LatticeFermion& b, 
	      int spin_index);

void FermToProp(const LatticeFermion& a, LatticePropagator& b, 
		int color_index, int spin_index);

void PropToFerm(const LatticePropagator& b, LatticeFermion& a, 
		int color_index, int spin_index);

#endif


