// -*- C++ -*-
// $Id: transf.h,v 1.1 2003-12-11 17:11:17 bjoo Exp $

#ifndef TRANSF_INCLUDE
#define TRANSF_INCLUDE

void CvToFerm(const LatticeColorVector& a, LatticeFermion& b, 
	      int spin_index);

void FermToProp(const LatticeFermion& a, LatticePropagator& b, 
		int color_index, int spin_index);

void PropToFerm(const LatticePropagator& b, LatticeFermion& a, 
		int color_index, int spin_index);

#endif


