// -*- C++ -*-
// $Id: laplacian.h,v 1.1 2003-05-27 17:45:29 ikuro Exp $

#ifndef __laplacian_h__
#define __laplacian_h__

void laplacian(const multi1d<LatticeColorMatrix>& u, 
	       LatticeColorVector& chi, 
	       int j_decay);

void laplacian(const multi1d<LatticeColorMatrix>& u, 
	       LatticePropagator& chi, 
	       int j_decay);


#endif
