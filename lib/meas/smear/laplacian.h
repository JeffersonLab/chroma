// -*- C++ -*-
// $Id: laplacian.h,v 1.2 2003-06-19 17:40:42 ikuro Exp $

#ifndef __laplacian_h__
#define __laplacian_h__

void laplacian(const multi1d<LatticeColorMatrix>& u, 
	       LatticeColorVector& chi, 
	       int j_decay,
	       int power);

void laplacian(const multi1d<LatticeColorMatrix>& u, 
	       LatticePropagator& chi, 
	       int j_decay,
	       int power);


#endif
