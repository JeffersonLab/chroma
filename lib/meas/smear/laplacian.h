// -*- C++ -*-
// $Id: laplacian.h,v 1.3 2003-10-10 03:46:47 edwards Exp $

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
