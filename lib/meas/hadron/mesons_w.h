// -*- C++ -*-
// $Id: mesons_w.h,v 1.1 2003-02-15 05:54:26 edwards Exp $

#ifndef MESONS_INCLUDE
#define MESONS_INCLUDE

void mesons(const LatticePropagator& quark_prop_1, const LatticePropagator& quark_prop_2, 
	    multi2d<Real>& meson_propagator, 
	    const multi1d<int>& t_source, int j_decay);

#endif
