// -*- C++ -*-
// $Id: baryon_w.h,v 1.1 2003-02-15 05:54:26 edwards Exp $

#ifndef BARYON_INCLUDE
#define BARYON_INCLUDE

void baryon(LatticePropagator& quark_propagator, 
	    multi2d<Complex>& barprop, 
	    const multi1d<int>& t_source, int j_decay, int bc_spec);

#endif
