// -*- C++ -*-
// $Id: barcomp_w.h,v 1.1 2003-03-08 03:57:47 edwards Exp $

#ifndef BARCOMP_INCLUDE
#define BARCOMP_INCLUDE

void barcomp(const LatticePropagator& quark_propagator_1, 
	     const LatticePropagator& quark_propagator_2, 
	     const LatticePropagator& quark_propagator_3, 
	     int t0, int j_decay, int bc_spec,
	     const string& nml_group,
	     NmlWriter& nml);

#endif
