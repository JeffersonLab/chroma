// -*- C++ -*-
// $Id: barcomp_w.h,v 1.2 2003-04-01 01:56:19 edwards Exp $

#ifndef __barcomp_h__
#define __barcomp_h__

void barcomp(const LatticePropagator& quark_propagator_1, 
	     const LatticePropagator& quark_propagator_2, 
	     const LatticePropagator& quark_propagator_3, 
	     int t0, int j_decay, int bc_spec,
	     const string& nml_group,
	     NmlWriter& nml);

#endif
