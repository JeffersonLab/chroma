// -*- C++ -*-
// $Id: barcomp_w.h,v 1.3 2003-04-17 20:08:20 dgr Exp $

#ifndef __barcomp_h__
#define __barcomp_h__

void barcomp(const LatticePropagator& quark_propagator_1, 
	     const PropHead& header_1,
	     const LatticePropagator& quark_propagator_2,
	     const PropHead& header_2,
	     const LatticePropagator& quark_propagator_3,
	     const PropHead& header_3,
	     int t0, int j_decay, int bc_spec,
	     const char file[],
	     NmlWriter& nml);

#endif
