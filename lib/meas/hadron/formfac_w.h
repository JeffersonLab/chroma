// -*- C++ -*-
// $Id: formfac_w.h,v 1.1 2003-02-25 20:25:28 edwards Exp $

#ifndef FORMFAC_INCLUDE
#define FORMFAC_INCLUDE

void FormFac(const multi1d<LatticeColorMatrix>& u, 
	     const LatticePropagator& quark_propagator,
	     const LatticePropagator& seq_quark_prop, 
	     const multi1d<int>& t_source, 
	     int t_sink, int j_decay,
	     NmlWriter& nml);

#endif
