// -*- C++ -*-
// $Id: formfac_w.h,v 1.2 2003-03-02 02:30:31 flemingg Exp $

#ifndef FORMFAC_INCLUDE
#define FORMFAC_INCLUDE

void FormFac(const multi1d<LatticeColorMatrix>& u, 
	     const LatticePropagator& quark_propagator,
	     const LatticePropagator& seq_quark_prop, 
	     const multi1d<int>& t_source, 
	     int source_mom2_max,
	     int t_sink,
	     const multi1d<int>& sink_mom,
	     int j_decay,
	     NmlWriter& nml);

#endif
