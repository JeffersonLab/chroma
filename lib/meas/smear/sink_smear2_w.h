// -*- C++ -*-
// $Id: sink_smear2_w.h,v 1.1 2003-03-07 05:42:01 edwards Exp $

#ifndef SINK_SMEAR2_INCLUDE
#define SINK_SMEAR2_INCLUDE

void sink_smear2(const multi1d<LatticeColorMatrix>& u, 
		 LatticePropagator& quark_propagator, 
		 const Real& wvf_type, int wvf_param, int WvfIntPar, int j_decay);

#endif
