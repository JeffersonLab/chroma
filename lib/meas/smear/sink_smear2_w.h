// -*- C++ -*-
// $Id: sink_smear2_w.h,v 1.2 2003-04-01 01:56:19 edwards Exp $

#ifndef __sink_smear2_h__
#define __sink_smear2_h__

void sink_smear2(const multi1d<LatticeColorMatrix>& u, 
		 LatticePropagator& quark_propagator, 
		 const Real& wvf_type, int wvf_param, int WvfIntPar, int j_decay);

#endif
