// -*- C++ -*-
// $Id: formfac_w.h,v 1.6 2003-06-25 16:11:44 edwards Exp $

#ifndef __formfac_h__
#define __formfac_h__

// void FormFac(const multi1d<LatticeColorMatrix>& u, 
//              const LatticePropagator& quark_propagator,
//              const LatticePropagator& seq_quark_prop, 
//              const multi1d<int>& t_source, 
//              int source_mom2_max,
//              int t_sink,
//              const multi1d<int>& sink_mom,
//              int j_decay,
//              NmlWriter& nml);

void FormFac(const multi1d<LatticeColorMatrix>& u, 
             const LatticePropagator& quark_propagator,
             const LatticePropagator& seq_quark_prop, 
             const SftMom& phases,
             int t0,
             XMLWriter& nml);

#endif
