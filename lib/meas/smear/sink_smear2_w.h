// -*- C++ -*-
// $Id: sink_smear2_w.h,v 1.3 2003-04-25 20:46:31 flemingg Exp $

#ifndef __sink_smear2_h__
#define __sink_smear2_h__

void
sink_smear2(const multi1d<LatticeColorMatrix>& u,
            LatticePropagator& quark_propagator, int wvf_type,
            const Real& wvf_param, int WvfIntPar, int j_decay) ;

#endif
