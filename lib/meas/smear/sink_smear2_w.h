// -*- C++ -*-
// $Id: sink_smear2_w.h,v 1.5 2003-10-10 03:46:47 edwards Exp $

#ifndef __sink_smear2_h__
#define __sink_smear2_h__

#include "meas/smear/gaus_smear.h"

enum WvfKind {
  WVF_KIND_GAUGE_INV_GAUSSIAN,
  WVF_KIND_WUPPERTAL,
  WVF_KIND_UNKNOWN
} ;

void
sink_smear2(const multi1d<LatticeColorMatrix>& u,
            LatticePropagator& quark_propagator, WvfKind Wvf_kind,
            const Real& wvf_param, int WvfIntPar, int j_decay) ;

#endif
