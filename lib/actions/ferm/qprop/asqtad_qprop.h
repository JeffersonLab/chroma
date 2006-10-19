// -*- C++ -*-
// $Id: asqtad_qprop.h,v 3.2 2006-10-19 17:36:07 edwards Exp $
/*! \file
 *  \brief Asqtad propagator wrapper
 *
 *  \ingroup qprop
 */
#ifndef ASQTAD_QPROP_H
#define ASQTAD_QPROP_H


// #else 
#include "chromabase.h"
#include "chroma_config.h"



#if defined(BUILD_CPS_ASQTAD_INVERTER)

#include "actions/ferm/qprop/asqtad_cps_wrapper_qprop.h"

namespace Chroma { 
  typedef AsqtadCPSWrapperQprop  AsqtadQprop;
}

#else 

#include "actions/ferm/qprop/eoprec_staggered_qprop.h"

namespace Chroma {
  typedef EvenOddFermActQprop<LatticeStaggeredFermion,
    multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > AsqtadQprop;

}

#endif 
// #endif

#endif
