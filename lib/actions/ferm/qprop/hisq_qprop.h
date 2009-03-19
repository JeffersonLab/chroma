// -*- C++ -*-
// $Id: hisq_qprop.h,v 1.2 2009-03-19 13:08:05 mcneile Exp $
/*! \file
 *  \brief Hisq propagator wrapper
 *
 *  \ingroup qprop
 */
#ifndef HISQ_QPROP_H
#define HISQ_QPROP_H

/**
  Add the level3 code to the HISQ inverter.
  
  The fat and Naik links are created outide these
  routines and passed in here.
 
**/


// #else 
#include "chromabase.h"
#include "chroma_config.h"

#if defined(BUILD_CPS_ASQTAD_INVERTER)

#include "actions/ferm/qprop/asqtad_cps_wrapper_qprop.h"

namespace Chroma { 
  typedef AsqtadCPSWrapperQprop  HisqQprop ;
}

#else 


#include "actions/ferm/qprop/eoprec_staggered_qprop.h"

namespace Chroma {
  typedef EvenOddFermActQprop<LatticeStaggeredFermion,
    multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > HisqQprop;

}

#endif 

#endif
