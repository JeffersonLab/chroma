// -*- C++ -*-
// $Id: hisq_qprop.h,v 1.1 2007-05-09 12:45:31 mcneile Exp $
/*! \file
 *  \brief Hisq propagator wrapper
 *
 *  \ingroup qprop
 */
#ifndef HISQ_QPROP_H
#define HISQ_QPROP_H


// #else 
#include "chromabase.h"
#include "chroma_config.h"



#include "actions/ferm/qprop/eoprec_staggered_qprop.h"

namespace Chroma {
  typedef EvenOddFermActQprop<LatticeStaggeredFermion,
    multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > HisqQprop;

}


#endif
