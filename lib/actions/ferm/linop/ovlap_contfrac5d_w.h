// -*- C++ -*-
// $Id: ovlap_contfrac5d_w.h,v 3.1 2006-10-19 16:01:31 edwards Exp $
/*! \file
 *  \brief Include possibly optimized partfrac5d
 */

#ifndef __ovlap_contfrac5d_w_h__
#define __ovlap_contfrac5d_w_h__

#include "chroma_config.h"

// The QDP naive dslash class: QDPWilsonDslash
#include "eoprec_ovlap_contfrac5d_linop_array_w.h"

// The following is an ifdef lis that switches in optimised
// Dslash-es.

#ifdef BUILD_CFZ_CFZ_LINOP
#include "eoprec_ovlap_contfrac5d_linop_array_opt_w.h"
namespace Chroma {
typedef OptEvenOddPrecOvlapContFrac5DLinOpArray EvenOddPrecOvlapContFrac5DLinOpArray;
}  // end namespace Chroma


#else

// Bottom line, if no optimised Dslash-s exist then the naive QDP Dslash
// becomes the EvenOddPrecOvlapContFrac5DLinOpArray
namespace Chroma {
typedef QDPEvenOddPrecOvlapContFrac5DLinOpArray EvenOddPrecOvlapContFrac5DLinOpArray;
}  // end namespace Chroma
#endif


#endif
