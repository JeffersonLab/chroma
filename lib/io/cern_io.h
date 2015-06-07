// $Id: cppacs_io.h,v 3.0 2006/04/03 04:58:55 edwards Exp $

/*! \file
 *  \brief CPPACS gauge format routines
 */

#ifndef __cern_io_h__
#define __cern_io_h__

#include "chromabase.h"

namespace Chroma {


//! CERN gauge field header
struct CERNGauge_t
{
  multi1d<int> nrow;      // Lattice size
  double  plaq;           // plaquette value
};


 void readCERN(multi1d<LatticeColorMatrix>&, const std::string& );

}  // end namespace Chroma

#endif
