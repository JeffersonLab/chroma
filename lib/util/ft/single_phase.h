// -*- C++ -*-
//  $Id: single_phase.h,v 3.1 2006-11-27 20:07:58 edwards Exp $
/*! \file
 *  \brief Compute a single phase factor
 */

#ifndef __single_phase_h__
#define __single_phase_h__

#include "chromabase.h"

namespace Chroma 
{
  //! A single exp(ip.x) phase used in hadron construction
  /*! @ingroup ft */
  LatticeComplex singlePhase(const multi1d<int>& t_srce, 
			     const multi1d<int>& sink_mom, 
			     int j_decay);

}  // end namespace Chroma

#endif
