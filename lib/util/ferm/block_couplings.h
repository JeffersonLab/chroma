// -*- C++ -*-
// $Id: block_couplings.h,v 1.2 2009-01-30 21:07:21 kostas Exp $
/*! \file
 *  \brief Caclulates the couplings between neighboring blocks given a displacement path
 */

#ifndef __block_couplings_h__
#define __block_couplings_h__
#include <vector> 

#include "chromabase.h"

namespace Chroma
{
  vector<int> block_couplings(const int b,
			      const Set& S,const multi1d<int>& disp,
			      const int len);

} // namespace Chroma

#endif
