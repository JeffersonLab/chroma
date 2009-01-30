// -*- C++ -*-
// $Id: block_couplings.h,v 1.1 2009-01-30 20:52:47 kostas Exp $
/*! \file
 *  \brief Caclulates the couplings between neighboring blocks given a displacement path
 */

#ifndef __block_couplings_h__
#define __block_couplings_h__
#include <vector> 

#include "chromabase.h"

namespace Chroma
{
  vector<int> block_coublings(const int b,
			      const Set& S,const multi1d<int>& disp,
			      const int len);

} // namespace Chroma

#endif
