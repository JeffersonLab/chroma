// -*- C++ -*-
// $Id: block_couplings.h,v 1.5 2009-02-23 19:50:15 edwards Exp $
/*! \file
 *  \brief Caclulates the couplings between neighboring blocks given a displacement path
 */

#ifndef __block_couplings_h__
#define __block_couplings_h__

#include <vector> 
#include "chromabase.h"

namespace Chroma
{
  //! Holds info on the coupling of sites that are displaced within a block
  /*!
   * \ingroup ferm
   * @{
   */
  class DisplacedBlock
  {
  public:
    int blk;
    multi1d<int> disp;
    DisplacedBlock(int b, multi1d<int> d) : blk(b),disp(d) {}
    DisplacedBlock() : blk(0) {}
  };

  vector<int> block_couplings(const int b,
			      const Set& S,const multi1d<int>& disp,
			      const int len);

  bool blocks_couple(const multi1d<DisplacedBlock>& b,
		     const Set& S,
		     const int len,int f);

  bool blocks_couple(const multi1d<DisplacedBlock>& b,
                     const Set& S,
                     const int len) ;

  /*! @} */  // end of group ferm
} // namespace Chroma

#endif
