// -*- C++ -*-
// $Id: block.h,v 3.1 2006-08-24 02:33:52 edwards Exp $
/*! \file
 *  \brief Construct "block" links
 */

#ifndef __block_h__
#define __block_h__

namespace Chroma 
{

  //! Construct block links
  /*!
   * \ingroup glue
   *
   * Construct block links from:
   *      x           x                       x------x
   *      |           |           _______            |
   *      |           |           \                  |
   *      |           |            \                 |
   *      |     =     x     +       \                x
   *      |           |             /                |
   *      |           |            /                 |
   *      |           |           /______            |
   *      x           x                       x------x

   * projected back onto SU(Nc)

   * Warning: this works only for Nc = 2 and 3 !

   * \param u_block    blocked gauge field ( Write )
   * \param u          gauge field ( Read )
   * \param mu         direction of blocked gauge field ( Read )
   * \param bl_level   blocking level (of the u's) ( Read )
   * \param BlkAccu    accuracy in fuzzy link projection ( Read )
   * \param BlkMax     maximum number of iterations in fuzzy link projection ( Read )
   * \param j_decay    no staple in direction j_decay ( Read ) 
   */

  void block(LatticeColorMatrix& u_block, 
	     const multi1d<LatticeColorMatrix>& u, 
	     int mu, int bl_level, 
	     const Real& BlkAccu, int BlkMax, int j_decay);

}  // end namespace Chroma

#endif
