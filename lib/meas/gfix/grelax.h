// -*- C++ -*-
// $Id: grelax.h,v 3.1 2009-10-09 15:33:43 bjoo Exp $
/*! \file
 *  \brief Perform a single gauge fixing iteration
 */

#ifndef __grelax_h__
#define __grelax_h__

namespace Chroma {

//! Perform a single gauge fixing iteration
/*!
 * \ingroup gfix
 *
 * Performs one gauge fixing 'iteration', one checkerboard and SU(2)
 * subgroup only, for gauge fixing to Coulomb gauge in slices perpendicular
 * to the direction "j_decay".
 *
 * \param g          Current (global) gauge transformation matrices ( Modify )
 * \param u          original gauge field ( Read )
 * \param j_decay    direction perpendicular to slices to be gauge fixed ( Read )
 * \param su2_index  SU(2) subgroup index ( Read )
 * \param cb         checkerboard index ( Read )
 * \param ordo       use overrelaxation or not ( Read )
 * \param orpara     overrelaxation parameter ( Read ) 
 */

void grelax(LatticeColorMatrix& g,
	    const multi1d<LatticeColorMatrix>& u, 
	    int j_decay, int su2_index, int cb, bool ordo,
	    const Real& orpara);

}  // end namespace Chroma

#endif
