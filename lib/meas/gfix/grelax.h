// -*- C++ -*-
// $Id: grelax.h,v 1.3 2004-01-02 22:48:02 edwards Exp $
/*! \file
 *  \brief Perform a single gauge fixing iteration
 */

#ifndef __grelax_h__
#define __grelax_h__

//! Perform a single gauge fixing iteration
/*!
 * \ingroup gfix
 *
 * Performs one gauge fixing 'iteration', one checkerboard and SU(2)
 * subgroup only, for gauge fixing to Coulomb gauge in slices perpendicular
 * to the direction "j_decay".
 *
 * \param ug         (gauge fixed) gauge field ( Modify )
 * \param u_neg      (gauge fixed) gauge field, negative links ( Read )
 * \param j_decay    direction perpendicular to slices to be gauge fixed ( Read )
 * \param su2_index  SU(2) subgroup index ( Read )
 * \param cb         checkerboard index ( Read )
 * \param ordo       use overrelaxation or not ( Read )
 * \param orpara     overrelaxation parameter ( Read ) 
 */

void grelax(multi1d<LatticeColorMatrix>& ug, 
	    const multi1d<LatticeColorMatrix>& u_neg,
	    int j_decay, int su2_index, int cb, bool ordo,
	    const Real& orpara);

#endif
