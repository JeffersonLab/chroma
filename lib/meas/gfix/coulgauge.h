// -*- C++ -*-
// $Id: coulgauge.h,v 1.1 2003-12-06 20:56:56 edwards Exp $
/*! \file
 *  \brief Coulomb (and Landau) gauge fixing 
 */

#ifndef __coulgauge_h__
#define __coulgauge_h__

//! Coulomb (and Landau) gauge fixing
/*!
 * \ingroup gfix
 *
 * Driver for gauge fixing to Coulomb gauge in slices perpendicular
 * to the direction "j_decay".
 * If j_decay >= Nd: fix to Landau gauge.
 * Note: as written this works only for SU(2) and SU(3)!

 * \param ug       (gauge fixed) gauge field ( Modify )
 * \param n_gf     number of gauge fixing iterations ( Write )
 * \param j_decay  direction perpendicular to slices to be gauge fixed ( Read )
 * \param GFAccu   desired accuracy for gauge fixing ( Read )
 * \param GFMax    maximal number of gauge fixing iterations ( Read )
 * \param OrDo     use overrelaxation or not ( Read )
 * \param OrPara   overrelaxation parameter ( Read )
 */

void gfix(multi1d<LatticeColorMatrix>& ug, int& n_gf, 
	  int j_decay, const Real& GFAccu, int GFMax, 
	  bool OrDo, const Real& OrPara);

#endif
