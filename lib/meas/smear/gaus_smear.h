// -*- C++ -*-
// $Id: gaus_smear.h,v 1.5 2003-12-27 04:27:55 edwards Exp $
/*! \file
 *  \brief Gaussian smearing of color vector
 */

#ifndef __gaus_smear_h__
#define __gaus_smear_h__

//! Do a covariant Gaussian smearing of a lattice color vector field
/*! This is a wrapper over the template definition
 *
 * \ingroup smear
 *
 * Arguments:
 *
 *  \param u        gauge field ( Read )
 *  \param chi      color vector field ( Modify )
 *  \param width    width of "shell" wave function ( Read )
 *  \param ItrGaus  number of iterations to approximate Gaussian ( Read )
 *  \param j_decay  direction of decay ( Read )
 */
void gausSmear(const multi1d<LatticeColorMatrix>& u, 
	       LatticeColorVector& chi, 
	       const Real& width, int ItrGaus, int j_decay);


//! Do a covariant Gaussian smearing of a lattice propagator field
/*! This is a wrapper over the template definition
 *
 * \ingroup smear
 *
 * Arguments:
 *
 *  \param u        gauge field ( Read )
 *  \param chi      propagator field ( Modify )
 *  \param width    width of "shell" wave function ( Read )
 *  \param ItrGaus  number of iterations to approximate Gaussian ( Read )
 *  \param j_decay  direction of decay ( Read )
 */
void gausSmear(const multi1d<LatticeColorMatrix>& u, 
	       LatticePropagator& chi, 
	       const Real& width, int ItrGaus, int j_decay);

#endif
