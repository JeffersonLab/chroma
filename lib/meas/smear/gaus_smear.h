// -*- C++ -*-
// $Id: gaus_smear.h,v 1.6 2005-01-14 18:42:37 edwards Exp $
/*! \file
 *  \brief Gaussian smearing of color vector
 */

#ifndef __gaus_smear_h__
#define __gaus_smear_h__

namespace Chroma {

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

}  // end namespace Chroma

#endif
