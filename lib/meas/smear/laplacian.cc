// $Id: laplacian.cc,v 1.1 2003-05-27 17:45:29 ikuro Exp $
/*! \file
 *  \brief Gaussian smearing of color vector
 */

#include "chromabase.h"
#include "meas/smear/laplacian.h"
#include "actions/boson/operator/klein_gord.h"

using namespace QDP;

//! Do a covariant Gaussian smearing of a lattice field
/*!
 * Arguments:
 *
 *  \param u        gauge field ( Read )
 *  \param chi      color vector field ( Modify )
 *  \param j_decay  direction of decay ( Read )
 */

template<typename T>
void laplacian(const multi1d<LatticeColorMatrix>& u, 
	       T& chi, 
	       int j_decay)
{
  T psi = -1 * chi;

  /* hit with laplacian (Klein-Gordon with m=0) */
  klein_gord(u, psi, chi, 0, j_decay);
}

void laplacian(const multi1d<LatticeColorMatrix>& u, 
	       LatticeColorVector& chi, 
	       int j_decay)
{
  laplacian<LatticeColorVector>(u, chi, j_decay);
}

void laplacian(const multi1d<LatticeColorMatrix>& u, 
	       LatticePropagator& chi, 
	       int j_decay)
{
  laplacian<LatticePropagator>(u, chi, j_decay);
}


