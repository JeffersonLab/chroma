// $Id: laplacian.cc,v 1.3 2004-01-09 03:01:26 edwards Exp $
/*! \file
 *  \brief Laplacian smearing of a source
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
	       int j_decay,
	       int power)
{
  T psi;

  for (int p=0; p<power; p++) {
    psi = -1 * chi;

    /* hit with laplacian (Klein-Gordon with m=0) */
    klein_gord(u, psi, chi, 0, j_decay);
  }
}

void laplacian(const multi1d<LatticeColorMatrix>& u, 
	       LatticeColorVector& chi, 
	       int j_decay,
	       int power)
{
  laplacian<LatticeColorVector>(u, chi, j_decay, power);
}

void laplacian(const multi1d<LatticeColorMatrix>& u, 
	       LatticePropagator& chi, 
	       int j_decay,
	       int power)
{
  laplacian<LatticePropagator>(u, chi, j_decay, power);
}


