// $Id: gaus_smear.cc,v 1.3 2003-02-16 04:14:37 edwards Exp $
/*! \file
 *  \brief Gaussian smearing of color vector
 */

#include "chromabase.h"
#include "actions/dirac/klein_gord.h"

using namespace QDP;

//! Do a covariant Gaussian smearing of a color vector field
/*!
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
	       const Real& width, int ItrGaus, int j_decay)
{
  LatticeColorVector psi;

  Real ftmp = - (width*width) / Real(4*ItrGaus);
  /* The Klein-Gordon operator is (Lapl + mass_sq), where Lapl = -d^2/dx^2.. */
  /* We want (1 + ftmp * Lapl ) = (Lapl + 1/ftmp)*ftmp */
  Real ftmpi = Real(1) / ftmp;
  
  for(int n = 0; n < ItrGaus; ++n)
  {
    psi = chi * ftmp;
    klein_gord(u, psi, chi, ftmpi, j_decay);
  }
}
