// $Id: diractodr.cc,v 1.1 2004-05-14 00:23:41 edwards Exp $
/*! \file
 *  \brief Basis rotation matrix from Dirac to Degrand-Rossi (and reverse)
 */

#include "chromabase.h"
#include "util/ferm/transf.h"

using namespace QDP;

//! The Dirac to Degrand-Rossi spin transformation matrix
/*!
 * \ingroup ferm
 *
 * Return the similarity transformation matrix from 
 * Euclidean Dirac to Euclidean Degrand-Rossi basis
 *
 * \returns The U in   Gamma_{Degrand-Rossi} = U Gamma_Dirac U^dag
 */

SpinMatrix DiracToDRMat()
{
  /*
   * The magic basis transf is found from
   *
   * NOTE: DR = Degrand-Rossi - the spin basis of QDP
   *
   *  psi_DR = U psi_Dirac
   *  psibar_DR Gamma_DR psi_DR = psibar_Dirac Gamma_Dirac psi_Dirac
   *
   * implies
   *  Gamma_DR = U Gamma_Dirac U^dag
   *
   * and the magic formula is
   *
   *   U = (1/sqrt(2)) | i*sigma_2    i*sigma_2 |
   *                   | i*sigma_2   -i*sigma_2 |
   *     = (1/sqrt(2)) |   0   1        0   1   |
   *                   |  -1   0       -1   0   |
   *                   |   0   1        0  -1   |
   *                   |  -1   0        1   0   |
   *
   *   U^dag = -U = U^transpose
   */
  /*
   * NOTE: I do not see some really short combination of 
   * QDP Gamma matrices that can make this beasty, 
   * so I'll just hardwire it...
   */
  SpinMatrix U = zero;
  Real     foo = 1 / sqrt(Real(2));
  Complex  one = cmplx( foo,Real(0));
  Complex mone = cmplx(-foo,Real(0));

  pokeSpin(U,  one, 0, 1);
  pokeSpin(U,  one, 0, 3);
  pokeSpin(U, mone, 1, 0);
  pokeSpin(U, mone, 1, 2);
  pokeSpin(U,  one, 2, 1);
  pokeSpin(U, mone, 2, 3);
  pokeSpin(U, mone, 3, 0);
  pokeSpin(U,  one, 3, 2);

  return U;
}
