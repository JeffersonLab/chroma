// $Id: paulitodr.cc,v 3.0 2006-04-03 04:59:11 edwards Exp $
/*! \file
 *  \brief Basis rotation matrix from Pauli-Schwinger (Euclidean Sakurai) to Degrand-Rossi (and reverse)
 */

#include "chromabase.h"
#include "util/ferm/paulitodr.h"

namespace Chroma {

//! The Pauli-Schwinger (Euclidean Sakurai) to Degrand-Rossi spin transformation matrix
/*!
 * \ingroup ferm
 *
 * Return the similarity transformation matrix from 
 * Euclidean Pauli-Schwinger to Euclidean Degrand-Rossi basis
 *
 * \returns The U in   Gamma_{Degrand-Rossi} = U Gamma_PS U^dag
 */

SpinMatrixD PauliToDRMat()
{
  /*
   * The magic basis transf is found from
   *
   * NOTE: DR = Degrand-Rossi - the spin basis of QDP
   *
   *  psi_DR = U psi_Pauli
   *  psibar_DR Gamma_DR psi_DR = psibar_Pauli Gamma_Pauli psi_Pauli
   *
   * implies
   *  Gamma_DR = U Gamma_Pauli U^dag
   *
   * and the magic formula is
   *
   *   U = (1/sqrt(2)) | i*sigma_2   -i*sigma_2 |
   *                   | i*sigma_2    i*sigma_2 |
   *
   *     = (1/sqrt(2)) |   0   1        0  -1   |
   *                   |  -1   0        1   0   |
   *                   |   0   1        0   1   |
   *                   |  -1   0       -1   0   |
   */
  /*
   * NOTE: I do not see some really short combination of 
   * QDP Gamma matrices that can make this beasty, 
   * so I'll just hardwire it...
   */
  SpinMatrixD U = zero;
  RealD     foo = RealD(1) / sqrt(RealD(2));
  ComplexD  one = cmplx( foo,RealD(0));
  ComplexD mone = cmplx(-foo,RealD(0));

  pokeSpin(U,  one, 0, 1);
  pokeSpin(U, mone, 0, 3);
  pokeSpin(U, mone, 1, 0);
  pokeSpin(U,  one, 1, 2);
  pokeSpin(U,  one, 2, 1);
  pokeSpin(U,  one, 2, 3);
  pokeSpin(U, mone, 3, 0);
  pokeSpin(U, mone, 3, 2);

  return U;
}

}  // end namespace Chroma
