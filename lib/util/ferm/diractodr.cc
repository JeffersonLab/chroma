// $Id: diractodr.cc,v 3.0 2006-04-03 04:59:11 edwards Exp $
/*! \file
 *  \brief Basis rotation matrix from Dirac to Degrand-Rossi (and reverse)
 */

#include "chromabase.h"
#include "util/ferm/diractodr.h"


namespace Chroma {

//! The Dirac to Degrand-Rossi spin transformation matrix (and reverse)
/*!
 * \ingroup ferm
 *
 * Return the similarity transformation matrix from 
 * Euclidean Dirac to Euclidean Degrand-Rossi basis (or reverse)
 *
 * \returns the U such that  Gamma_{Dirac} = U^\dagger Gamma_{DeGrand-Rossi} U
 *            (or such that  Gamma_{DeGrand-Rossi} = U Gamma_{Dirac} U^\dagger )
 */

SpinMatrixD DiracToDRMat()
{
  /*
   * Following is the definition of unitary matrix which transforms
   * the basis of gamma matrices from DeGrand-Rossi (DR) to Dirac-Pauli (DP).
   * (Note: the default gamma matrices in chroma/qdp/szin is DeGrand-Rossi)
   *
   * \gamma_\mu(DP)       = U^\dagger \gamma_\mu(DR) U
   * \psi(DP)             = U^\dagger \psi(DR)
   * quark_propagator(DP) = U^\dagger quark_propagator(DR) U
   *
   * 
   *        i     /  s2   -s2 \        1     /  0   1   0  -1 \
   * U = ------- |             | =  ------- |  -1   0   1   0  |
   *     sqrt(2) |             |    sqrt(2) |   0   1   0   1  |
   *              \  s2    s2 /              \ -1   0  -1   0 /
   *
   * where s2 is Pauli sigma_2 matrix.
   * 
   * U_inverse = U^\dagger
   *
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
