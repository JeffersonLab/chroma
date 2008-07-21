// $Id: tdiractodr.cc,v 3.1 2008-07-21 02:33:28 edwards Exp $
/*! \file
 *  \brief Basis rotation matrix from Dirac to Degrand-Rossi (and reverse)
 */

#include "util/ferm/tdiractodr.h"

namespace Chroma 
{

  //! The Dirac to Degrand-Rossi spin transformation matrix
  /*!
   * \ingroup ferm
   *
   * Return the similarity transformation matrix from 
   * Euclidean Dirac to Euclidean Degrand-Rossi basis
   *
   * \returns The U in   Gamma_{Degrand-Rossi} = U Gamma_Dirac U^dag
   */

  SpinMatrixD TDiracToDRMat()
  {
    /*
     * The magic basis transf is found from
     *
     * NOTE: DR = Degrand-Rossi - the spin basis of QDP
     * NOTE: TDiract = "twisted" Dirac, in which sign of spatial cpts reversed
     *  psi_DR = U psi_TDirac
     *  psibar_DR Gamma_DR psi_DR = psibar_Dirac Gamma_TDirac psi_Dirac
     *
     * implies
     *  Gamma_DR = U Gamma_TDirac U^dag
     *
     * and the magic formula is
     *
     *   U = (1/sqrt(2)) | i*sigma_2    i*sigma_2 |
     *                   | i*sigma_2   -i*sigma_2 |
     *
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
    SpinMatrixD U = zero;
    RealD     foo = RealD(1) / sqrt(RealD(2));
    ComplexD  one = cmplx( foo,RealD(0));
    ComplexD mone = cmplx(-foo,RealD(0));

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

}  // end namespace Chroma
