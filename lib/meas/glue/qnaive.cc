// $Id: qnaive.cc,v 3.3 2008-04-24 14:17:14 edwards Exp $
/*! \file
 *  \brief Calculate the topological charge from the gluonic definition
 *
 *  Conventions are according to Bilson-Thompson et al., hep-lat/0203008
 *
 * Author: Christian Hagen
 */

#include "chromabase.h"
#include "meas/glue/qnaive.h"

namespace Chroma 
{

  //! Compute topological charge
  /*!
   * \ingroup glue
   *
   * \param u          gauge field (Read)
   * \param k5         improvement parameter (Read)
   * \param qtop       topological charge (Write) 
   */

  void qtop_naive(const multi1d<LatticeColorMatrix>& u, const Real k5, Double& qtop)
  {
    START_CODE();

    /* Local Variables */
    LatticeColorMatrix u_clov_1;
    LatticeColorMatrix u_clov_2;
    LatticeColorMatrix tmp_1;
    LatticeColorMatrix tmp_2;
    LatticeReal qtop_tmp;

    Real k1,k2,k3,k4,kk5;

//     Double tmp;

    if( Nd != 4 )
      QDP_error_exit("Nd for the topological charge has to be 4 but: ", Nd);

    k1 = 19.0/9.0 - 55.0 * k5;
    k1 *= 2.0;
    k2 = 1.0/36.0 - 16.0 * k5;
    k2 *= 2.0;
    k3 = 64.0 * k5 - 32.0/45.0;
    k4 = 1.0/15.0 - 6.0 * k5;
    kk5 = k5;
    kk5 *= 2.0;

    qtop = 0;

    int mu1 = 0;
    for(int nu1=1; nu1<Nd; ++nu1)
    {

      if( toBool(k1 != 0) ) {

//       QDPIO::cout << "k1 = " << k1 << endl;

      /* First "plus-plus" 1x1 */
      tmp_1 = u[mu1] * shift(u[nu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(adj(u[mu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 = k1 * tmp_1;

      /* First "plus-minus" 1x1 */
      tmp_1 = u[mu1] * shift(shift(adj(u[nu1]), FORWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[mu1]), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 -= k1 * tmp_1;

      /* First "minus-minus" 1x1 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 += k1 * tmp_1;

      /* First "minus-plus" 1x1 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(u[nu1], BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(u[mu1], BACKWARD, mu1), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 -= k1 * tmp_1;

      }

      if( toBool(k2!=0) ) {

//       QDPIO::cout << "k2 = " << k2 << endl;

      /* First "plus-plus" 2x2 */
      tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], FORWARD, mu1), FORWARD, mu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu1], FORWARD, mu1), FORWARD, mu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), FORWARD, mu1), FORWARD, nu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(adj(u[mu1]), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 += k2 * tmp_1;

      /* First "plus-minus" 2x2 */
      tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu1]), FORWARD, mu1), FORWARD, mu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu1]), FORWARD, mu1), FORWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), FORWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(adj(u[mu1]), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 -= k2 * tmp_1;

      /* First "minus-minus" 2x2 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 += k2 * tmp_1;

      /* First "minus-plus" 2x2 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, mu1), BACKWARD, mu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu1], BACKWARD, mu1), BACKWARD, mu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), FORWARD, nu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[mu1], BACKWARD, mu1), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 -= k2 * tmp_1;

      }

      if( toBool(k3!=0) ) {

//       QDPIO::cout << "k3 = " << k3 << endl;

      /* First "plus-plus" 2x1 */
      tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], FORWARD, mu1), FORWARD, mu1);
      tmp_1 = tmp_2 * shift(shift(adj(u[mu1]), FORWARD, mu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[mu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 += k3 * tmp_1;

      /* First "plus-minus" 2x1 */
      tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu1]), FORWARD, mu1), FORWARD, mu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(adj(u[mu1]), FORWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[mu1]), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 -= k3 * tmp_1;

      /* First "minus-minus" 2x1 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu1]),BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[mu1],BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 += k3 * tmp_1;

      /* First "minus-plus" 2x1 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, mu1), BACKWARD, mu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(u[mu1], BACKWARD, mu1), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 -= k3 * tmp_1;


      /* First "plus-plus" 1x2 */
      tmp_1 = u[mu1] * shift(u[nu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], FORWARD, mu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(adj(u[mu1]), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 += k3 * tmp_1;

      /* First "plus-minus" 1x2 */
      tmp_1 = u[mu1] * shift(shift(adj(u[nu1]), FORWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu1]), FORWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(adj(u[mu1]), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 -= k3 * tmp_1;

      /* First "minus-minus" 1x2 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 += k3 * tmp_1;

      /* First "minus-plus" 1x2 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(u[nu1], BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, mu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[mu1], BACKWARD, mu1), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 -= k3 * tmp_1;

      }

      if( toBool(k4!=0) ) {

//       QDPIO::cout << "k4 = " << k4 << endl;

      /* First "plus-plus" 3x1 */
      tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(u[mu1], FORWARD, mu1), FORWARD, mu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu1], FORWARD, mu1), FORWARD, mu1), FORWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), FORWARD, mu1), FORWARD, mu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(adj(u[mu1]), FORWARD, mu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[mu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 += k4 * tmp_1;

      /* First "plus-minus" 3x1 */
      tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(u[mu1], FORWARD, mu1), FORWARD, mu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu1]), FORWARD, mu1), FORWARD, mu1), FORWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), FORWARD, mu1), FORWARD, mu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(adj(u[mu1]), FORWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[mu1]), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 -= k4 * tmp_1;

      /* First "minus-minus" 3x1 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 += k4 * tmp_1;

      /* First "minus-plus" 3x1 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(u[mu1], BACKWARD, mu1), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 -= k4 * tmp_1;


      /* First "plus-plus" 1x3 */
      tmp_1 = u[mu1] * shift(u[nu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], FORWARD, mu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu1], FORWARD, mu1), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), FORWARD, nu1), FORWARD, nu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(adj(u[nu1]), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 += k4 * tmp_1;

      /* First "plus-minus" 1x3 */
      tmp_1 = u[mu1] * shift(shift(adj(u[nu1]), FORWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu1]), FORWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu1]), FORWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 -= k4 * tmp_1;

      /* First "minus-minus" 1x3 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu1]),BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 += k4 * tmp_1;

      /* First "minus-plus" 1x3 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(u[nu1], BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, mu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu1], BACKWARD, mu1), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu1], BACKWARD, mu1), FORWARD, nu1), FORWARD, nu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(adj(u[nu1]), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 -= k4 * tmp_1;

      }

      if( toBool(kk5!=0) ) {

//       QDPIO::cout << "k5 = " << kk5 << endl;

      /* First "plus-plus" 3x3 */
      tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(u[mu1], FORWARD, mu1), FORWARD, mu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu1], FORWARD, mu1), FORWARD, mu1), FORWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[nu1], FORWARD, mu1), FORWARD, mu1), FORWARD, mu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(shift(u[nu1], FORWARD, mu1), FORWARD, mu1), FORWARD, mu1), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(shift(adj(u[mu1]), FORWARD, mu1), FORWARD, mu1), FORWARD, nu1), FORWARD, nu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[mu1]), FORWARD, mu1), FORWARD, nu1), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), FORWARD, nu1), FORWARD, nu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(adj(u[nu1]), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 += kk5 * tmp_1;

      /* First "plus-minus" 3x3 */
      tmp_1 = u[mu1] * shift(u[mu1], FORWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(u[mu1], FORWARD, mu1), FORWARD, mu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu1]), FORWARD, mu1), FORWARD, mu1), FORWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(shift(adj(u[nu1]), FORWARD, mu1), FORWARD, mu1), FORWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(shift(shift(adj(u[nu1]), FORWARD, mu1), FORWARD, mu1), FORWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(shift(adj(u[mu1]), FORWARD, mu1), FORWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[mu1]), FORWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 -= kk5 * tmp_1;

      /* First "minus-minus" 3x3 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(shift(shift(adj(u[nu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1), BACKWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(u[nu1], BACKWARD, nu1), BACKWARD, nu1);
      tmp_1 = tmp_2 * shift(u[nu1], BACKWARD, nu1);
      u_clov_1 += kk5 * tmp_1;

      /* First "minus-plus" 3x3 */
      tmp_1 = shift(adj(u[mu1]), BACKWARD, mu1) * shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu1]), BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[nu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(shift(u[nu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), BACKWARD, mu1), FORWARD, nu1), FORWARD, nu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(shift(shift(shift(u[mu1], BACKWARD, mu1), BACKWARD, mu1), FORWARD, nu1), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu1], BACKWARD, mu1), FORWARD, nu1), FORWARD, nu1), FORWARD, nu1);
      tmp_1 = tmp_2 * shift(shift(adj(u[nu1]), FORWARD, nu1), FORWARD, nu1);
      tmp_2 = tmp_1 * shift(adj(u[nu1]), FORWARD, nu1);
      tmp_1 = tmp_2 * adj(u[nu1]);
      u_clov_1 -= kk5 * tmp_1;

      }

      int mu2 = nu1%3;
      mu2++;
      int nu2 = mu2%3;
      nu2++;


      if( toBool(k1!=0) ) {

      /* Second "plus-plus" 1x1 */
      tmp_1 = u[mu2] * shift(u[nu2], FORWARD, mu2);
      tmp_2 = tmp_1 * shift(adj(u[mu2]), FORWARD, nu2);
      tmp_1 = tmp_2 * adj(u[nu2]);
      u_clov_2 = k1 * tmp_1;

      /* Second "plus-minus" 1x1 */
      tmp_1 = u[mu2] * shift(shift(adj(u[nu2]), FORWARD, mu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(adj(u[mu2]), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(u[nu2], BACKWARD, nu2);
      u_clov_2 -= k1 * tmp_1;

      /* Second "minus-minus" 1x1 */
      tmp_1 = shift(adj(u[mu2]), BACKWARD, mu2) * shift(shift(adj(u[nu2]), BACKWARD, mu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(u[mu2], BACKWARD, mu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(u[nu2], BACKWARD, nu2);
      u_clov_2 += k1 * tmp_1;

      /* Second "minus-plus" 1x1 */
      tmp_1 = shift(adj(u[mu2]), BACKWARD, mu2) * shift(u[nu2], BACKWARD, mu2);
      tmp_2 = tmp_1 * shift(shift(u[mu2], BACKWARD, mu2), FORWARD, nu2);
      tmp_1 = tmp_2 * adj(u[nu2]);
      u_clov_2 -= k1 * tmp_1;

      }

      if( toBool(k2!=0) ) {

      /* Second "plus-plus" 2x2 */
      tmp_1 = u[mu2] * shift(u[mu2], FORWARD, mu2);
      tmp_2 = tmp_1 * shift(shift(u[nu2], FORWARD, mu2), FORWARD, mu2);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu2], FORWARD, mu2), FORWARD, mu2), FORWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu2]), FORWARD, mu2), FORWARD, nu2), FORWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(adj(u[mu2]), FORWARD, nu2), FORWARD, nu2);
      tmp_2 = tmp_1 * shift(adj(u[nu2]), FORWARD, nu2);
      tmp_1 = tmp_2 * adj(u[nu2]);
      u_clov_2 += k2 * tmp_1;

      /* Second "plus-minus" 2x2 */
      tmp_1 = u[mu2] * shift(u[mu2], FORWARD, mu2);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu2]), FORWARD, mu2), FORWARD, mu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu2]), FORWARD, mu2), FORWARD, mu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu2]), FORWARD, mu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(adj(u[mu2]), BACKWARD, nu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(u[nu2], BACKWARD, nu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(u[nu2], BACKWARD, nu2);
      u_clov_2 -= k2 * tmp_1;

      /* Second "minus-minus" 2x2 */
      tmp_1 = shift(adj(u[mu2]), BACKWARD, mu2) * shift(shift(adj(u[mu2]), BACKWARD, mu2), BACKWARD, mu2);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu2]), BACKWARD, mu2), BACKWARD, mu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu2]), BACKWARD, mu2), BACKWARD, mu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu2], BACKWARD, mu2), BACKWARD, mu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(shift(u[mu2], BACKWARD, mu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(u[nu2], BACKWARD, nu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(u[nu2], BACKWARD, nu2);
      u_clov_2 += k2 * tmp_1;

      /* Second "minus-plus" 2x2 */
      tmp_1 = shift(adj(u[mu2]), BACKWARD, mu2) * shift(shift(adj(u[mu2]), BACKWARD, mu2), BACKWARD, mu2);
      tmp_2 = tmp_1 * shift(shift(u[nu2], BACKWARD, mu2), BACKWARD, mu2);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu2], BACKWARD, mu2), BACKWARD, mu2), FORWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu2], BACKWARD, mu2), BACKWARD, mu2), FORWARD, nu2), FORWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(shift(u[mu2], BACKWARD, mu2), FORWARD, nu2), FORWARD, nu2);
      tmp_2 = tmp_1 * shift(adj(u[nu2]), FORWARD, nu2);
      tmp_1 = tmp_2 * adj(u[nu2]);
      u_clov_2 -= k2 * tmp_1;

      }

      if( toBool(k3!=0) ) {

      /* First "plus-plus" 2x1 */
      tmp_1 = u[mu2] * shift(u[mu2], FORWARD, mu2);
      tmp_2 = tmp_1 * shift(shift(u[nu2], FORWARD, mu2), FORWARD, mu2);
      tmp_1 = tmp_2 * shift(shift(adj(u[mu2]), FORWARD, mu2), FORWARD, nu2);
      tmp_2 = tmp_1 * shift(adj(u[mu2]), FORWARD, nu2);
      tmp_1 = tmp_2 * adj(u[nu2]);
      u_clov_2 += k3 * tmp_1;

      /* First "plus-minus" 2x1 */
      tmp_1 = u[mu2] * shift(u[mu2], FORWARD, mu2);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu2]), FORWARD, mu2), FORWARD, mu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(adj(u[mu2]), FORWARD, mu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(adj(u[mu2]), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(u[nu2], BACKWARD, nu2);
      u_clov_2 -= k3 * tmp_1;

      /* First "minus-minus" 2x1 */
      tmp_1 = shift(adj(u[mu2]), BACKWARD, mu2) * shift(shift(adj(u[mu2]), BACKWARD, mu2), BACKWARD, mu2);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu2]),BACKWARD, mu2), BACKWARD, mu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(shift(u[mu2],BACKWARD, mu2), BACKWARD, mu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(u[mu2], BACKWARD, mu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(u[nu2], BACKWARD, nu2);
      u_clov_2 += k3 * tmp_1;

      /* First "minus-plus" 2x1 */
      tmp_1 = shift(adj(u[mu2]), BACKWARD, mu2) * shift(shift(adj(u[mu2]), BACKWARD, mu2), BACKWARD, mu2);
      tmp_2 = tmp_1 * shift(shift(u[nu2], BACKWARD, mu2), BACKWARD, mu2);
      tmp_1 = tmp_2 * shift(shift(shift(u[mu2], BACKWARD, mu2), BACKWARD, mu2), FORWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(u[mu2], BACKWARD, mu2), FORWARD, nu2);
      tmp_1 = tmp_2 * adj(u[nu2]);
      u_clov_2 -= k3 * tmp_1;


      /* First "plus-plus" 1x2 */
      tmp_1 = u[mu2] * shift(u[nu2], FORWARD, mu2);
      tmp_2 = tmp_1 * shift(shift(u[nu2], FORWARD, mu2), FORWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(adj(u[mu2]), FORWARD, nu2), FORWARD, nu2);
      tmp_2 = tmp_1 * shift(adj(u[nu2]), FORWARD, nu2);
      tmp_1 = tmp_2 * adj(u[nu2]);
      u_clov_2 += k3 * tmp_1;

      /* First "plus-minus" 1x2 */
      tmp_1 = u[mu2] * shift(shift(adj(u[nu2]), FORWARD, mu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu2]), FORWARD, mu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(adj(u[mu2]), BACKWARD, nu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(u[nu2], BACKWARD, nu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(u[nu2], BACKWARD, nu2);
      u_clov_2 -= k3 * tmp_1;

      /* First "minus-minus" 1x2 */
      tmp_1 = shift(adj(u[mu2]), BACKWARD, mu2) * shift(shift(adj(u[nu2]), BACKWARD, mu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu2]), BACKWARD, mu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(shift(u[mu2], BACKWARD, mu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(u[nu2], BACKWARD, nu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(u[nu2], BACKWARD, nu2);
      u_clov_2 += k3 * tmp_1;

      /* First "minus-plus" 1x2 */
      tmp_1 = shift(adj(u[mu2]), BACKWARD, mu2) * shift(u[nu2], BACKWARD, mu2);
      tmp_2 = tmp_1 * shift(shift(u[nu2], BACKWARD, mu2), FORWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(shift(u[mu2], BACKWARD, mu2), FORWARD, nu2), FORWARD, nu2);
      tmp_2 = tmp_1 * shift(adj(u[nu2]), FORWARD, nu2);
      tmp_1 = tmp_2 * adj(u[nu2]);
      u_clov_2 -= k3 * tmp_1;

      }

      if( toBool(k4!=0) ) {

      /* Second "plus-plus" 3x1 */
      tmp_1 = u[mu2] * shift(u[mu2], FORWARD, mu2);
      tmp_2 = tmp_1 * shift(shift(u[mu2], FORWARD, mu2), FORWARD, mu2);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu2], FORWARD, mu2), FORWARD, mu2), FORWARD, mu2);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu2]), FORWARD, mu2), FORWARD, mu2), FORWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(adj(u[mu2]), FORWARD, mu2), FORWARD, nu2);
      tmp_2 = tmp_1 * shift(adj(u[mu2]), FORWARD, nu2);
      tmp_1 = tmp_2 * adj(u[nu2]);
      u_clov_2 += k4 * tmp_1;

      /* Second "plus-minus" 3x1 */
      tmp_1 = u[mu2] * shift(u[mu2], FORWARD, mu2);
      tmp_2 = tmp_1 * shift(shift(u[mu2], FORWARD, mu2), FORWARD, mu2);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu2]), FORWARD, mu2), FORWARD, mu2), FORWARD, mu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu2]), FORWARD, mu2), FORWARD, mu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(adj(u[mu2]), FORWARD, mu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(adj(u[mu2]), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(u[nu2], BACKWARD, nu2);
      u_clov_2 -= k4 * tmp_1;

      /* Second "minus-minus" 3x1 */
      tmp_1 = shift(adj(u[mu2]), BACKWARD, mu2) * shift(shift(adj(u[mu2]), BACKWARD, mu2), BACKWARD, mu2);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu2]), BACKWARD, mu2), BACKWARD, mu2), BACKWARD, mu2);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu2]), BACKWARD, mu2), BACKWARD, mu2), BACKWARD, mu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu2], BACKWARD, mu2), BACKWARD, mu2), BACKWARD, mu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(shift(u[mu2], BACKWARD, mu2), BACKWARD, mu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(u[mu2], BACKWARD, mu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(u[nu2], BACKWARD, nu2);
      u_clov_2 += k4 * tmp_1;

      /* Second "minus-plus" 3x1 */
      tmp_1 = shift(adj(u[mu2]), BACKWARD, mu2) * shift(shift(adj(u[mu2]), BACKWARD, mu2), BACKWARD, mu2);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu2]), BACKWARD, mu2), BACKWARD, mu2), BACKWARD, mu2);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu2], BACKWARD, mu2), BACKWARD, mu2), BACKWARD, mu2);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu2], BACKWARD, mu2), BACKWARD, mu2), BACKWARD, mu2), FORWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(shift(u[mu2], BACKWARD, mu2), BACKWARD, mu2), FORWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(u[mu2], BACKWARD, mu2), FORWARD, nu2);
      tmp_1 = tmp_2 * adj(u[nu2]);
      u_clov_2 -= k4 * tmp_1;


      /* Second "plus-plus" 1x3 */
      tmp_1 = u[mu2] * shift(u[nu2], FORWARD, mu2);
      tmp_2 = tmp_1 * shift(shift(u[nu2], FORWARD, mu2), FORWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu2], FORWARD, mu2), FORWARD, nu2), FORWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu2]), FORWARD, nu2), FORWARD, nu2), FORWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(adj(u[nu2]), FORWARD, nu2), FORWARD, nu2);
      tmp_2 = tmp_1 * shift(adj(u[nu2]), FORWARD, nu2);
      tmp_1 = tmp_2 * adj(u[nu2]);
      u_clov_2 += k4 * tmp_1;

      /* Second "plus-minus" 1x3 */
      tmp_1 = u[mu2] * shift(shift(adj(u[nu2]), FORWARD, mu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu2]), FORWARD, mu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu2]), FORWARD, mu2), BACKWARD, nu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu2]), BACKWARD, nu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu2], BACKWARD, nu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(u[nu2], BACKWARD, nu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(u[nu2], BACKWARD, nu2);
      u_clov_2 -= k4 * tmp_1;

      /* Second "minus-minus" 1x3 */
      tmp_1 = shift(adj(u[mu2]), BACKWARD, mu2) * shift(shift(adj(u[nu2]), BACKWARD, mu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[nu2]), BACKWARD, mu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu2]),BACKWARD, mu2), BACKWARD, nu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu2], BACKWARD, mu2), BACKWARD, nu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu2], BACKWARD, nu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(u[nu2], BACKWARD, nu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(u[nu2], BACKWARD, nu2);
      u_clov_2 += k4 * tmp_1;

      /* Second "minus-plus" 1x3 */
      tmp_1 = shift(adj(u[mu2]), BACKWARD, mu2) * shift(u[nu2], BACKWARD, mu2);
      tmp_2 = tmp_1 * shift(shift(u[nu2], BACKWARD, mu2), FORWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu2], BACKWARD, mu2), FORWARD, nu2), FORWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu2], BACKWARD, mu2), FORWARD, nu2), FORWARD, nu2), FORWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(adj(u[nu2]), FORWARD, nu2), FORWARD, nu2);
      tmp_2 = tmp_1 * shift(adj(u[nu2]), FORWARD, nu2);
      tmp_1 = tmp_2 * adj(u[nu2]);
      u_clov_2 -= k4 * tmp_1;

      }

      if( toBool(kk5!=0) ) {

      /* First "plus-plus" 3x3 */
      tmp_1 = u[mu2] * shift(u[mu2], FORWARD, mu2);
      tmp_2 = tmp_1 * shift(shift(u[mu2], FORWARD, mu2), FORWARD, mu2);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu2], FORWARD, mu2), FORWARD, mu2), FORWARD, mu2);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[nu2], FORWARD, mu2), FORWARD, mu2), FORWARD, mu2), FORWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(shift(shift(shift(u[nu2], FORWARD, mu2), FORWARD, mu2), FORWARD, mu2), FORWARD, nu2), FORWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(shift(shift(shift(adj(u[mu2]), FORWARD, mu2), FORWARD, mu2), FORWARD, nu2), FORWARD, nu2), FORWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[mu2]), FORWARD, mu2), FORWARD, nu2), FORWARD, nu2), FORWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu2]), FORWARD, nu2), FORWARD, nu2), FORWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(adj(u[nu2]), FORWARD, nu2), FORWARD, nu2);
      tmp_2 = tmp_1 * shift(adj(u[nu2]), FORWARD, nu2);
      tmp_1 = tmp_2 * adj(u[nu2]);
      u_clov_2 += kk5 * tmp_1;

      /* First "plus-minus" 3x3 */
      tmp_1 = u[mu2] * shift(u[mu2], FORWARD, mu2);
      tmp_2 = tmp_1 * shift(shift(u[mu2], FORWARD, mu2), FORWARD, mu2);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu2]), FORWARD, mu2), FORWARD, mu2), FORWARD, mu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(shift(shift(shift(adj(u[nu2]), FORWARD, mu2), FORWARD, mu2), FORWARD, mu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(shift(shift(shift(shift(adj(u[nu2]), FORWARD, mu2), FORWARD, mu2), FORWARD, mu2), BACKWARD, nu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(shift(shift(shift(adj(u[mu2]), FORWARD, mu2), FORWARD, mu2), BACKWARD, nu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[mu2]), FORWARD, mu2), BACKWARD, nu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu2]), BACKWARD, nu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu2], BACKWARD, nu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(u[nu2], BACKWARD, nu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(u[nu2], BACKWARD, nu2);
      u_clov_2 -= kk5 * tmp_1;

      /* First "minus-minus" 3x3 */
      tmp_1 = shift(adj(u[mu2]), BACKWARD, mu2) * shift(shift(adj(u[mu2]), BACKWARD, mu2), BACKWARD, mu2);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu2]), BACKWARD, mu2), BACKWARD, mu2), BACKWARD, mu2);
      tmp_1 = tmp_2 * shift(shift(shift(shift(adj(u[nu2]), BACKWARD, mu2), BACKWARD, mu2), BACKWARD, mu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(shift(shift(shift(adj(u[nu2]), BACKWARD, mu2), BACKWARD, mu2), BACKWARD, mu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(shift(shift(shift(shift(adj(u[nu2]), BACKWARD, mu2), BACKWARD, mu2), BACKWARD, mu2), BACKWARD, nu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(shift(shift(shift(shift(u[mu2], BACKWARD, mu2), BACKWARD, mu2), BACKWARD, mu2), BACKWARD, nu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(shift(shift(shift(u[mu2], BACKWARD, mu2), BACKWARD, mu2), BACKWARD, nu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu2], BACKWARD, mu2), BACKWARD, nu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu2], BACKWARD, nu2), BACKWARD, nu2), BACKWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(u[nu2], BACKWARD, nu2), BACKWARD, nu2);
      tmp_1 = tmp_2 * shift(u[nu2], BACKWARD, nu2);
      u_clov_2 += kk5 * tmp_1;

      /* First "minus-plus" 3x3 */
      tmp_1 = shift(adj(u[mu2]), BACKWARD, mu2) * shift(shift(adj(u[mu2]), BACKWARD, mu2), BACKWARD, mu2);
      tmp_2 = tmp_1 * shift(shift(shift(adj(u[mu2]), BACKWARD, mu2), BACKWARD, mu2), BACKWARD, mu2);
      tmp_1 = tmp_2 * shift(shift(shift(u[nu2], BACKWARD, mu2), BACKWARD, mu2), BACKWARD, mu2);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[nu2], BACKWARD, mu2), BACKWARD, mu2), BACKWARD, mu2), FORWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(shift(shift(shift(u[nu2], BACKWARD, mu2), BACKWARD, mu2), BACKWARD, mu2), FORWARD, nu2), FORWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(shift(shift(shift(shift(u[mu2], BACKWARD, mu2), BACKWARD, mu2), BACKWARD, mu2), FORWARD, nu2), FORWARD, nu2), FORWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(shift(shift(shift(u[mu2], BACKWARD, mu2), BACKWARD, mu2), FORWARD, nu2), FORWARD, nu2), FORWARD, nu2);
      tmp_2 = tmp_1 * shift(shift(shift(shift(u[mu2], BACKWARD, mu2), FORWARD, nu2), FORWARD, nu2), FORWARD, nu2);
      tmp_1 = tmp_2 * shift(shift(adj(u[nu2]), FORWARD, nu2), FORWARD, nu2);
      tmp_2 = tmp_1 * shift(adj(u[nu2]), FORWARD, nu2);
      tmp_1 = tmp_2 * adj(u[nu2]);
      u_clov_2 -= kk5 * tmp_1;

      }

      /* Now comes the contribution to the topological charge */      

      tmp_1=1;

      tmp_2 = adj(u_clov_1);
      u_clov_1 -= tmp_2;
      tmp_2 = tmp_1 * trace(u_clov_1)/Nc;
      u_clov_1 -= tmp_2;
      tmp_2 = adj(u_clov_2);
      u_clov_2 -= tmp_2;
      tmp_2 = tmp_1 * trace(u_clov_2)/Nc;
      u_clov_2 -= tmp_2;
      tmp_2 = u_clov_1 * u_clov_2;

      qtop_tmp = real(trace(tmp_2));
      qtop -= sum(qtop_tmp);

    }

    /* Topological charge */
    qtop /= ( 16*16*twopi*twopi );
    QDPIO::cout << "qtop = " << qtop << endl;

    END_CODE();
  }

}  // end namespace Chroma
