// $Id: qnaive.h,v 3.2 2008-04-24 14:17:14 edwards Exp $
/*! \file
 *  \brief Calculate the topological charge from the gluonic definition
 *
 *  Conventions are according to Bilson-Thompson et al., hep-lat/0203008
 *
 * Author: Christian Hagen
 */

#ifndef __qnaive_h__
#define __qnaive_h__

namespace Chroma 
{

  //! Compute top charge
  /*!
   * \ingroup glue
   *
   * \param u          gauge field (Read)
   * \param k5         improvement parameter (Read)
   * \param qtop       topological charge (Write) 
   */

  void qtop_naive(const multi1d<LatticeColorMatrix>& u, const Real k5, Double& qtop);

}  // end namespace Chroma

#endif
