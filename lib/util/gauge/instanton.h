// -*- C++ -*-
/*! \file
 *  \brief Construct an instanton or anti-instanton configuration in singular gauge. 
*/

#ifndef __instanton_h__
#define __instanton_h__

namespace Chroma
{

  //! Instanton construction
  /*!
   * \ingroup gauge
   *
   * Arguments:
   *
   * \param u_inst          Gauge field                   (Write)
   * \param center          location of instanton center ( Read )
   * \param rho             size of instanton ( Read )
   * \param su2_index       SU(2) subgroup index ( Read )
   * \param sign            instanton (1) or anti-instanton ( Read ) 
   */

  void instanton(multi1d<LatticeColorMatrix>& u_inst, const multi1d<Real>& center, Real rho, int su2_index, int sign);
}

#endif
