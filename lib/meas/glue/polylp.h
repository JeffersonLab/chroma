// -*- C++ -*-
/*! \file
 *  \brief Calculate the global normalized sum of the Polyakov loop
 */

#ifndef __polylp_h__
#define __polylp_h__

namespace Chroma 
{

  //! Compute Polyakov loop
  /*!
   * \ingroup glue
   *
   * \param u          gauge field (Read)
   * \param poly_loop  Polyakov loop average in direction mu (Write) 
   * \param mu         direction of Polyakov loop (Read)
   */

  template<typename Q> 
  void polylp(const multi1d<Q>& u, DComplex& poly_loop, int mu);

  //! Compute Polyakov loop
  /*!
   * \ingroup glue
   *
   * \param u          gauge field (Read)
   * \param poly_loop  Polyakov loop average (Write) 
   */

  template<typename Q>
  void polylp(const multi1d<Q>& u, multi1d<DComplex>& poly_loop);

}  // end namespace Chroma

#endif
