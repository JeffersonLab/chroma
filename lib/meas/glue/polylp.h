// -*- C++ -*-
// $Id: polylp.h,v 3.1 2006-08-24 21:04:31 edwards Exp $
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

  void polylp(const multi1d<LatticeColorMatrixF3>& u, DComplex& poly_loop, int mu);

  void polylp(const multi1d<LatticeColorMatrixD3>& u, DComplex& poly_loop, int mu);

  //! Compute Polyakov loop
  /*!
   * \ingroup glue
   *
   * \param u          gauge field (Read)
   * \param poly_loop  Polyakov loop average (Write) 
   */

  void polylp(const multi1d<LatticeColorMatrixF3>& u, multi1d<DComplex>& poly_loop);
  void polylp(const multi1d<LatticeColorMatrixD3>& u, multi1d<DComplex>& poly_loop);

}  // end namespace Chroma

#endif
