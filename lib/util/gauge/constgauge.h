// -*- C++ -*-
// $Id: constgauge.h,v 3.0 2006-04-03 04:59:12 edwards Exp $
/*! \file
 *  \brief Take the complex conjugate of a gauge field
 */

#ifndef __constgauge_h__
#define __constgauge_h__

namespace Chroma
{


  //! Const diagonal gauge field
  /*!
   * \ingroup gauge
   *
   * Arguments:
   *
   *  \param u          Gauge field                   (Modify)
   *  \param theta      Angles                        (Read)
   */

  void constgauge(multi1d<LatticeColorMatrix>& u, const multi2d<Real> theta);
}

#endif
