// -*- C++ -*-
// $Id: constgauge.h,v 3.1 2006-08-25 23:46:37 edwards Exp $
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
