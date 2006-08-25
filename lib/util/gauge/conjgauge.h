// -*- C++ -*-
// $Id: conjgauge.h,v 3.1 2006-08-25 23:46:37 edwards Exp $
/*! \file
 *  \brief Take the complex conjugate of a gauge field
 */

#ifndef __conjgauge_h__
#define __conjgauge_h__

namespace Chroma
{

  //! Take complex conjugate of gauge field
  /*!
   * \ingroup gauge
   *
   * Arguments:
   *
   *  \param u          Gauge field                   (Modify)
   */

  void conjgauge(multi1d<LatticeColorMatrix>& u);
}

#endif
