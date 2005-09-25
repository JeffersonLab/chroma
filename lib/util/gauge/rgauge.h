// -*- C++ -*-
// $Id: rgauge.h,v 2.0 2005-09-25 21:04:44 edwards Exp $
/*! \file
 *  \brief Do a random gauge transformation on the u fields
 */

#ifndef __rgauge_h__
#define __rgauge_h__

namespace Chroma
{

  //! Do a random gauge transformation on the u fields
  /*!
   * \ingroup gauge
   *
   * Arguments:
   *
   *  \param u          Gauge field                   (Modify)
   *  \param g          Gauge transformation matrices (Write)
   */

  void rgauge(multi1d<LatticeColorMatrix>& u, LatticeColorMatrix& g);


  //! Do a random gauge transformation on the u fields
  /*!
   * \ingroup gauge
   *
   * Convenience function: does not return gauge transformation matrices
   *
   * Arguments:
   *
   *  \param u          Gauge field                   (Modify)
   */

  void rgauge(multi1d<LatticeColorMatrix>& u);

}

#endif
