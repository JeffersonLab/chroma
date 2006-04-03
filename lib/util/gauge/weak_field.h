// -*- C++ -*-
// $Id: weak_field.h,v 3.0 2006-04-03 04:59:13 edwards Exp $
/*! \file
 *  \brief Construct a weak field
 */

#ifndef __weakfield_h__
#define __weakfield_h__

namespace Chroma
{
  //! Construct a weak field
  /*!
   * \ingroup gauge
   *
   * Arguments:
   *
   *  \param u          Gauge field                   (Modify)
   */

  void weakField(multi1d<LatticeColorMatrix>& u);
}

#endif
