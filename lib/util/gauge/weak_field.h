// -*- C++ -*-
// $Id: weak_field.h,v 2.1 2005-12-18 03:51:10 edwards Exp $
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
