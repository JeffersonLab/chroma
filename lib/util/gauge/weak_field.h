// -*- C++ -*-
// $Id: weak_field.h,v 3.1 2006-08-25 23:46:37 edwards Exp $
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
