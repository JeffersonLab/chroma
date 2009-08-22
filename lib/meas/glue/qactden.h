// -*- C++ -*-
// $Id: qactden.h,v 3.2 2009-08-22 20:05:49 edwards Exp $
/*! \file
 *  \brief Measure the lattice density of the lattice energy and the naive topological charge.
 */

#ifndef __qactden_h__
#define __qactden_h__

namespace Chroma 
{

  //! Measure the lattice density of the lattice energy and the naive topological charge.
  /*!
   * \ingroup glue
   *
   * \param lrqtop  topological charge density (Write)
   * \param lract   action to continuum instanton action density (Write) 
   * \param u       gauge field (Read)
   */

  void qactden(LatticeReal& lract, LatticeReal& lrqtop, const multi1d<LatticeColorMatrix>& u);

}  // end namespace Chroma

#endif
