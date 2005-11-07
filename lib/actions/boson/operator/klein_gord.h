// -*- C++ -*-
// $Id: klein_gord.h,v 2.1 2005-11-07 21:26:37 edwards Exp $

/*! \file
 * \brief Klein-Gordon operator
 */

#ifndef KLEIN_GORD_INCLUDE
#define KLEIN_GORD_INCLUDE

namespace Chroma
{
  //! Klein-Gordon operator
  /*! @ingroup boson */
  void klein_gord(const multi1d<LatticeColorMatrix>& u, 
		  const LatticeColorVector& psi, 
		  LatticeColorVector& chi, 
		  const Real& mass_sq, int j_decay);

  //! Klein-Gordon operator
  /*! @ingroup boson */
  void klein_gord(const multi1d<LatticeColorMatrix>& u, 
		  const LatticeFermion& psi, 
		  LatticeFermion& chi, 
		  const Real& mass_sq, int j_decay);

  //! Klein-Gordon operator
  /*! @ingroup boson */
  void klein_gord(const multi1d<LatticeColorMatrix>& u, 
		  const LatticePropagator& psi, 
		  LatticePropagator& chi, 
		  const Real& mass_sq, int j_decay);
}


#endif
