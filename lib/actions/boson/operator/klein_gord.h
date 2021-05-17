// -*- C++ -*-

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
		  const LatticeColorMatrix& psi, 
		  LatticeColorMatrix& chi, 
		  const Real& mass_sq, int j_decay);

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
		  const LatticeStaggeredPropagator& psi, 
		  LatticeStaggeredPropagator& chi, 
		  const Real& mass_sq, int j_decay);

  //! Klein-Gordon operator
  /*! @ingroup boson */
  void klein_gord(const multi1d<LatticeColorMatrix>& u, 
		  const LatticePropagator& psi, 
		  LatticePropagator& chi, 
		  const Real& mass_sq, int j_decay);
}


#endif
