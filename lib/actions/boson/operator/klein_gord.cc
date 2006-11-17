// $Id: klein_gord.cc,v 3.1 2006-11-17 02:17:31 edwards Exp $
/*! \file
 *  \brief Apply Klein-Gordon operator
 */

#include "chromabase.h"
#include "actions/boson/operator/klein_gord.h"


namespace Chroma
{
  //! Compute the covariant Klein-Gordon operator
  /*!
   *  For 0 <= j_decay < Nd, the laplacian is only taken in the directions
   * other than j_decay.
   *
   *  Chi  :=  (mass_sq + 2*(Nd[-1]) * Psi -
   *           sum_{mu ne j_decay} [ U_mu(x) * Psi(x+mu) +
   *                                 U^dagger_mu(x-mu) * Psi(x-mu) ]  
   *
   *
   * Arguments:
   *
   *  \param u          Gauge field               (Read)
   *  \param psi        Lattice field             (Read)
   *  \param chi        Lattice field             (Write)
   *  \param mass_sq    square of the mass term   (Read)
   *  \param j_decay    'left out' direction      (Read) 
   */

  template<typename T>
  void klein_gord(const multi1d<LatticeColorMatrix>& u, 
		  const T& psi, 
		  T& chi, 
		  const Real& mass_sq, int j_decay)
  {
    Real ftmp;

    if( j_decay < Nd )
      ftmp = Real(2*Nd-2) + mass_sq;
    else
      ftmp = Real(2*Nd) + mass_sq;

    chi = psi * ftmp;

    for(int mu = 0; mu < Nd; ++mu )
      if( mu != j_decay )
      {
	chi -= u[mu]*shift(psi, FORWARD, mu) + shift(adj(u[mu])*psi, BACKWARD, mu);
      }
  }


  //! Compute the covariant Klein-Gordon operator on a color vector
  void klein_gord(const multi1d<LatticeColorMatrix>& u, 
		  const LatticeColorVector& psi, 
		  LatticeColorVector& chi, 
		  const Real& mass_sq, int j_decay)
  {
    klein_gord<LatticeColorVector>(u, psi, chi, mass_sq, j_decay);
  }

  //! Compute the covariant Klein-Gordon operator on a color vector
  void klein_gord(const multi1d<LatticeColorMatrix>& u, 
		  const LatticeFermion& psi, 
		  LatticeFermion& chi, 
		  const Real& mass_sq, int j_decay)
  {
    klein_gord<LatticeFermion>(u, psi, chi, mass_sq, j_decay);
  }

  //! Compute the covariant Klein-Gordon operator on a propagator
  void klein_gord(const multi1d<LatticeColorMatrix>& u, 
		  const LatticeStaggeredPropagator& psi, 
		  LatticeStaggeredPropagator& chi, 
		  const Real& mass_sq, int j_decay)
  {
    klein_gord<LatticeStaggeredPropagator>(u, psi, chi, mass_sq, j_decay);
  }

  //! Compute the covariant Klein-Gordon operator on a propagator
  void klein_gord(const multi1d<LatticeColorMatrix>& u, 
		  const LatticePropagator& psi, 
		  LatticePropagator& chi, 
		  const Real& mass_sq, int j_decay)
  {
    klein_gord<LatticePropagator>(u, psi, chi, mass_sq, j_decay);
  }

}
