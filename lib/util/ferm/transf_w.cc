// $Id: transf_w.cc,v 1.1 2003-02-15 05:54:26 edwards Exp $
/*! \file
 *  \brief Insertion/Extraction utilities for vectors/fermions/propagators
 */

#include "qdp.h"
#include "szin.h"

using namespace QDP;

//! Convert (insert) a LatticeColorVector into a LatticeFermion
/*!
 * \param a      Source color vector
 * \param b      Destination fermion
 * \param spin_index   Spin index
 */
void CvToFerm(const LatticeColorVector& a, LatticeFermion& b, 
	      int spin_index)
{
  pokeSpin(b, a, spin_index);
}


//! Insert a LatticeFermion into a LatticeFermion
/*!
 * \param a      Source fermion
 * \param b      Destination propagator
 * \param color_index  Color index
 * \param spin_index   Spin index
 */
void FermToProp(const LatticeFermion& a, LatticePropagator& b, 
		int color_index, int spin_index)
{
  for(int j = 0; j < Ns; ++j)
  {
    LatticeColorMatrix bb = peekSpin(b, j, spin_index);
    LatticeColorVector aa = peekSpin(a, j);

    for(int i = 0; i < Nc; ++i)
      pokeColor(bb, peekColor(aa, i), i, color_index);

    pokeSpin(b, bb, j, spin_index);
  }
}


//! Extract a LatticeFermion from a LatticePropagator
/*!
 * \param a      Source propagator
 * \param b      Destination fermion
 * \param color_index  Color index
 * \param spin_index   Spin index
 */
void PropToFerm(const LatticePropagator& b, LatticeFermion& a, 
		int color_index, int spin_index)
{
  for(int j = 0; j < Ns; ++j)
  {
    LatticeColorMatrix bb = peekSpin(b, j, spin_index);
    LatticeColorVector aa = peekSpin(a, j);

    for(int i = 0; i < Nc; ++i)
      pokeColor(aa, peekColor(bb, i, color_index), i);

    pokeSpin(a, aa, j);
  }
}


