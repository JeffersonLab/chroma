// $Id: transf.cc,v 1.2 2004-05-23 21:44:09 edwards Exp $
/*! \file
 *  \brief Insertion/Extraction utilities for vectors/fermions/propagators
 */

#include "chromabase.h"
#include "util/ferm/transf.h"

using namespace QDP;

//! Convert (insert) a LatticeColorVector into a LatticeFermion
/*!
 * \ingroup ferm
 *
 * \param a      Source color vector
 * \param b      Destination fermion
 * \param spin_index   Spin index
 */
void CvToFerm(const LatticeColorVectorF& a, LatticeFermionF& b, 
	      int spin_index)
{
  pokeSpin(b, a, spin_index);
}

//! Convert (insert) a LatticeColorVector into a LatticeFermion
/*!
 * \ingroup ferm
 *
 * \param a      Source color vector
 * \param b      Destination fermion
 * \param spin_index   Spin index
 */
void CvToFerm(const LatticeColorVectorD& a, LatticeFermionD& b, 
	      int spin_index)
{
  pokeSpin(b, a, spin_index);
}


//! Insert a LatticeFermion into a LatticeFermion
/*!
 * \ingroup ferm
 *
 * \param a      Source fermion
 * \param b      Destination propagator
 * \param color_index  Color index
 * \param spin_index   Spin index
 */
void FermToProp(const LatticeFermionF& a, LatticePropagatorF& b, 
		int color_index, int spin_index)
{
  for(int j = 0; j < Ns; ++j)
  {
    LatticeColorMatrixF bb = peekSpin(b, j, spin_index);
    LatticeColorVectorF aa = peekSpin(a, j);

    for(int i = 0; i < Nc; ++i)
      pokeColor(bb, peekColor(aa, i), i, color_index);

    pokeSpin(b, bb, j, spin_index);
  }
}

//! Insert a LatticeFermion into a LatticeFermion
/*!
 * \ingroup ferm
 *
 * \param a      Source fermion
 * \param b      Destination propagator
 * \param color_index  Color index
 * \param spin_index   Spin index
 */
void FermToProp(const LatticeFermionD& a, LatticePropagatorD& b, 
		int color_index, int spin_index)
{
  for(int j = 0; j < Ns; ++j)
  {
    LatticeColorMatrixD bb = peekSpin(b, j, spin_index);
    LatticeColorVectorD aa = peekSpin(a, j);

    for(int i = 0; i < Nc; ++i)
      pokeColor(bb, peekColor(aa, i), i, color_index);

    pokeSpin(b, bb, j, spin_index);
  }
}


//! Extract a LatticeFermion from a LatticePropagator
/*!
 * \ingroup ferm
 *
 * \param a      Source propagator
 * \param b      Destination fermion
 * \param color_index  Color index
 * \param spin_index   Spin index
 */
void PropToFerm(const LatticePropagatorF& b, LatticeFermionF& a, 
		int color_index, int spin_index)
{
  for(int j = 0; j < Ns; ++j)
  {
    LatticeColorMatrixF bb = peekSpin(b, j, spin_index);
    LatticeColorVectorF aa = peekSpin(a, j);

    for(int i = 0; i < Nc; ++i)
      pokeColor(aa, peekColor(bb, i, color_index), i);

    pokeSpin(a, aa, j);
  }
}

//! Extract a LatticeFermion from a LatticePropagator
/*!
 * \ingroup ferm
 *
 * \param a      Source propagator
 * \param b      Destination fermion
 * \param color_index  Color index
 * \param spin_index   Spin index
 */
void PropToFerm(const LatticePropagatorD& b, LatticeFermionD& a, 
		int color_index, int spin_index)
{
  for(int j = 0; j < Ns; ++j)
  {
    LatticeColorMatrixD bb = peekSpin(b, j, spin_index);
    LatticeColorVectorD aa = peekSpin(a, j);

    for(int i = 0; i < Nc; ++i)
      pokeColor(aa, peekColor(bb, i, color_index), i);

    pokeSpin(a, aa, j);
  }
}


