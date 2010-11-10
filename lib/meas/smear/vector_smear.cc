// $Id: vector_smear.cc,v 3.2 2008-11-04 17:26:16 edwards Exp $
/*! \file
 *  \brief vector smearing of color vector
 */

#include "chromabase.h"
#include "meas/smear/vector_smear.h"
#include "util/ft/sftmom.h"

namespace Chroma 
{

  //! Do a vector smearing of a lattice colorvector field
  /*!
   * Arguments:
   *
   *  \param chi      color vector field ( Modify )
   *  \param vecs     vectors for the smearing ( Read )
   *  \param sigma    exponential smearing parameter ( Read )
   *  \param j_decay  direction of decay ( Read )
   */

  void vectorSmear(LatticeColorVector& chi, 
		   const MapObject<int,EVPair<LatticeColorVector> >& vecs,
		   const Real& sigma, const int& j_decay)
  {
    LatticeColorVector psi = zero;

    int num_vecs = vecs.size();

    SftMom phases(0, false, j_decay);
    int nt = phases.numSubsets();

    for (int n = 0 ; n < num_vecs ; ++n)
    {
      EVPair<LatticeColorVector> ev;
      vecs.get(n, ev);

      const multi1d<Real>& evals = ev.eigenValue.weights;
      LatticeColorVector cvec = ev.eigenVector;
      LatticeComplex tmp = localInnerProduct( cvec, chi );
		
      multi2d<DComplex> t_sum = phases.sft(tmp);

      if (nt != evals.size())
      {
	QDPIO::cerr << __func__ << ": number of evalues not decay extent\n";
	QDP_abort(1);
      }

      for (int t = 0 ; t < nt ; ++t)
      {
	psi[phases.getSet()[t]] += cvec * t_sum[0][t] * 
	  exp( Real(-1.0) * sigma * sigma / Real(4.0) * evals[t]);
      }
    }
	
    chi = psi;
  }


  //! Do a vector smearing of a color matrix
  /*!
   * Arguments:
   *
   *  \param chi      color matrix field ( Modify )
   *  \param vecs     vectors for the smearing ( Read )
   *  \param sigma    exponential smearing parameter ( Read )
   *  \param j_decay  direction of decay ( Read )
   */

  void vectorSmear(LatticeColorMatrix& chi, 
		   const MapObject<int,EVPair<LatticeColorVector> >& vecs,
		   const Real& sigma, const int& j_decay)
  {
    LatticeColorMatrix col_out = zero;

    // Smear a color matrix
    for (int cc = 0 ; cc < Nc ; ++cc)
    {
      LatticeColorVector vec_in = zero;

      // Extract column of color matrix
      for (int cr = 0 ; cr < Nc ; ++cr)
      {
	pokeColor(vec_in, peekColor(chi, cr, cc), cr);
      }

      // Smear that column
      vectorSmear(vec_in, vecs, sigma, j_decay);
		
      // Insert it back in
      for (int cr = 0 ; cr < Nc ; ++cr)
      {
	pokeColor(col_out, peekColor(vec_in, cr), cr, cc);
      }
    } // cc

    chi = col_out;
  }


  //! Do a vector smearing of a lattice fermion field
  /*!
   * Arguments:
   *
   *  \param chi      fermion field ( Modify )
   *  \param vecs     vectors for the smearing ( Read )
   *  \param sigma    exponential smearing parameter ( Read )
   *  \param j_decay  direction of decay ( Read )
   */

  void vectorSmear(LatticeFermion& chi, 
		   const MapObject<int,EVPair<LatticeColorVector> >& vecs,
		   const Real& sigma, const int& j_decay)
  {
    LatticeFermion psi = zero;

    for (int s = 0 ; s < Ns ; ++s)
    {
      LatticeColorVector chi_s = peekSpin(chi, s);
      vectorSmear(chi_s, vecs, sigma, j_decay);
		
      pokeSpin(psi, chi_s, s);
    }
	
    chi = psi;
  }

  //! Do a vector smearing of a lattice fermion field
  /*!
   * Arguments:
   *
   *  \param chi      fermion field ( Modify )
   *  \param vecs     vectors for the smearing ( Read )
   *  \param sigma    exponential smearing parameter ( Read )
   *  \param j_decay  direction of decay ( Read )
   */

  void vectorSmear(LatticeStaggeredPropagator& chi, 
		   const MapObject<int,EVPair<LatticeColorVector> >& vecs,
		   const Real& sigma, const int& j_decay)
  {
    LatticeStaggeredPropagator psi = zero;

    QDPIO::cerr << __func__ << ": not implemented\n";
    QDP_abort(1);

    chi = psi;
  }

  //! Do a vector smearing of a lattice fermion field
  /*!
   * Arguments:
   *
   *  \param chi      fermion field ( Modify )
   *  \param vecs     vectors for the smearing ( Read )
   *  \param sigma    exponential smearing parameter ( Read )
   *  \param j_decay  direction of decay ( Read )
   */

  void vectorSmear(LatticePropagator& chi, 
		   const MapObject<int,EVPair<LatticeColorVector> >& vecs,
		   const Real& sigma, const int& j_decay)
  {
    LatticePropagator psi = zero;

    for (int sc = 0 ; sc < Ns ; ++sc)
    {
      for (int sr = 0 ; sr < Ns ; ++sr)
      {
	LatticeColorMatrix col = peekSpin(chi, sr, sc);

	vectorSmear(col, vecs, sigma, j_decay);

	pokeSpin(psi, col, sr, sc);
      } // sr
    } // sc

    chi = psi;
  }

}  // end namespace Chroma
