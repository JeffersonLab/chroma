/*! \file
 * \brief Contraction operators for two quarks
 */

#include "util/ferm/twoquark_contract_ops.h"

namespace Chroma
{
  //----------------------------------------------------------------------------
  //! Take transpose of a matrix in (explicit) spin space
  void transpose(multi2d<LatticeColorVector>& dist_rep, 
		 const multi2d<LatticeColorVector>& prop_rep)
  { 
    // Sanity check
    if (prop_rep.nrows() != Ns || prop_rep.ncols() != Ns)
    {
      QDPIO::cerr << __func__ << ": invalid size of prop_rep\n";
      QDP_abort(1);
    }

    dist_rep.resize(Ns,Ns);

    // spin row
    for(int sl = 0; sl < Ns; ++sl)
    {
      // spin col
      for(int sr = 0; sr < Ns; ++sr)
      {
	dist_rep(sl,sr) = prop_rep(sr,sl);    // WHAT ABOUT adj() ??
      } // sr
    } // sl
  } // void


  //----------------------------------------------------------------------------
  // Use gamma_5 hermiticity on a prop
  void gamma5Herm(multi2d<LatticeColorVector>& prop_time)
  {
    // Convert Gamma(15) into a representation of a spin matrix
    // NOTE: here, I'm use DP basis
    std::vector<MatrixSpinRep_t> g5(convertTwoQuarkSpinDP(15));

    // Load it
    // Form a temporary that holds the adjoint of the prop, but only for a time-slice.
    // Otherwise, we need several temporaries that blows up the storage
    multi2d<LatticeColorVector> tmp1_rep;
    transpose(tmp1_rep, prop_time);

    // Use gamma_5 hermiticity. 
    // Need to form  Gamma(15) * adj(prop) * Gamma(15) = Gamma(15)*tmp1*Gamma(15);
    multi2d<LatticeColorVector> tmp2_rep;
    multiplyRep(tmp2_rep, g5, tmp1_rep);
    multiplyRep(prop_time, tmp2_rep, g5);
  }


  //----------------------------------------------------------------------------
  // Dist(t2) = SpinMatrix*Prop(t2)
  void multiplyRep(multi2d<LatticeColorVector>& dist_rep, 
		   const std::vector<MatrixSpinRep_t>& spin, const multi2d<LatticeColorVector>& prop_rep)
  {
    // Sanity check
    if (prop_rep.nrows() != Ns || prop_rep.ncols() != Ns)
    {
      QDPIO::cerr << __func__ << ": invalid size of prop_rep\n";
      QDP_abort(1);
    }

    // Set the size and initialize
    dist_rep.resize(Ns,Ns);
    dist_rep = zero;

//    for(int sr = 0; sr < Ns; ++sr)
//    {
//      for(int sl = 0; sl < Ns; ++sl)
//      {
//	dist_rep(sl,sr) = zero;
//      }
//    }

    // sparse version of the row and intermediate
    for(int aks = 0; aks < spin.size(); ++aks)
    {
      int left  = spin[aks].left;
      int right = spin[aks].right;
      ComplexD spin_weight(spin[aks].op);
	    
      // Target will be initialize upon reference
      // spin column
      for(int sr = 0; sr < Ns; ++sr)
      {
	dist_rep(left,sr) += spin_weight * prop_rep(right,sr);
      } // sr
    } // aks
  } // void


  //----------------------------------------------------------------------------
  // Dist(t2) = Prop(t2)*SpinMatrix
  void multiplyRep(multi2d<LatticeColorVector>& dist_rep, 
		   const multi2d<LatticeColorVector>& prop_rep, const std::vector<MatrixSpinRep_t>& spin)
  {
    // Sanity check
    if (prop_rep.nrows() != Ns || prop_rep.ncols() != Ns)
    {
      QDPIO::cerr << __func__ << ": invalid size of prop_rep\n";
      QDP_abort(1);
    }

    // Set the size
    dist_rep.resize(Ns,Ns);
    dist_rep = zero;

//    for(int sr = 0; sr < Ns; ++sr)
//    {
//      for(int sl = 0; sl < Ns; ++sl)
//      {
//	dist_rep(sl,sr) = zero;
//      }
//    }

    // sparse version of the column and intermediate
    for(int bks = 0; bks < spin.size(); ++bks)
    {
      int left  = spin[bks].left;
      int right = spin[bks].right;
      ComplexD spin_weight(spin[bks].op);
	    
      // spin row
      for(int sl = 0; sl < Ns; ++sl)
      {
	dist_rep(sl,right) += prop_rep(sl,left) * spin_weight;
      }
    } // bks
  } // void

} // namespace Chroma
