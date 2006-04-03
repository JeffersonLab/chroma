// $Id: z2_src.cc,v 3.0 2006-04-03 04:59:06 edwards Exp $
/*! \file
 *  \brief Variety of Z2 noise sources
 */

#include "chromabase.h"
#include "meas/sources/z2_src.h"

namespace Chroma {

//! Volume source of complex Z2 noise
/*!
 * \ingroup sources
 *
 * This routine is specific to Wilson fermions! 
 *
 * \param a      Source fermion
 *
 *
 *  This type of source is required to compute disconnected
 *  diagrams. The source is complex Z2 noise, hence there
 *  is an additional normalization factor of 1/sqrt(2).
 *
 *
 */

  template<typename T>
  void z2_src_t(T & a)
  {

    a = zero ; 
    LatticeReal rnd ;
    LatticeReal ar ,ai ;
    LatticeColorVector colorvec = zero;

    for(int spin_index= 0 ; spin_index < Ns ; ++spin_index)
      for(int color_index= 0 ; color_index < Nc ; ++color_index)
      {
	random(rnd) ; 
	ar = where( rnd > 0.5 , LatticeReal(1) , LatticeReal(-1) );
	random(rnd) ; 
	ai = where( rnd > 0.5 , LatticeReal(1) , LatticeReal(-1) );
	LatticeComplex c = cmplx(ar,ai) ;

	colorvec = peekSpin(a,spin_index);

	pokeSpin(a,pokeColor(colorvec,c,color_index),spin_index);
      }

  }


  // template<>
  void z2_src(LatticeFermion& a)
  {
    z2_src_t<LatticeFermion>(a) ; 
  }


  // template<>
  void z2_src(LatticeStaggeredFermion& a)
  {
    z2_src_t<LatticeStaggeredFermion>(a) ; 
  }


  //! Timeslice source of complex Z2 noise
  /*!
   * \ingroup sources
   *
   * This routine is specific to Wilson fermions! 
   *
   * \param a      Source fermion
   * \param slice        time slice
   * \param mu           direction of slice
   *
   *
   *  This type of source is useful for computing hadronic
   *  decay like diagrams. 
   *
   *  QUARK MASS DEPENDENCE OF HADRON MASSES FROM LATTICE QCD.
   * By UKQCD Collaboration (M. Foster et al.)
   * Published in Phys.Rev.D59:074503,1999
   * e-Print Archive: hep-lat/9810021 
   *
   */
  void z2_src(LatticeFermion& a, int slice, int mu)
  {

    // compute a volume source of z2 noise
    LatticeFermion tmp; 
    z2_src(tmp) ;

    a = where(Layout::latticeCoordinate(mu) == slice, tmp, LatticeFermion(zero));
  }

}  // end namespace Chroma

