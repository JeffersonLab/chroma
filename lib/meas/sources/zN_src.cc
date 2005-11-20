// $Id: zN_src.cc,v 2.1 2005-11-20 18:28:38 edwards Exp $
/*! \file
 *  \brief Variety of Z(N) noise sources
 */

#include "chromabase.h"
#include "meas/sources/zN_src.h"

namespace Chroma 
{

  //! Volume source of Z(N) noise
  /*!
   * \ingroup sources
   *
   * This routine is specific to Wilson fermions! 
   *
   * \param a      Source fermion
   * \param N      The N in Z(N)
   *
   *
   *  This type of source is required to compute disconnected
   *  diagrams. The source is complex Z(N) noise, hence there
   *  is an additional normalization factor of 1/sqrt(N) or somesuch.
   *
   *
   */

  void zN_src(LatticeFermion & a, int N)
  {
    a = zero; 
    LatticeReal rnd1, theta;
    Real twopiN = Chroma::twopi / N;   // twopi defined in chroma/lib/chromabase.h
    LatticeComplex c;
    LatticeColorVector colorvec = zero;

    for(int spin_index= 0; spin_index < Ns; ++spin_index)
      for(int color_index= 0; color_index < Nc; ++color_index)
      {
	random(rnd1); 
	theta = twopiN * floor(N*rnd1);
	c = cmplx(cos(theta),sin(theta));

	colorvec = peekSpin(a,spin_index);

	pokeSpin(a,pokeColor(colorvec,c,color_index),spin_index);
      }

  }


}  // end namespace Chroma

