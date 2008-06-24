// $Id: sunfill.cc,v 3.4 2008-06-24 20:57:24 edwards Exp $
/*! \file
 *  \brief  Fill an SU(Nc) matrix with an SU(2) submatrix
 */

#include "chromabase.h"
#include "util/gauge/sunfill.h"

namespace Chroma 
{

  //! Fill a dest(su2_index) <- r_0,r_1,r_2,r_3  under a subset
  /*!
   * \ingroup gauge
   *
   * Fill an SU(Nc) matrix V with the SU(2) submatrix su2_index
   * paramtrized by b_k in the sigma matrix basis.
   *
   * Fill in B from B_SU(2) = b0 + i sum_k bk sigma_k
   *
   * There are Nc*(Nc-1)/2 unique SU(2) submatrices in an SU(Nc) matrix.
   * The user does not need to know exactly which one is which, just that
   * they are unique.
   *
   * Arguments:
   *
   * \param dest       su(n) matrix
   * \param r          su2 matrix represented in the O(4) rep. - an array of LatticeReal 
   * \param su2_index  int lying in [0, Nc*(Nc-1)/2)
   * \param s          subset for operations       (Read)
   */

  void
  sunFill(LatticeColorMatrix& dest,
	  const multi1d<LatticeReal>& r,
	  int su2_index,
	  const Subset& s)
  {
    START_CODE();

    /* Determine the SU(N) indices corresponding to the SU(2) indices */
    /* of the SU(2) subgroup $3 */
    int i1, i2;
    int found = 0;
    int del_i = 0;
    int index = -1;

    while ( del_i < (Nc-1) && found == 0 )
    {
      del_i++;
      for ( i1 = 0; i1 < (Nc-del_i); i1++ )
      {
	index++;
	if ( index == su2_index )
	{
	  found = 1;
	  break;
	}
      }
    }
    i2 = i1 + del_i;

    if ( found == 0 )
    {
      QDPIO::cerr << __func__ << ": trouble with SU2 subgroup index" << endl;
      QDP_abort(1);
    }

    /* 
     * Insert the b(k) of A_SU(2) = b0 + i sum_k bk sigma_k 
     * back into the SU(N) matrix
     */ 
    dest[s] = 1.0;

    pokeColor(dest[s], cmplx( r[0], r[3]), i1, i1);
    pokeColor(dest[s], cmplx( r[2], r[1]), i1, i2);
    pokeColor(dest[s], cmplx(-r[2], r[1]), i2, i1);
    pokeColor(dest[s], cmplx( r[0],-r[3]), i2, i2);

    END_CODE();
  }

}  // end namespace Chroma

