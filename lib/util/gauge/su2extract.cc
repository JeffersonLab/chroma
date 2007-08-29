// $Id: su2extract.cc,v 3.2 2007-08-29 13:31:47 edwards Exp $
/*! \file
 *  \brief  Extract an unnormalized SU(2) matrix from a GL(3,C) matrix
 */

#include "chromabase.h"
#include "util/gauge/su2extract.h"

namespace Chroma  
{

  //! Su2_extract: r_0,r_1,r_2,r_3 <- source(su2_index)  [SU(N) field]  under a subset
  /*! 
   * \ingroup gauge
   *
   * Extract components r_k proportional to SU(2) submatrix su2_index
   * from the "SU(Nc)" matrix V. The SU(2) matrix is parametrized in the
   * sigma matrix basis.
   *
   * There are Nc*(Nc-1)/2 unique SU(2) submatrices in an SU(Nc) matrix.
   * The user does not need to know exactly which one is which, just that
   * they are unique.
   *
   * Arguments:
   *
   * \param r           su(2) matrix represented in the O(4) rep. - an array of LatticeReal 
   * \param source      extract the su2 matrix from this su(n) gauge field
   * \param su2_index   int lying in [0, Nc*(Nc-1)/2)
   * \param s           subset for operations       (Read)
   */

  template<typename S>
  inline
  void
  su2Extract_t(multi1d<LatticeReal>& r,
	       const LatticeColorMatrix& source, 
	       int su2_index, 
	       const S& s)
  {
    START_CODE();

    if (r.size() != 4)
    {
      QDPIO::cerr << "su2Extract: return result invalid size" << endl;
      QDP_abort(1);
    }

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

    /* Compute the b(k) of A_SU(2) = b0 + i sum_k bk sigma_k */ 
    r[0][s] = real(peekColor(source,i1,i1)) + real(peekColor(source,i2,i2));
    r[1][s] = imag(peekColor(source,i1,i2)) + imag(peekColor(source,i2,i1));
    r[2][s] = real(peekColor(source,i1,i2)) - real(peekColor(source,i2,i1));
    r[3][s] = imag(peekColor(source,i1,i1)) - imag(peekColor(source,i2,i2));

    END_CODE();
  }


  void
  su2Extract(multi1d<LatticeReal>& r,
	     const LatticeColorMatrix& source, 
	     int su2_index, 
	     const Subset& s)
  {
    su2Extract_t(r, source, su2_index, s);
  }


}  // end namespace Chroma
