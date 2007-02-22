// $Id: gramschm_array.cc,v 3.1 2007-02-22 21:11:48 bjoo Exp $
/*! \file
 *  \brief Gramm-Schmidt orthogonolization
 */

#include "meas/eig/gramschm_array.h"

namespace Chroma {

//! Gramm-Schmidt orthogonolization
/*!
 * \ingroup eig
 *
 * Templated version
 *
 * Arguments:
 *  \param psi         5D Pseudofermion field     	       (Modify)
 *  \param vec         5D subspace wrt orthog     	       (Read)
 *  \param Nvec        Number of vectors               (Read)
 *  \param Npsi        Number of source vectors        (Read) 
 *  \param sub         Subset to use                   (Read) 
 */

// The multi2d is an annoyance.
// I anticipate that it is called as 
// multi2d T(Nvec, N5)  -- Nvec is "size2" and N5 is "size1"
// T[i] is a vector of size N5
template< typename T>
void GramSchmArray_T(multi2d<T>& psi, const int Npsi,
		     const multi2d<T>& vec, const int Nvec,
		     const Subset& sub) 
{

  START_CODE();

  // if I was paranoid I'd assert that Npsi and Nvec are 
  // reasonable
#if 0
  // Make sure the No of vecs is OK
  if (  Npsi > psi.size2() || Nvec > vec.size2() ) { 
    QDP_error_exit("Npsi and Nvec are out of range in GramSchm. Npsi=%d, Nvec=%d, psi.size2()= %d, vec.size2() = %d", Npsi,Nvec, psi.size2(), vec.size2());
  }

  if( psi.size1() != vec.size1() ) {
    QDP_error_exit("psi and vec have different lengths in 5th Dim. psi.size1() = %d vec.size1() = %d\n", psi.size1(), vec.size1() );
  }
#endif
  int N5 = vec.size1();

  for(int s = 0; s < Npsi; ++s)  
  {
    for(int i = 0; i < Nvec; ++i)   
    {
      // Accumulate the 5D inner product into xp
      Complex xp = innerProduct(vec[i][0], psi[s][0], sub);
      for(int n = 1; n < N5; n++) {
	xp += innerProduct(vec[i][n], psi[s][n], sub);
      }
      
      // Now orthogonalise all N5 components
      for(int n=0; n < N5; n++) { 
	psi[s][n][sub] -= vec[i][n] * xp;
      }
    }
  }
   
  END_CODE();
}

//! Gram Schmidt rothogonalisation
/*!
 * \ingroup eig
 * 
 * Orthogonalise single vector psi against 
 * the first Nvec vectors of vec
 *
 * Templated version
 * Arguments:
 *  \param psi         5D Pseudofermion field     	       (Modify)
 *  \param vec         5D subspace wrt orthog     	       (Read)
 *  \param Nvec        no of vectors to orthog against (Read)
 *  \param sub         Subset to use                   (Read) 
 */

// The multi2d is an annoyance.
// I anticipate that it is called as 
// multi2d T(Nvec, N5)  -- Nvec is "size2" and N5 is "size1"
// T[i] is a vector of size N5
template <typename T>
void GramSchmArray_T(multi1d<T>& psi, 
		     const multi2d<T>& vec, 
		     const int Nvec,
		     const Subset& sub) 
{

  START_CODE();
#if 0 
  if ( Nvec > vec.size2() ) { 
    QDP_error_exit("Nvec out of range in GramSchm. Nvec=%d vec.size2()=%d\n",
		   Nvec, vec.size2());
  }

  if( psi.size() != vec.size1() ) { 
    QDP_error_exit("psi and vec have different lengths in the 5th dimension. psi.size()=%d vec.size1()=%d", psi.size(), vec.size1());
  }
#endif
  int N5 = psi.size();

  for(int i = 0; i < Nvec; ++i)   
  {
    Complex xp = innerProduct(vec[i][0], psi[0], sub);
    for(int n = 0; n < N5; n++) { 
      xp += innerProduct(vec[i][n], psi[n], sub);
    }

    for(int n=0; n < N5; n++)
      psi[n][sub] -= vec[i][n] * xp;
  }

  END_CODE();
}

//! Gram Schmidt rothogonalisation
/*!
 * \ingroup eig
 * 
 * Orthogonalise single vector psi against 
 * vec
 *
 * Templated version
 * Arguments:
 *  \param psi         5D Pseudofermion field     	       (Modify)
 *  \param vec         5D vector to orthog. against     	       (Read)
 *  \param sub         Subset to use                   (Read) 
 */
template <typename T>
void GramSchmArray_T(multi1d<T>& psi, 
		     const multi1d<T>& vec,
		     const Subset& sub) 
{

  START_CODE();

#if 0
  if( psi.size() != vec.size() ) { 
    QDP_error_exit("psi and vec have different lengths in the 5th dimension: psi.size() = %d, vec.size()=%d\n", psi.size(), vec.size());
  }
#endif
   
  int N5 = psi.size();

  Complex xp = innerProduct(vec[0], psi[0], sub);
  for(int i = 1; i < N5; i++) { 
    xp += innerProduct(vec[i], psi[i], sub);
  }

  for(int i=0; i < N5; i++) { 
    psi[i][sub] -= vec[i] * xp;
  }

  END_CODE();
}



//! Gramm-Schmidt orthogonolization
/*!
 * \ingroup eig
 *
 * Arguments:
 *  \param psi         5D Pseudofermion field     	       (Modify)
 *  \param vec         5D subspace wrt orthog     	       (Read)
 *  \param Nvec        Number of vectors               (Read)
 *  \param Npsi        Number of source vectors        (Read) 
 */

void GramSchmArray(multi2d<LatticeFermion>& psi, const int Npsi, 
		   const multi2d<LatticeFermion>& vec, const int Nvec,
		   const Subset& sub) 
{
  GramSchmArray_T(psi, Npsi, vec, Nvec, sub);
}

//! Gram Schmidt rothogonalisation
/*!
 * \ingroup eig
 * 
 * Orthogonalise single vector psi against 
 * the first Nvec vectors of vec
 *
 * Arguments:
 *  \param psi         5D Pseudofermion field     	       (Modify)
 *  \param vec         5D subspace wrt orthog     	       (Read)
 *  \param Nvec        no of vectors to orthog against (Read)
 *  \param sub         Subset to use                   (Read) 
 */
void GramSchmArray(multi1d<LatticeFermion>& psi, 
		   const multi2d<LatticeFermion>& vec, 
		   const int Nvec,
		   const Subset& sub) 
{
  START_CODE();

  GramSchmArray_T(psi, vec, Nvec, sub);

  END_CODE();
}

//! Gram Schmidt rothogonalisation
/*!
 * \ingroup eig
 * 
 * Convenience function: Orthogonalise all vectors of psi against 
 * the first Nvec vectors of vec
 *
 * Arguments:
 *  \param psi         5D Pseudofermion field     	       (Modify)
 *  \param vec         5D subspace wrt orthog     	       (Read)
 *  \param Nvec        Number of vectors               (Read)
 *  \param sub         Subset to use                   (Read) 
 */
void GramSchmArray(multi2d<LatticeFermion>& psi, 
		   const multi2d<LatticeFermion>& vec, const int Nvec,
		   const Subset& sub) 
{
  GramSchmArray_T(psi, psi.size2(), vec, Nvec, sub);
}



//! Gram Schmidt rothogonalisation
/*!
 * \ingroup eig
 * 
 * Convenience function: Orthogonalise all vectors of psi against 
 * the all the vectors of vec
 *
 * Arguments:
 *  \param psi         5D Pseudofermion fields    	       (Modify)
 *  \param vec         5D subspace wrt orthog     	       (Read)
 *  \param sub         Subset to use                   (Read) 
 */
void GramSchmArray(multi2d<LatticeFermion>& psi, 
		   const multi2d<LatticeFermion>& vec,
		   const Subset& sub) 
{
  GramSchmArray_T(psi, psi.size2(), vec, vec.size2(), sub);
}

//! Gram Schmidt rothogonalisation
/*!
 * \ingroup eig
 * 
 * Convenience function: Orthogonalise single vector psi against 
 * all the vectors of vec
 *
 * Arguments:
 *  \param psi         5D Pseudofermion field     	       (Modify)
 *  \param vec         5D subspace wrt orthog     	       (Read)
 *  \param Nvec        no of vectors to orthog against (Read)
 *  \param sub         Subset to use                   (Read) 
 */
void GramSchmArray(multi1d<LatticeFermion>& psi, 
		   const multi2d<LatticeFermion>& vec,
		   const Subset& sub) 
{
  GramSchmArray_T(psi, vec, vec.size2(), sub);
}


//! Gram Schmidt rothogonalisation
/*!
 * \ingroup eig
 * 
 * Orthogonalise single vector psi against 
 * a single vector vec
 *
 * Templated version
 * Arguments:
 *  \param psi         5D Pseudofermion field     	       (Modify)
 *  \param vec         5D subspace wrt orthog     	       (Read)
 *  \param Nvec        no of vectors to orthog against (Read)
 *  \param sub         Subset to use                   (Read) 
 */
void GramSchmArray(multi1d<LatticeFermion>& psi, 
		   const multi1d<LatticeFermion>& vec,
		   const Subset& sub) 
{
  GramSchmArray_T(psi,vec,sub);
}

}  // end namespace Chroma
