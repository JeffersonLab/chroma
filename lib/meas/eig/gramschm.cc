// $Id: gramschm.cc,v 1.6 2005-01-14 18:42:35 edwards Exp $
/*! \file
 *  \brief Gramm-Schmidt orthogonolization
 */


#include "meas/eig/gramschm.h"

namespace Chroma {


//! Gramm-Schmidt orthogonolization
/*!
 * \ingroup eig
 *
 * Templated version
 *
 * Arguments:
 *  \param psi         Pseudofermion field     	       (Modify)
 *  \param vec         subspace wrt orthog     	       (Read)
 *  \param Nvec        Number of vectors               (Read)
 *  \param Npsi        Number of source vectors        (Read) 
 */

template< typename T>
void GramSchm_T(multi1d<T>& psi, const int Npsi,
	      const multi1d<T>& vec, const int Nvec) {

  START_CODE();

  // if I was paranoid I'd assert that Npsi and Nvec are 
  // reasonable
#if 0
  if (  Npsi > psi.size() || Nvec > vec.size() ) { 
    QDP_error_exit("Npsi and Nvec are out of range in GramSchm");
  }
#endif

  for(int s = 0; s < Npsi; ++s)  {
    for(int i = 0; i < Nvec; ++i)   {
      Complex xp = innerProduct(vec[i], psi[s]);
      psi[s] -= vec[i] * xp;
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
 *  \param psi         Pseudofermion field     	       (Modify)
 *  \param vec         subspace wrt orthog     	       (Read)
 *  \param Nvec        no of vectors to orthog against (Read)
 */
template <typename T>
void GramSchm_T(T& psi, 
		const multi1d<T>& vec, 
		const int Nvec)
{

  START_CODE();
#if 0 
  if ( Nvec > vec.size() ) { 
    QDP_error_exit("Nvec out of range in GramSchm");
  }
#endif
  for(int i = 0; i < Nvec; ++i)   {
    Complex xp = innerProduct(vec[i], psi);
    psi -= vec[i] * xp;
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
 *  \param psi         Pseudofermion field     	       (Modify)
 *  \param vec         subspace wrt orthog     	       (Read)
 *  \param Nvec        no of vectors to orthog against (Read)
 */
template <typename T>
void GramSchm_T(T& psi, 
		const T& vec)
{

  START_CODE();

  Complex xp = innerProduct(vec, psi);
  psi -= vec * xp;
  
  END_CODE();
}



//! Gramm-Schmidt orthogonolization
/*!
 * \ingroup eig
 *
 * Arguments:
 *  \param psi         Pseudofermion field     	       (Modify)
 *  \param vec         subspace wrt orthog     	       (Read)
 *  \param Nvec        Number of vectors               (Read)
 *  \param Npsi        Number of source vectors        (Read) 
 */

void GramSchm(multi1d<LatticeFermion>& psi, const int Npsi, 
	      const multi1d<LatticeFermion>& vec, const int Nvec)
{
  GramSchm_T(psi, Npsi, vec, Nvec);
}

//! Gram Schmidt rothogonalisation
/*!
 * \ingroup eig
 * 
 * Orthogonalise single vector psi against 
 * the first Nvec vectors of vec
 *
 * Arguments:
 *  \param psi         Pseudofermion field     	       (Modify)
 *  \param vec         subspace wrt orthog     	       (Read)
 *  \param Nvec        no of vectors to orthog against (Read)
 */
void GramSchm(LatticeFermion& psi, 
	      const multi1d<LatticeFermion>& vec, 
	      const int Nvec)
{

  GramSchm_T(psi, vec, Nvec);

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
 *  \param psi         Pseudofermion field     	       (Modify)
 *  \param vec         subspace wrt orthog     	       (Read)
 *  \param Nvec        Number of vectors               (Read)
 */
void GramSchm(multi1d<LatticeFermion>& psi, 
	      const multi1d<LatticeFermion>& vec, const int Nvec)
{
  GramSchm_T(psi, psi.size(), vec, Nvec);
}



//! Gram Schmidt rothogonalisation
/*!
 * \ingroup eig
 * 
 * Convenience function: Orthogonalise all vectors of psi against 
 * the all the vectors of vec
 *
 * Arguments:
 *  \param psi         Pseudofermion field     	       (Modify)
 *  \param vec         subspace wrt orthog     	       (Read)
 */
void GramSchm(multi1d<LatticeFermion>& psi, 
	      const multi1d<LatticeFermion>& vec)
{
  GramSchm_T(psi, psi.size(), vec, vec.size());
}

//! Gram Schmidt rothogonalisation
/*!
 * \ingroup eig
 * 
 * Convenience function: Orthogonalise single vector psi against 
 * all the vectors of vec
 *
 * Arguments:
 *  \param psi         Pseudofermion field     	       (Modify)
 *  \param vec         subspace wrt orthog     	       (Read)
 *  \param Nvec        no of vectors to orthog against (Read)
 */
void GramSchm(LatticeFermion& psi, 
	      const multi1d<LatticeFermion>& vec)
{
  GramSchm_T(psi, vec, vec.size());
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
 *  \param psi         Pseudofermion field     	       (Modify)
 *  \param vec         subspace wrt orthog     	       (Read)
 *  \param Nvec        no of vectors to orthog against (Read)
 */
void GramSchm(LatticeFermion& psi, 
		const LatticeFermion& vec)
{

  START_CODE();

  Complex xp = innerProduct(vec, psi);
  psi -= vec * xp;
  
  END_CODE();
}

}  // end namespace Chroma
