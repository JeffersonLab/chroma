// $Id: gramschm.cc,v 3.2 2008-06-29 20:33:23 edwards Exp $
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
 *  \param sub         Subset to use                   (Read) 
 */

template< typename T>
void GramSchm_T(multi1d<T>& psi, const int Npsi,
		const multi1d<T>& vec, const int Nvec,
		const Subset& sub)
{
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
      Complex xp = innerProduct(vec[i], psi[s], sub);
      psi[s][sub] -= vec[i] * xp;
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
 *  \param sub         Subset to use                   (Read) 
 */
template <typename T>
void GramSchm_T(T& psi, 
		const multi1d<T>& vec, 
		const int Nvec,
		const Subset& sub) 
{

  START_CODE();
#if 0 
  if ( Nvec > vec.size() ) { 
    QDP_error_exit("Nvec out of range in GramSchm");
  }
#endif
  for(int i = 0; i < Nvec; ++i)   {
    Complex xp = innerProduct(vec[i], psi, sub);
    psi[sub] -= vec[i] * xp;
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
 *  \param sub         Subset to use                   (Read) 
 */
template <typename T>
void GramSchm_T(T& psi, 
		const T& vec,
		const Subset& sub) 
{

  START_CODE();

  Complex xp = innerProduct(vec, psi, sub);
  psi[sub] -= vec * xp;
  
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
 *  \param sub         Subset to use                   (Read) 
 */

void GramSchm(multi1d<LatticeFermion>& psi, const int Npsi, 
	      const multi1d<LatticeFermion>& vec, const int Nvec,
	      const Subset& sub) 
{
  GramSchm_T(psi, Npsi, vec, Nvec, sub);
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
 *  \param sub         Subset to use                   (Read) 
 */
void GramSchm(LatticeFermion& psi, 
	      const multi1d<LatticeFermion>& vec, 
	      const int Nvec,
	      const Subset& sub) 
{
  START_CODE();

  GramSchm_T(psi, vec, Nvec, sub);

  END_CODE();
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
 *  \param sub         Subset to use                   (Read) 
 */
void GramSchm(LatticeColorVector& psi, 
	      const multi1d<LatticeColorVector>& vec, 
	      const int Nvec,
	      const Subset& sub) 
{
  START_CODE();

  GramSchm_T(psi, vec, Nvec, sub);

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
 *  \param sub         Subset to use                   (Read) 
 */
void GramSchm(multi1d<LatticeFermion>& psi, 
	      const multi1d<LatticeFermion>& vec, const int Nvec,
	      const Subset& sub) 
{
  GramSchm_T(psi, psi.size(), vec, Nvec, sub);
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
 *  \param sub         Subset to use                   (Read) 
 */
void GramSchm(multi1d<LatticeFermion>& psi, 
	      const multi1d<LatticeFermion>& vec,
	      const Subset& sub) 
{
  GramSchm_T(psi, psi.size(), vec, vec.size(), sub);
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
 *  \param sub         Subset to use                   (Read) 
 */
void GramSchm(LatticeFermion& psi, 
	      const multi1d<LatticeFermion>& vec,
	      const Subset& sub)
{
  GramSchm_T(psi, vec, vec.size(), sub);
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
 *  \param sub         Subset to use                   (Read) 
 */
void GramSchm(LatticeFermion& psi, 
	      const LatticeFermion& vec,
	      const Subset& sub) 
{

  START_CODE();

  Complex xp = innerProduct(vec, psi, sub);
  psi[sub] -= vec * xp;
  
  END_CODE();
}

}  // end namespace Chroma
