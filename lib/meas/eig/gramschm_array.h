// $Id: gramschm_array.h,v 3.1 2007-02-22 21:11:48 bjoo Exp $
/*! \file
 *  \brief Gramm-Schmidt orthogonolization
 */
#ifndef __gramschm_array_w__
#define __gramschm_array_w__

#include "chromabase.h"

namespace Chroma {

//! Gram Schmidt rothogonalisation
/*!
 * \ingroup eig
 * Orthogonalise the first Npsi vectors from psi, against the first
 * Nvec vectors from vec
 *
 * Arguments:
 *  \param psi         5D Pseudofermion field  	       (Modify)
 *  \param vec         5D subspace wrt orthog  	       (Read)
 *  \param Nvec        Number of vectors               (Read)
 *  \param Npsi        Number of source vectors        (Read) 
 *  \param s           Subset to use                   (Read) 
 */
void GramSchmArray(multi2d<LatticeFermion>& psi, 
		   const int Npsi,
		   const multi2d<LatticeFermion>& vec, 
		   const int Nvec,
		   const Subset& s);


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
 *  \param s           Subset to use                   (Read) 
 */
void GramSchmArray(multi2d<LatticeFermion>& psi, 
		   const multi2d<LatticeFermion>& vec, 
		   const int Nvec,
		   const Subset& s);


//! Gram Schmidt rothogonalisation
/*!
 * \ingroup eig
 * 
 * Convenience function: Orthogonalise all vectors of psi against 
 * the all the vectors of vec
 *
 * Arguments:
 *  \param psi         5D Pseudofermion field     	       (Modify)
 *  \param vec         5D subspace wrt orthog     	       (Read)
 *  \param s           Subset to use                   (Read) 
 */
void GramSchmArray(multi2d<LatticeFermion>& psi, 
		   const multi2d<LatticeFermion>& vec,
		   const Subset& s);


//! Gram Schmidt rothogonalisation
/*!
 * \ingroup eig
 * 
 * Convenience function: Orthogonalise single vector psi against 
 * the first Nvec vectors of vec
 *
 * Arguments:
 *  \param psi         5D Pseudofermion field     	       (Modify)
 *  \param vec         5D subspace wrt orthog     	       (Read)
 *  \param Nvec        no of vectors to orthog against (Read)
 *  \param s           Subset to use                   (Read) 
 */
void GramSchmArray(multi1d<LatticeFermion>& psi, 
		   const multi2d<LatticeFermion>& vec, 
		   const int Nvec,
		   const Subset& s);


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
 *  \param s           Subset to use                   (Read) 
 */
void GramSchmArray(multi1d<LatticeFermion>& psi, 
		   const multi2d<LatticeFermion>& vec,
		   const Subset& s);

//! Gram Schmidt rothogonalisation
/*!
 * \ingroup eig
 * 
 * Convenience function: Orthogonalise single vector psi against 
 * a single vector vec
 *
 * Arguments:
 *  \param psi         5D Pseudofermion field          (Modify)
 *  \param vec         5D vector  wrt orthog           (Read)
 *  \param s           Subset to use                   (Read) 
 */
void GramSchmArray(multi1d<LatticeFermion>& psi,
		   const multi1d<LatticeFermion>& vec,
		   const Subset& s);

}  // end namespace Chroma

#endif
