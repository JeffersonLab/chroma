#ifndef __gramschm_w__
#define __gramschm_w__

#include "chromabase.h"
/*! \file
 *  \brief Gramm-Schmidt orthogonolization
 */

//! Gram Schmidt rothogonalisation
/*!
 * \ingroup eig
 * Orthogonalise the first Npsi vectors from psi, against the first
 * Nvec vectors from vec
 *
 * Arguments:
 *  \param psi         Pseudofermion field     	       (Modify)
 *  \param vec         subspace wrt orthog     	       (Read)
 *  \param Nvec        Number of vectors               (Read)
 *  \param Npsi        Number of source vectors        (Read) 
 */
void GramSchm(multi1d<LatticeFermion>& psi, 
	      const int Npsi,
	      const multi1d<LatticeFermion>& vec, 
	      const int Nvec);


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
	      const multi1d<LatticeFermion>& vec, 
	      const int Nvec);


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
	      const multi1d<LatticeFermion>& vec);

//! Gram Schmidt rothogonalisation
/*!
 * \ingroup eig
 * 
 * Convenience function: Orthogonalise single vector psi against 
 * the first Nvec vectors of vec
 *
 * Arguments:
 *  \param psi         Pseudofermion field     	       (Modify)
 *  \param vec         subspace wrt orthog     	       (Read)
 *  \param Nvec        no of vectors to orthog against (Read)
 */
void GramSchm(LatticeFermion& psi, 
	      const multi1d<LatticeFermion>& vec, 
	      const int Nvec);


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
	      const multi1d<LatticeFermion>& vec);

void GramSchm(LatticeFermion& psi,
	      const LatticeFermion& vec);
#endif
