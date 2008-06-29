// $Id: gramschm.h,v 3.2 2008-06-29 20:33:23 edwards Exp $
/*! \file
 *  \brief Gramm-Schmidt orthogonolization
 */

#ifndef __gramschm_w__
#define __gramschm_w__

#include "chromabase.h"

namespace Chroma {

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
 *  \param sub         Subset to use                   (Read) 
 */
void GramSchm(multi1d<LatticeFermion>& psi, 
	      const int Npsi,
	      const multi1d<LatticeFermion>& vec, 
	      const int Nvec,
	      const Subset& sub);


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
	      const multi1d<LatticeFermion>& vec, 
	      const int Nvec,
	      const Subset& sub);


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
	      const Subset& sub);

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
 *  \param sub         Subset to use                   (Read) 
 */
void GramSchm(LatticeFermion& psi, 
	      const multi1d<LatticeFermion>& vec, 
	      const int Nvec,
	      const Subset& sub);


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
 *  \param sub         Subset to use                   (Read) 
 */
void GramSchm(LatticeColorVector& psi, 
	      const multi1d<LatticeColorVector>& vec, 
	      const int Nvec,
	      const Subset& sub);


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
	      const Subset& sub);

void GramSchm(LatticeFermion& psi,
	      const LatticeFermion& vec,
	      const Subset& sub);

}  // end namespace Chroma

#endif
