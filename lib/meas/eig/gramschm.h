// $Id: gramschm.h,v 1.3 2005-02-10 22:22:42 edwards Exp $
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
	      const OrderedSubset& sub);


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
	      const OrderedSubset& sub);


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
	      const OrderedSubset& sub);

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
	      const OrderedSubset& sub);


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
	      const OrderedSubset& sub);

void GramSchm(LatticeFermion& psi,
	      const LatticeFermion& vec,
	      const OrderedSubset& sub);

}  // end namespace Chroma

#endif
