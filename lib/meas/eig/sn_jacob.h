// $Id: sn_jacob.h,v 1.2 2005-01-14 18:42:35 edwards Exp $
#ifndef __sn_jacob_h__
#define __sn_jacob_h__

namespace Chroma {

//! Single-node Jacobi rotation
/*!
 * \ingroup eig
 *
 * This subroutine contains a "single node" Jacobi routine
 * to be used with the Ritz functional eigenvialue/vector finder.
 *
 *
 *  Psi		Eigenvectors			(Modify)
 *  N_eig	Eigenvalue number 		(Read)
 *  lambda	Diagonals / Eigenvalues		(Modify)
 *  off_diag	Upper triang off-diag matrix elems	(Modify)
 *  Toler	Tolerance for off-diag elems	(Read)
 *  N_max	Maximal number of Jacobi iters	(Read)
 *  Ncb		Number of sublattices		(Read)
 *  N_Count	Number of Jacobi iters		(Write) 
 */
template <typename T>
void SN_Jacob_t(multi1d<T>& psi, 
	      const int N_eig, 
	      multi1d<Real>& lambda, 
	      multi1d<Complex>& off_diag, 
	      Real tolerance, 
	      int N_max,
	      int& n_count);


//! Single-node Jacobi rotation
/*!
 * \ingroup eig
 *
 * This subroutine contains a "single node" Jacobi routine
 * to be used with the Ritz functional eigenvialue/vector finder.
 *
 *
 *  Psi		Eigenvectors			(Modify)
 *  N_eig	Eigenvalue number 		(Read)
 *  lambda	Diagonals / Eigenvalues		(Modify)
 *  off_diag	Upper triang off-diag matrix elems	(Modify)
 *  Toler	Tolerance for off-diag elems	(Read)
 *  N_max	Maximal number of Jacobi iters	(Read)
 *  Ncb		Number of sublattices		(Read)
 *  N_Count	Number of Jacobi iters		(Write) 
 */
void SN_Jacob(multi1d<LatticeFermion>& psi, 
	      const int N_eig, 
	      multi1d<Real>& lambda, 
	      multi1d<Complex>& off_diag, 
	      Real tolerance, 
	      int N_max,
	      int& n_count);

}  // end namespace Chroma

#endif
