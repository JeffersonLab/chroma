// $Id: sn_jacob_array.h,v 3.1 2007-02-22 21:11:49 bjoo Exp $

#ifndef __sn_jacob_array_h__
#define __sn_jacob_array_h__

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
 *  \param sub         Subset to use                   (Read) 
 *
 * \return 	Number of Jacobi iters		(Write) 
 */
int SN_Jacob_Array(multi2d<LatticeFermion>& psi, 
		   const int N_eig, 
		   multi1d<Real>& lambda, 
		   multi1d<Complex>& off_diag, 
		   Real tolerance, 
		   int N_max,
		   const Subset& sub);

}  // end namespace Chroma

#endif
