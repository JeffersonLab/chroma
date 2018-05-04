// -*- C++ -*-
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#ifndef __inv_rel_sumr__
#define __inv_rel_sumr__

#include "chromabase.h"
#include "linearop.h"

namespace Chroma {

// Solve a shifted unitary system
//
// A x = b
//
// Where X is of the form: A = zeta I + rho U
//
// rho > 0 and zeta are complex, and U is unitary
//
// We solve with the method described in:
//
// "A Fast Minimal Residual Algorithm for Shifted Unitary Matrices"
// by Carl F. Jagels, and Lothar Reichel
// Numerical Linear Algebra with Applications, Vol 1(6), 555-570(1994)
//
// This paper is referenced by and applied to the Overlap Dirac Operator
// by G. Arnold, N. Cundy, J. van den Eshof, A Frommer, S. Krieg, T. Lippert,
// K. Schaefer "Numerical Methods for the QCD Overlap Operator: II. 
// Optimal Krylov Subspace Methods" -- hep-lat/0311025
// which is where the name SUMR was coined.
//
template<typename T>
void InvRelSUMR(const LinearOperator<T>& U,
		const T& b,
		T& x,
		const Complex& zeta,
		const Real& rho,
		const Real& epsilon, 
		int MaxSUMR, 
		int& n_count);

}  // end namespace Chroma

#endif
