// $Id: invbicgstab.cc,v 1.1 2004-05-19 00:21:23 bjoo Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#include "chromabase.h"
#include "actions/ferm/invert/invbicgstab.h"


template<typename T>
void InvBiCGStab_a(const LinearOperator<T>& A,
		   const T& chi,
		   T& psi,
		   const Real& RsdCG, 
		   int MaxCG, 
		   int& n_count)
{
  const OrderedSubset& s = A.subset();

  Real rsd_sq = (RsdCG * RsdCG) * Real(norm2(chi,s));

  QDPIO::cerr << "Not yet implemented " << endl;
  QDP_abort(1);
}


// Fix here for now
template<>
void InvBiCGStab(const LinearOperator<LatticeFermion>& A,
		 const LatticeFermion& chi,
		 LatticeFermion& psi,
		 const Real& RsdCG, 
		 int MaxCG, 
		 int& n_count)
{
  InvBiCGStab_a(A, chi, psi, RsdCG, MaxCG, n_count);
}

