// $Id: invbicgstab_array.cc,v 1.1 2004-05-19 00:21:23 bjoo Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#include "chromabase.h"
#include "actions/ferm/invert/invbicgstab_array.h"


template<typename T>
void InvBiCGStab_a(const LinearOperator< multi1d<T> >& A,
		   const multi1d<T> & chi,
		   multi1d<T>& psi,
		   const Real& RsdCG, 
		   int MaxCG, 
		   int& n_count)
{
  const int N = psi.size();
  const OrderedSubset& s = A.subset();

  Real chi_sq =  Real(norm2(chi,s));
  QDPIO::cerr << "Not yet implemented " << endl;
  QDP_abort(1);
}


// Fix here for now
template<>
void InvBiCGStab(const LinearOperator< multi1d<LatticeFermion> >& A,
		 const multi1d<LatticeFermion>& chi,
		 multi1d<LatticeFermion>& psi,
		 const Real& RsdCG, 
		 int MaxCG, 
		 int& n_count)
{
  InvBiCGStab_a(A, chi, psi, RsdCG, MaxCG, n_count);
}

