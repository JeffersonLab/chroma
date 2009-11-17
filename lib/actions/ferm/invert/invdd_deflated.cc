// $Id: invdd_deflated.cc,v 1.1 2009/11/11 16:39:50 jbulava Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm for a generic Linear Operator
 */

#include "chromabase.h"
#include "actions/ferm/invert/invdd_deflated.h"


namespace Chroma {

	namespace InvDDDeflatedEnv {


		/*****************************************************
		 * A : The linear operator which is to be inverted 
		 * chi : The source on which to invert
		 * psi : the result of the inversion
		 * vecs : The vectors which span the deflation subspace
		 * blk : The definition of the domain decompostion
		 * RsdCG : The maximim allowed residual
		 * MaxCG : The maximum number of iterations
		 * *******************************************************/

  SystemSolverResults_t 
ddDeflInv(const LinearOperator<LatticeFermion>& A,
       const LatticeFermion& chi,
       LatticeFermion& psi,
			 const mutli1d<LatticeFermion>& vecs,
			 const QDP::Set blk,
       const Real& RsdCG, 
       int MaxCG) {

		START_CODE();

		//The struct which contains inversion information
  SystemSolverResults_t  res;

	//First, form the projection of chi onto the deflation subspace
	int nvecs_per_block = vecs.size();
	int nblocks = blk.getNumSubsets();
	multi1d<Complex> chi_little(nblocks * nvecs_per_block);
	int cntr = 0;
	for (int v = 0 ; v < nvecs_per_block ; ++v) {
		LatticeComplex prod = localInnerProduct( vecs[v] , chi);
		multi1d<Complex> temp = sumMulti( prod, blk);
		for (int b = 0 ; b < nblocks ; ++b) {
			chi_little[cntr] = temp[b];
			cntr++;
		}
	}

	//Next, solve the little dirac equation with chi_little as a source
	multi1d<Complex> ainv_chi_little;
	SystemSolverResults_t lres = solve_little( A, chi_little, ainv_chi_little, 
			vecs, blk, RsdCG, MaxCG);

	//Now Form P_L * \chi
	///!!!!Could be sped up???
	LatticeFermion pl_chi = chi;
	cntr = 0;
	for (int v = 0 ; v < nvecs_per_block ; ++v) 
			for (int b = 0 ; b < nblocks ; ++b) {
				LatticeFermion temp = zero;
				//Need to restrict the action of A to a subspace, I hope this works
				A(temp[blk[b]], vecs[v], Plus); 
				pl_chi[blk[b]] -= temp * ainv_chi_little[cntr];
				cntr++;
			}

	//Now send P_L * \chi to a normal solver
	SystemSolverResults_t dres = any_old_solve(A, psi, pl_chi);
	
	//Act with P_R on psi
	LatticeFermion dpsi = zero;
	A(dpsi, psi, Plus);
	multi1d<Complex> dpsi_little(nblocks * nvecs_per_block);
	cntr = 0;
	for (int v = 0 ; v < nvecs_per_block ; ++v) {
		LatticeComplex prod = localInnerProduct( vecs[v] , dpsi);
		multi1d<Complex> temp = sumMulti( prod, blk);
		for (int b = 0 ; b < nblocks ; ++b) { 
			dpsi_little[cntr] = temp[b];	
			cntr++;
		}
	}

	cntr = 0; 
	for (int v = 0 ; v < nvecs_per_block ; ++v) 
		for (int b = 0 ; b < nblocks ; ++b) {
			psi[blk[b]] -= vecs[v] * dpsi_little[cntr];
			cntr++;
		}
		
	//Add the little solution to the orthogonal compliment solution
	cntr = 0;
	for (int v = 0 ; v < nvecs_per_block ; ++v) 
		for (int b = 0 ; b < nblocks ; ++b)
		{
			psi[blk[b]] += ainv_chi_little[cntr] * vecs[v];
			cntr++;
		}

   END_CODE();
  return res;
}

} //end namespace
}  // end namespace Chroma
