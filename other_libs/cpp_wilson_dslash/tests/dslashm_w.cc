#include "dslashm_w.h"

using namespace QDP;


template<typename T, typename H, typename U>
void dslash_a(T& chi, 
	    const multi1d<U>& u, 
	    const T& psi,
	    int isign, int cb)
{

    /*     F 
     *   a2  (x)  :=  U  (x) (1 - isign gamma  ) psi(x)
     *     mu          mu                    mu
     */
    /*     B           +
     *   a2  (x)  :=  U  (x-mu) (1 + isign gamma  ) psi(x-mu)
     *     mu          mu                       mu
     */
    // Recontruct the bottom two spinor components from the top two
    /*                        F           B
     *   chi(x) :=  sum_mu  a2  (x)  +  a2  (x)
     *                        mu          mu
     */

    /* Why are these lines split? An array syntax would help, but the problem is deeper.
     * The expression templates require NO variable args (like int's) to a function
     * and all args must be known at compile time. Hence, the function names carry
     * (as functions usually do) the meaning (and implicit args) to a function.
     */

    int otherCB = (cb == 0 ? 1 : 0);

    switch (isign)
    {
    case 1: 
      {


	H tmp;
	H tmp2;

	tmp[rb[otherCB]]  = spinProjectDir0Minus(psi);
	tmp2[rb[cb]] = shift(tmp, FORWARD, 0);
	chi[rb[cb]] = spinReconstructDir0Minus(u[0]*tmp2);
	
	
	tmp[rb[otherCB]]  = spinProjectDir1Minus(psi);
	tmp2[rb[cb]] = shift(tmp, FORWARD, 1);
	chi[rb[cb]] += spinReconstructDir1Minus(u[1]*tmp2);

	tmp[rb[otherCB]]  = spinProjectDir2Minus(psi);
	tmp2[rb[cb]] = shift(tmp, FORWARD, 2);
	chi[rb[cb]] += spinReconstructDir2Minus(u[2]*tmp2);
	
	tmp[rb[otherCB]]  = spinProjectDir3Minus(psi);
	tmp2[rb[cb]] = shift(tmp, FORWARD, 3);
	chi[rb[cb]] += spinReconstructDir3Minus(u[3]*tmp2);
	
	
	tmp[rb[otherCB]]  = adj(u[0])*spinProjectDir0Plus(psi);
	tmp2[rb[cb]] = shift(tmp, BACKWARD, 0);
	chi[rb[cb]] += spinReconstructDir0Plus(tmp2);
	
	tmp[rb[otherCB]]  = adj(u[1])*spinProjectDir1Plus(psi);
	tmp2[rb[cb]] = shift(tmp, BACKWARD, 1);
	chi[rb[cb]] += spinReconstructDir1Plus(tmp2);
	
	tmp[rb[otherCB]]  = adj(u[2])*spinProjectDir2Plus(psi);
	tmp2[rb[cb]] = shift(tmp, BACKWARD, 2);
	chi[rb[cb]] += spinReconstructDir2Plus(tmp2);

	tmp[rb[otherCB]]  = adj(u[3])*spinProjectDir3Plus(psi);
	tmp2[rb[cb]] = shift(tmp, BACKWARD, 3);
	chi[rb[cb]] += spinReconstructDir3Plus(tmp2);	
	
      }


      break;

    case (-1):
      {


	H tmp;
	H tmp2;

	tmp[rb[otherCB]]  = spinProjectDir0Plus(psi);
	tmp2[rb[cb]] = shift(tmp, FORWARD, 0);
	chi[rb[cb]] = spinReconstructDir0Plus(u[0]*tmp2);
	
	tmp[rb[otherCB]]  = spinProjectDir1Plus(psi);
	tmp2[rb[cb]] = shift(tmp, FORWARD, 1);
	chi[rb[cb]] += spinReconstructDir1Plus(u[1]*tmp2);

	tmp[rb[otherCB]] = spinProjectDir2Plus(psi);
	tmp2[rb[cb]] = shift(tmp, FORWARD, 2);
	chi[rb[cb]] += spinReconstructDir2Plus(u[2]*tmp2);
	
	tmp[rb[otherCB]]  = spinProjectDir3Plus(psi);
	tmp2[rb[cb]] = shift(tmp, FORWARD, 3);
	chi[rb[cb]] += spinReconstructDir3Plus(u[3]*tmp2);
	
	
	tmp[rb[otherCB]]  = adj(u[0])*spinProjectDir0Minus(psi);
	tmp2[rb[cb]] = shift(tmp, BACKWARD, 0);
	chi[rb[cb]] += spinReconstructDir0Minus(tmp2);
	
	tmp[rb[otherCB]]  = adj(u[1])*spinProjectDir1Minus(psi);
	tmp2[rb[cb]] = shift(tmp, BACKWARD, 1);
	chi[rb[cb]] += spinReconstructDir1Minus(tmp2);
	
	tmp[rb[otherCB]]  = adj(u[2])*spinProjectDir2Minus(psi);
	tmp2[rb[cb]] = shift(tmp, BACKWARD, 2);
	chi[rb[cb]] += spinReconstructDir2Minus(tmp2);

	tmp[rb[otherCB]]  = adj(u[3])*spinProjectDir3Minus(psi);    
	tmp2 = shift(tmp, BACKWARD, 3);
	chi[rb[cb]] += spinReconstructDir3Minus(tmp2);	
	
      }

      break;
    }

}

void dslash(LatticeFermionF& chi, 
	    const multi1d<LatticeColorMatrixF>& u, 
	    const LatticeFermionF& psi,
	    int isign, int cb)
{
  dslash_a<LatticeFermionF, LatticeHalfFermionF, LatticeColorMatrixF>(chi, u, psi, isign, cb);
}


void dslash(LatticeFermionD& chi, 
	    const multi1d<LatticeColorMatrixD>& u, 
	    const LatticeFermionD& psi,
	    int isign, int cb)
{
  dslash_a<LatticeFermionD, LatticeHalfFermionD, LatticeColorMatrixD>(chi, u, psi, isign, cb);
}


