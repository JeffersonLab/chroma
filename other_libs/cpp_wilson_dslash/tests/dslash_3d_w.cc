#include "dslash_3d_w.h"

#include "dslashm_w.h"

using namespace QDP;


template<typename T, typename U>
void dslash_3d_T(T& chi, 
	       const multi1d<U>& u, 
	       const T& psi,
	       int isign, int cb3d)
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
    int otherCB = (cb3d == 0 ? 1 : 0);

    switch (isign)
    {
    case 1: 
      {


	LatticeHalfFermion tmp;
	LatticeHalfFermion tmp2;

	tmp[rb3[otherCB]]  = spinProjectDir0Minus(psi);
	tmp2[rb3[cb3d]] = shift(tmp, FORWARD, 0);
	chi[rb3[cb3d]] = spinReconstructDir0Minus(u[0]*tmp2);
	
	
	tmp[rb3[otherCB]]  = spinProjectDir1Minus(psi);
	tmp2[rb3[cb3d]] = shift(tmp, FORWARD, 1);
	chi[rb3[cb3d]] += spinReconstructDir1Minus(u[1]*tmp2);

	tmp[rb3[otherCB]]  = spinProjectDir2Minus(psi);
	tmp2[rb3[cb3d]] = shift(tmp, FORWARD, 2);
	chi[rb3[cb3d]] += spinReconstructDir2Minus(u[2]*tmp2);

	
	tmp[rb3[otherCB]]  = adj(u[0])*spinProjectDir0Plus(psi);
	tmp2[rb3[cb3d]] = shift(tmp, BACKWARD, 0);
	chi[rb3[cb3d]] += spinReconstructDir0Plus(tmp2);
	
	tmp[rb3[otherCB]]  = adj(u[1])*spinProjectDir1Plus(psi);
	tmp2[rb3[cb3d]] = shift(tmp, BACKWARD, 1);
	chi[rb3[cb3d]] += spinReconstructDir1Plus(tmp2);
	
	tmp[rb3[otherCB]]  = adj(u[2])*spinProjectDir2Plus(psi);
	tmp2[rb3[cb3d]] = shift(tmp, BACKWARD, 2);
	chi[rb3[cb3d]] += spinReconstructDir2Plus(tmp2);

      }


      break;

    case (-1):
      {


	LatticeHalfFermion tmp;
	LatticeHalfFermion tmp2;

	tmp[rb3[otherCB]]  = spinProjectDir0Plus(psi);
	tmp2[rb3[cb3d]] = shift(tmp, FORWARD, 0);
	chi[rb3[cb3d]] = spinReconstructDir0Plus(u[0]*tmp2);
	
	tmp[rb3[otherCB]]  = spinProjectDir1Plus(psi);
	tmp2[rb3[cb3d]] = shift(tmp, FORWARD, 1);
	chi[rb3[cb3d]] += spinReconstructDir1Plus(u[1]*tmp2);

	tmp[rb3[otherCB]] = spinProjectDir2Plus(psi);
	tmp2[rb3[cb3d]] = shift(tmp, FORWARD, 2);
	chi[rb3[cb3d]] += spinReconstructDir2Plus(u[2]*tmp2);
	
	
	
	tmp[rb3[otherCB]]  = adj(u[0])*spinProjectDir0Minus(psi);
	tmp2[rb3[cb3d]] = shift(tmp, BACKWARD, 0);
	chi[rb3[cb3d]] += spinReconstructDir0Minus(tmp2);
	
	tmp[rb3[otherCB]]  = adj(u[1])*spinProjectDir1Minus(psi);
	tmp2[rb3[cb3d]] = shift(tmp, BACKWARD, 1);
	chi[rb3[cb3d]] += spinReconstructDir1Minus(tmp2);
	
	tmp[rb3[otherCB]]  = adj(u[2])*spinProjectDir2Minus(psi);
	tmp2[rb3[cb3d]] = shift(tmp, BACKWARD, 2);
	chi[rb3[cb3d]] += spinReconstructDir2Minus(tmp2);

      }
      
      break;
    }
}

void dslash_3d(LatticeFermionF& chi, 
	       const multi1d<LatticeColorMatrixF>& u,
	       const LatticeFermionF& psi,
	       int isign, int cb3d)
{
  dslash_3d_T(chi, u, psi, isign, cb3d);
}

void dslash_3d(LatticeFermionD& chi, 
	       const multi1d<LatticeColorMatrixD>& u,
	       const LatticeFermionD& psi,
	       int isign, int cb3d)
{
 dslash_3d_T(chi, u, psi, isign, cb3d);
}
