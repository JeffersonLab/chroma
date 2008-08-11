// -*- C++ -*-
// $Id: eoprec_staggered_qprop.h,v 3.5 2008-08-11 10:45:35 mcneile Exp $
/*! \file
 *  \brief Propagator solver for an even-odd non-preconditioned fermion operator
 *
 *  Solve for the propagator of an even-odd non-preconditioned fermion operator
 */

#ifndef PREC_STAGGERED_QPROP_H
#define PREC_STAGGERED_QPROP_H

#include "stagtype_fermact_s.h"
#include "actions/ferm/invert/invcg1.h"
#include "actions/ferm/invert/syssolver_cg_params.h"


namespace Chroma
{
  //! Propagator of a generic even-odd fermion linear operator
  /*! \ingroup qprop
   *
   * This routine is actually generic to all even-odd fermions
   */
  template<typename T, typename P, typename Q>
  class EvenOddFermActQprop : public SystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param A_        M^dag*M operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    EvenOddFermActQprop(Handle< EvenOddLinearOperator<T,P,Q> > M_,
			Handle< LinearOperator<T> > A_,
			const Real& Mass_,
			const SysSolverCGParams& invParam_) : 
      M(M_), A(A_), Mass(Mass_), invParam(invParam_) 
      {}

    EvenOddFermActQprop(const EvenOddStaggeredTypeFermAct<T,P,Q>& S_,
			Handle< FermState<T,P,Q> > state,
			const SysSolverCGParams& invParam_) :
      M(S_.linOp(state)), A(S_.lMdagM(state)),
      Mass(S_.getQuarkMass()), invParam(invParam_) {}

    //! Destructor is automatic
    ~EvenOddFermActQprop() {}

    //! Return the subset on which the operator acts
    const Subset& subset() const {return all;}

    //! Solver the linear system
    /*!
     * \param psi      quark propagator ( Modify )
     * \param chi      source ( Read )
     * \return number of CG iterations
     */
    SystemSolverResults_t operator() (T& psi, const T& chi) const
    {
      START_CODE();

      LatticeStaggeredFermion tmp, tmp1, tmp2;
      tmp = tmp1 = tmp2 = zero;
      Real invm;

      //  switch(invType)
      //{
      //case CG_INVERTER: 

      // Make preconditioned source:  tmp_1_e = M_ee chi_e + M_eo^{dag} chi_o

      M->evenEvenLinOp(tmp, chi, PLUS);
      M->evenOddLinOp(tmp2, chi, MINUS);
      tmp[rb[0]] += tmp2;
    

      /* psi = (M^dag * M)^(-1) chi  = A^{-1} chi*/
      SystemSolverResults_t res = InvCG1(*A, tmp, psi, 
					 invParam.RsdCG, 
					 invParam.MaxCG,
					 invParam.MinCG);
      
      // psi[rb[0]] is returned, so reconstruct psi[rb[1]]
      invm = Real(1)/(2*Mass);
    
      // tmp_1_o = D_oe psi_e 
      M->oddEvenLinOp(tmp1, psi, PLUS);

      // tmp_1_o = (1/2m) D_oe psi_e
      tmp1[rb[1]] *= invm;

      // tmp_2_o = (1/2m) chi_o
      tmp2[rb[1]]  = invm * chi;

      // psi_o = (1/2m) chi_o - (1/2m) D_oe psi_e 
      psi[rb[1]] = tmp2 - tmp1;
      //    break;  
  
      //  default:
      // QDP_error_exit("Unknown inverter type", invType);
      //}

      if ( res.n_count == invParam.MaxCG )
	QDP_error_exit("no convergence in the inverter", res.n_count);

      // Compute residual
      {
	T  r;
	(*M)(r, psi, PLUS);
	r -= chi;
	res.resid = sqrt(norm2(r));
	 QDPIO::cout << "eoprec_staggered_qprop:  true residual:  " << res.resid << endl;
      }

      END_CODE();

      return res;
    }


  private:
    // Hide default constructor
    EvenOddFermActQprop() {}

    Handle< EvenOddLinearOperator<T,P,Q> > M;
    Handle< LinearOperator<T> > A;
    Real Mass;
    SysSolverCGParams invParam; 
  };

}; // End namespace

#endif 

