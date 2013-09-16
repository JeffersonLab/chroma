// $Id: eoprec_fermact_qprop.cc,v 3.3 2009-03-26 00:07:25 edwards Exp $
/*! \file
 *  \brief Propagator solver for a generic even-odd preconditioned fermion operator
 *
 *  Solve for the propagator of a even-odd non-preconditioned fermion operator
 */

#include "eoprec_wilstype_fermact_w.h"

namespace Chroma 
{ 
  //! Propagator of a generic even-odd preconditioned fermion linear operator
  /*! \ingroup qprop
   *
   * This routine is actually generic to all even-odd preconditioned fermions
   */
  template<typename T, typename P, typename Q>
  class PrecFermActQprop : public SystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param A_         Linear operator ( Read )
     * \param invParam_  inverter parameters ( Read )
     */
    PrecFermActQprop(Handle< EvenOddPrecLinearOperator<T,P,Q> > A_,
		     Handle< LinOpSystemSolver<T> > invA_) : A(A_), invA(invA_) 
      {}

    //! Destructor is automatic
    ~PrecFermActQprop() {}

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

      /* Step (i) */
      /* chi_tmp =  chi_o - D_oe * A_ee^-1 * chi_e */
      T chi_tmp;
      {
	T tmp1, tmp2;

	A->evenEvenInvLinOp(tmp1, chi, PLUS);
	A->oddEvenLinOp(tmp2, tmp1, PLUS);
	chi_tmp[rb[1]] = chi - tmp2;
      }

      // Call inverter
      SystemSolverResults_t res = (*invA)(psi, chi_tmp);

      /* Step (ii) */
      /* psi_e = A_ee^-1 * [chi_e  -  D_eo * psi_o] */
      {
	T tmp1, tmp2;

	A->evenOddLinOp(tmp1, psi, PLUS);
	tmp2[rb[0]] = chi - tmp1;
	A->evenEvenInvLinOp(psi, tmp2, PLUS);
      }
  
      // Compute residual
      {
	T  r;
	A->unprecLinOp(r, psi, PLUS);
	r -= chi;
	res.resid = sqrt(norm2(r));
      }

      END_CODE();

      return res;
    }

  private:
    // Hide default constructor
    PrecFermActQprop() {}

    Handle< EvenOddPrecLinearOperator<T,P,Q> > A;
    Handle< LinOpSystemSolver<T> > invA;
  };


  typedef LatticeFermion LF;
  typedef multi1d<LatticeColorMatrix> LCM;


  template<>
  SystemSolver<LF>* 
  EvenOddPrecWilsonTypeFermAct<LF,LCM,LCM>::qprop(Handle< FermState<LF,LCM,LCM> > state,
						  const GroupXML_t& invParam) const
  {
    return new PrecFermActQprop<LF,LCM,LCM>(Handle< EvenOddPrecLinearOperator<LF,LCM,LCM> >((*this).linOp(state)), 
					    Handle< LinOpSystemSolver<LF> >((*this).invLinOp(state,invParam)));
  }
  
} // namespace Chroma 
