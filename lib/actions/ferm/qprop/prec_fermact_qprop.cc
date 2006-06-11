// $Id: prec_fermact_qprop.cc,v 3.1 2006-06-11 06:30:32 edwards Exp $
/*! \file
 *  \brief Propagator solver for a generic even-odd preconditioned fermion operator
 *
 *  Solve for the propagator of a even-odd non-preconditioned fermion operator
 */

#include "fermact.h"
#include "actions/ferm/invert/invcg2.h"

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
		     const InvertParam_t& invParam_) : A(A_), invParam(invParam_) 
      {}

    //! Destructor is automatic
    ~PrecFermActQprop() {}

    //! Return the subset on which the operator acts
    const OrderedSubset& subset() const {return all;}

    //! Solver the linear system
    /*!
     * \param psi      quark propagator ( Modify )
     * \param chi      source ( Read )
     * \return number of CG iterations
     */
    SystemSolverResults_t operator() (T& psi, const T& chi) const
    {
      START_CODE();

      SystemSolverResults_t res;

      /* Step (i) */
      /* chi_tmp =  chi_o - D_oe * A_ee^-1 * chi_o */
      T chi_tmp;
      {
	T tmp1, tmp2;

	A->evenEvenInvLinOp(tmp1, chi, PLUS);
	A->oddEvenLinOp(tmp2, tmp1, PLUS);
	chi_tmp[rb[1]] = chi - tmp2;
      }

      switch(invParam.invType)
      {
      case CG_INVERTER: 
      {
	/* tmp = A_dag(u) * chi_tmp */
	T  tmp;
	(*A)(tmp, chi_tmp, MINUS);
    
	/* psi = (M^dag * M)^(-1) chi */
	InvCG2(*A, tmp, psi, invParam.RsdCG, invParam.MaxCG, res.n_count);
      }
      break;
  
#if 0
      case MR_INVERTER:
	/* psi = M^(-1) chi_tmp */
	InvMR(*A, chi_tmp, psi, 
	      invParam.MRover, 
	      invParam.RsdCG, 
	      invParam.MaxCG, res.n_count);
	break;

      case BICG_INVERTER:
	/* psi = M^(-1) chi_tmp */
	InvBiCG(*A, chi_tmp, psi, 
		invParam.RsdCG, 
		invParam.MaxCG, res.n_count);
	break;
#endif
  
      default:
	QDP_error_exit("Unknown inverter type", invParam.invType);
      }
  
      if ( res.n_count == invParam.MaxCG )
	QDP_error_exit("no convergence in the inverter", res.n_count);
  
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
    const InvertParam_t invParam;
  };


  typedef LatticeFermion LF;
  typedef multi1d<LatticeColorMatrix> LCM;


  template<>
  SystemSolver<LF>* 
  EvenOddPrecWilsonTypeFermAct<LF,LCM,LCM>::qprop(Handle< FermState<LF,LCM,LCM> > state,
						  const InvertParam_t& invParam) const
  {
    return new PrecFermActQprop<LF,LCM,LCM>(
      Handle< EvenOddPrecLinearOperator<LF,LCM,LCM> >(linOp(state)),
      invParam);
  }
  
} // namespace Chroma 
