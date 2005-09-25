// $Id: fermact_qprop.cc,v 2.0 2005-09-25 21:04:30 edwards Exp $
/*! \file
 *  \brief Propagator solver for a generic non-preconditioned fermion operator
 *
 *  Solve for the propagator of a generic non-preconditioned fermion operator
 */

#include "fermact.h"
#include "actions/ferm/invert/invcg2.h"


namespace Chroma 
{ 
  //! Propagator of a generic non-preconditioned fermion linear operator
  /*! \ingroup qprop
   *
   * This routine is actually generic to all non-preconditioned (not red/black) fermions
   *
   * Compute the lattice fermion for a generic non-red/black fermion
   * using the source in "chi" - so, the source can
   * be of any desired form. The result will appear in "psi", which on input
   * contains an initial guess for the solution.
   */
  template<typename T>
  class FermActQprop : public SystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param A_         Linear operator ( Read )
     * \param invParam_  inverter parameters ( Read )
     */
    FermActQprop(Handle< const LinearOperator<T> > A_,
		 const InvertParam_t& invParam_) : A(A_), invParam(invParam_) 
    {}

    //! Destructor is automatic
    ~FermActQprop() {}

    //! Return the subset on which the operator acts
    const OrderedSubset& subset() const {return all;}

    //! Solver the linear system
    /*!
     * \param psi      quark propagator ( Modify )
     * \param chi      source ( Read )
     * \return number of CG iterations
     */
    int operator() (T& psi, const T& chi) const
    {
      START_CODE();

      int n_count;
  
      switch(invParam.invType)
      {
      case CG_INVERTER: 
      {
	/* tmp = M_dag(u) * chi_1 */
	T  tmp;
	(*A)(tmp, chi, MINUS);
	
	/* psi = (M^dag * M)^(-1) chi */
	InvCG2 (*A, tmp, psi, invParam.RsdCG, invParam.MaxCG, n_count);
      }
      break;
  
#if 0
      case MR_INVERTER:
	/* psi = M^(-1) chi */
	InvMR (*A, chi, psi, invParam.MRover, invParam.RsdCG, invParam.MaxCG, n_count);
	break;

      case BICG_INVERTER:
	/* psi = M^(-1) chi */
	InvBiCG (*A, chi, psi, invParam.RsdCG, invParam.MaxCG, n_count);
	break;
#endif
  
      default:
	QDP_error_exit("Unknown inverter type", invParam.invType);
      }
      
      if ( n_count == invParam.MaxCG )
	QDP_error_exit("no convergence in the inverter", n_count);
  
      END_CODE();

      return n_count;
    }

  private:
    // Hide default constructor
    FermActQprop() {}

    Handle< const LinearOperator<T> > A;
    const InvertParam_t invParam;
  };


  template<>
  const SystemSolver<LatticeFermion>* 
  FermAct4D<LatticeFermion>::qprop(Handle<const ConnectState> state,
				   const InvertParam_t& invParam) const
  {
    return new FermActQprop<LatticeFermion>(Handle< const LinearOperator<LatticeFermion> >(linOp(state)), 
					    invParam);
  }


  template<>
  const SystemSolver<LatticeStaggeredFermion>* 
  FermAct4D<LatticeStaggeredFermion>::qprop(Handle<const ConnectState> state,
					    const InvertParam_t& invParam) const
  {
    return new FermActQprop<LatticeStaggeredFermion>(Handle< const LinearOperator<LatticeStaggeredFermion> >(linOp(state)), 
					    invParam);
  }

} // namespace Chroma

