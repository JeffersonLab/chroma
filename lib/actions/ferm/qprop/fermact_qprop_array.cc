// $Id: fermact_qprop_array.cc,v 2.0 2005-09-25 21:04:30 edwards Exp $
/*! \file
 *  \brief Propagator solver for a generic non-preconditioned fermion operator
 *
 *  Solve for the propagator of a generic non-preconditioned fermion operator
 */

#include "fermact.h"
#include "actions/ferm/invert/invcg2_array.h"

namespace Chroma 
{
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
  class FermAct5DQprop : public SystemSolver< multi1d<T> >
  {
  public:
    //! Constructor
    /*!
     * \param A_         Linear operator ( Read )
     * \param invParam_  inverter parameters ( Read )
     */
    FermAct5DQprop(Handle< const LinearOperator< multi1d<T> > > A_,
		   const InvertParam_t& invParam_) : A(A_), invParam(invParam_) 
    {}

    //! Destructor is automatic
    ~FermAct5DQprop() {}

    //! Expected length of array index
    int size() const {return A->size();}

    //! Return the subset on which the operator acts
    const OrderedSubset& subset() const {return all;}

    //! Solver the linear system
    /*!
     * \param psi      quark propagator ( Modify )
     * \param chi      source ( Read )
     * \return number of CG iterations
     */
    int operator() (multi1d<T>& psi, const multi1d<T>& chi) const
    {
      START_CODE();

      int n_count;
      QDPIO::cout << "Inv param type " << invParam.invType << endl;

      if (psi.size() != size() && chi.size() != size())
	QDP_error_exit("FA5DQprop: sizes wrong");

      switch(invParam.invType)
      {
      case CG_INVERTER: 
      {
	/* chi_1 = M_dag(u) * chi_1 */
	multi1d<T> tmp(size());
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
    FermAct5DQprop() {}

    Handle< const LinearOperator< multi1d<T> > > A;
    const InvertParam_t invParam;
  };


  template<>
  const SystemSolver< multi1d<LatticeFermion> >* 
  FermAct5D<LatticeFermion>::qpropT(Handle<const ConnectState> state,
				    const InvertParam_t& invParam) const
  {
    Handle< const LinearOperator< multi1d<LatticeFermion> > > A(linOp(state));

    return new FermAct5DQprop<LatticeFermion>(A, invParam);
  }

}; // Namespace Chroma
