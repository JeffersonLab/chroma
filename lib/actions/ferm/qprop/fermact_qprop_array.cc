// $Id: fermact_qprop_array.cc,v 3.3 2007-02-22 21:11:48 bjoo Exp $
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
  class FermAct5DQprop : public SystemSolverArray<T>
  {
  public:
    //! Constructor
    /*!
     * \param A_         Linear operator ( Read )
     * \param invParam_  inverter parameters ( Read )
     */
    FermAct5DQprop(Handle< LinearOperatorArray<T> > A_,
		   Handle< LinOpSystemSolverArray<T> > invA_) : A(A_), invA(invA_) 
    {}

    //! Destructor is automatic
    ~FermAct5DQprop() {}

    //! Expected length of array index
    int size() const {return A->size();}

    //! Return the subset on which the operator acts
    const Subset& subset() const {return all;}

    //! Solver the linear system
    /*!
     * \param psi      quark propagator ( Modify )
     * \param chi      source ( Read )
     * \return number of CG iterations
     */
    SystemSolverResults_t operator() (multi1d<T>& psi, const multi1d<T>& chi) const
    {
      START_CODE();

      if (psi.size() != size() && chi.size() != size())
	QDP_error_exit("FA5DQprop: sizes wrong");

      // Call inverter
      SystemSolverResults_t res = (*invA)(psi, chi);
  
      // Compute residual
      {
	multi1d<T>  r(size());
	(*A)(r, psi, PLUS);
	r -= chi;
	res.resid = sqrt(norm2(r));
      }

      END_CODE();

      return res;
    }

  private:
    // Hide default constructor
    FermAct5DQprop() {}

    Handle< LinearOperatorArray<T> > A;
    Handle< LinOpSystemSolverArray<T> > invA;
  };


  template<>
  SystemSolverArray<LatticeFermion>* 
  FermAct5D<LatticeFermion, 
	    multi1d<LatticeColorMatrix>,
	    multi1d<LatticeColorMatrix> >::qpropT(Handle< FermState< LatticeFermion,
						  multi1d<LatticeColorMatrix>,
						  multi1d<LatticeColorMatrix> > > state,
						  const GroupXML_t& invParam) const
  {
    return new FermAct5DQprop<LatticeFermion>(linOp(state),
					      invLinOp(state,invParam));
  }

} // Namespace Chroma
