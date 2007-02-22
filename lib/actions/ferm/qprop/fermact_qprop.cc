// $Id: fermact_qprop.cc,v 3.3 2007-02-22 21:11:48 bjoo Exp $
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
    FermActQprop(Handle< LinearOperator<T> > A_,
		 Handle< SystemSolver<T> > invA_) : A(A_), invA(invA_) 
    {}

    //! Destructor is automatic
    ~FermActQprop() {}

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

      // Call inverter
      SystemSolverResults_t res = (*invA)(psi, chi);
  
      // Compute residual
      {
	T  r;
	(*A)(r, psi, PLUS);
	r -= chi;
	res.resid = sqrt(norm2(r));
      }

      END_CODE();

      return res;
    }

  private:
    // Hide default constructor
    FermActQprop() {}

    Handle< LinearOperator<T> > A;
    Handle< SystemSolver<T> > invA;
  };


  /*! \ingroup qprop */
  template<>
  SystemSolver<LatticeFermion>*
  FermAct4D<LatticeFermion,
	    multi1d<LatticeColorMatrix>,
	    multi1d<LatticeColorMatrix> >::qprop(Handle< FermState<LatticeFermion, 
						 multi1d<LatticeColorMatrix>,
						 multi1d<LatticeColorMatrix> > > state,
						 const GroupXML_t& invParam) const
  {
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    return new FermActQprop<T>(linOp(state),
			       invLinOp(state,invParam));
  }


  /*! \ingroup qprop */
  template<>
  SystemSolver<LatticeStaggeredFermion>* 
  FermAct4D<LatticeStaggeredFermion,
	    multi1d<LatticeColorMatrix>,
	    multi1d<LatticeColorMatrix> >::qprop(Handle< FermState< LatticeStaggeredFermion,
						 multi1d<LatticeColorMatrix>,
						 multi1d<LatticeColorMatrix> > > state,
						 const GroupXML_t& invParam) const
  {
    // Typedefs to save typing
    typedef LatticeStaggeredFermion      T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    return new FermActQprop<T>(linOp(state),
			       invLinOp(state,invParam));
  }


} // namespace Chroma

