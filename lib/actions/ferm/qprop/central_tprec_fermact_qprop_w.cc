// $Id: central_tprec_fermact_qprop_w.cc,v 1.2 2007-02-22 21:11:47 bjoo Exp $
/*! \file
 *  \brief Propagator solver for a generic even-odd preconditioned fermion operator
 *
 *  Solve for the propagator of a even-odd non-preconditioned fermion operator
 */

#include "unprec_s_cprec_t_wilstype_fermact_w.h"
#include "iluprec_s_cprec_t_wilstype_fermact_w.h"

namespace Chroma 
{ 
  //! Propagator for unpreconditioned space centrally preconditioned time FermAct
  /*! \ingroup qprop
   *
   */
  template<typename T, typename P, typename Q, template <typename, typename, typename> class L>
  class CentralPrecTimeFermActQprop : public SystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param A_         Linear operator ( Read )
     * \param invParam_  inverter parameters ( Read )
     */
    CentralPrecTimeFermActQprop(Handle< L<T,P,Q> > A_,
					   Handle< LinOpSystemSolver<T> > invA_) : A(A_), invA(invA_) 
      {}

    //! Destructor is automatic
    ~CentralPrecTimeFermActQprop() {}

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

      // We need to solve    C_L^{1}( 1 + C_L D_s C_R )C_R^{-1} \psi = \chi
      //
      // We do this by solving ( 1 + C_L D_s C_R ) \psi' = \chi'
      // with \chi' = C_L \chi
      //
      // and then at the end C_R^{-1] \psi = \psi' => \psi = C_R \psi'
      //
      // First compute \chi'
      T chi_prime;
      QDPIO::cout << "Preparing LinOp" << endl;
      A->cLeftLinOp(chi_prime, chi, PLUS);


      T psi_prime = zero;
      // Call inverter
      QDPIO::cout << "Solving" << endl;
      SystemSolverResults_t res = (*invA)(psi_prime, chi_prime);


      QDPIO::cout <<"Reconstructing" << endl;
      // Reconstruct psi = C_R psi_prime
      A->cRightLinOp(psi, psi_prime, PLUS);

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
    CentralPrecTimeFermActQprop() {}

    Handle< L<T,P,Q> > A;
    Handle< LinOpSystemSolver<T> > invA;
  };


  typedef LatticeFermion LF;
  typedef multi1d<LatticeColorMatrix> LCM;


  template<>
  SystemSolver<LF>* 
  UnprecSpaceCentralPrecTimeWilsonTypeFermAct<LF,LCM,LCM>::qprop(Handle< FermState<LF,LCM,LCM> > state,
								 const GroupXML_t& invParam) const
  {
    return new CentralPrecTimeFermActQprop<LF,LCM,LCM,UnprecSpaceCentralPrecTimeLinearOperator>(Handle< UnprecSpaceCentralPrecTimeLinearOperator<LF,LCM,LCM> >(linOp(state)), 
												Handle< LinOpSystemSolver<LF> >(invLinOp(state,invParam)));
  }

  template<>
  SystemSolver<LF>* 
  ILUPrecSpaceCentralPrecTimeWilsonTypeFermAct<LF,LCM,LCM>::qprop(Handle< FermState<LF,LCM,LCM> > state,
								 const GroupXML_t& invParam) const
  {
    return new CentralPrecTimeFermActQprop<LF,LCM,LCM,ILUPrecSpaceCentralPrecTimeLinearOperator>(Handle< ILUPrecSpaceCentralPrecTimeLinearOperator<LF,LCM,LCM> >(linOp(state)), 
												 Handle< LinOpSystemSolver<LF> >(invLinOp(state,invParam)));
  }
  
} // namespace Chroma 
