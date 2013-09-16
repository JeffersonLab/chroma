// $Id: central_tprec_fermact_qprop_w.cc,v 1.6 2009-03-19 17:06:43 mcneile Exp $
/*! \file
 *  \brief Propagator solver for a generic even-odd preconditioned fermion operator
 *
 *  Solve for the propagator of a even-odd non-preconditioned fermion operator
 */

#include "qdp_config.h"
#if QDP_NS == 4
#if QDP_ND == 4
#if QDP_NC == 3


#include "unprec_s_cprec_t_wilstype_fermact_w.h"
#include "iluprec_s_cprec_t_wilstype_fermact_w.h"
#include "eo3dprec_s_cprec_t_wilstype_fermact_w.h"

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
      T chi_prime=zero;
      QDPIO::cout << "Preparing LinOp" << endl;
      A->cLeftLinOp(chi_prime, chi, PLUS);


      psi=zero;
      T psi_prime = zero;

      // Call inverter
      QDPIO::cout << "Solving" << endl;
      SystemSolverResults_t res = (*invA)(psi_prime, chi_prime);


      QDPIO::cout <<"Reconstructing" << endl;
      // Reconstruct psi = C_R psi_prime
      A->cRightLinOp(psi, psi_prime, PLUS);

      // Compute residual
      {
	T  r=zero;
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

  template<typename T, typename P, typename Q, template <typename, typename, typename> class L>
  class Central2PrecTimeFermActQprop : public SystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param A_         Linear operator ( Read )
     * \param invParam_  inverter parameters ( Read )
     */
    Central2PrecTimeFermActQprop(Handle< L<T,P,Q> > A_,
					   Handle< LinOpSystemSolver<T> > invA_) : A(A_), invA(invA_) 
      {}

    //! Destructor is automatic
    ~Central2PrecTimeFermActQprop() {}

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
      T chi_prime=zero;
      QDPIO::cout << "Preparing LinOp" << endl;
      A->leftLinOp(chi_prime, chi, PLUS);


      psi=zero;
      T psi_prime = zero;

      // Call inverter
      QDPIO::cout << "Solving" << endl;
      SystemSolverResults_t res = (*invA)(psi_prime, chi_prime);


      QDPIO::cout <<"Reconstructing" << endl;
      // Reconstruct psi = C_R psi_prime
      A->rightLinOp(psi, psi_prime, PLUS);

      // Compute residual
      {
	T  r=zero;
	A->unprecLinOp(r, psi, PLUS);
	r -= chi;
	res.resid = sqrt(norm2(r));
      }

      END_CODE();

      return res;
    }

  private:
    // Hide default constructor
    Central2PrecTimeFermActQprop() {}

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
    return new CentralPrecTimeFermActQprop<LF,LCM,LCM,UnprecSpaceCentralPrecTimeLinearOperator>(Handle< UnprecSpaceCentralPrecTimeLinearOperator<LF,LCM,LCM> >((*this).linOp(state)), 
												Handle< LinOpSystemSolver<LF> >((*this).invLinOp(state,invParam)));
  }

  template<>
  SystemSolver<LF>* 
  ILUPrecSpaceCentralPrecTimeWilsonTypeFermAct<LF,LCM,LCM>::qprop(Handle< FermState<LF,LCM,LCM> > state,
								 const GroupXML_t& invParam) const
  {
    return new CentralPrecTimeFermActQprop<LF,LCM,LCM,ILUPrecSpaceCentralPrecTimeLinearOperator>(Handle< ILUPrecSpaceCentralPrecTimeLinearOperator<LF,LCM,LCM> >((*this).linOp(state)), 
												 Handle< LinOpSystemSolver<LF> >((*this).invLinOp(state,invParam)));
  }
  

  template<>
  SystemSolver<LF>* 
  ILU2PrecSpaceCentralPrecTimeWilsonTypeFermAct<LF,LCM,LCM>::qprop(Handle< FermState<LF,LCM,LCM> > state,
								 const GroupXML_t& invParam) const
  {
    return new Central2PrecTimeFermActQprop<LF,LCM,LCM,ILU2PrecSpaceCentralPrecTimeLinearOperator>(Handle< ILU2PrecSpaceCentralPrecTimeLinearOperator<LF,LCM,LCM> >((*this).linOp(state)), 
												 Handle< LinOpSystemSolver<LF> >((*this).invLinOp(state,invParam)));
  }
  

  //! Propagator for unpreconditioned space centrally preconditioned time FermAct
  /*! \ingroup qprop
   *
   */
  template<typename T, typename P, typename Q, template <typename, typename, typename> class L>
  class EO3DPrecCentralPrecTimeFermActQprop : public SystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param A_         Linear operator ( Read )
     * \param invParam_  inverter parameters ( Read )
     */
    EO3DPrecCentralPrecTimeFermActQprop(Handle< L<T,P,Q> > A_,
					   Handle< LinOpSystemSolver<T> > invA_) : A(A_), invA(invA_) 
      {}

    //! Destructor is automatic
    ~EO3DPrecCentralPrecTimeFermActQprop() {}

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
      A->cLeftLinOp(chi_prime, chi, PLUS,0);
      A->cLeftLinOp(chi_prime, chi, PLUS,1);

      // Now I want chi'' = ( 1   0             ) ( chi'_e )
      //                    ( -M_oe M_ee^{-1} 1 ) ( chi'_o ) 
      T tmp1, tmp2;
      A->evenEvenInvLinOp(tmp1, chi_prime, PLUS);
      A->oddEvenLinOp(tmp2, tmp1, PLUS);
      chi_prime[rb3[1]] -= tmp2;
      

      T psi_prime = zero;
      // Call inverter
      QDPIO::cout << "Solving" << endl;
      SystemSolverResults_t res = (*invA)(psi_prime, chi_prime);
      A->evenEvenInvLinOp(psi_prime, chi_prime, PLUS);

      QDPIO::cout <<"Reconstructing" << endl;

      // ( 1 - M_ee^{-1} M_eo )
      // ( 0          1       )
      A->evenOddLinOp(tmp1, psi_prime, PLUS);
      A->evenEvenInvLinOp(tmp2, tmp1, PLUS);
      psi_prime[rb3[0]] -= tmp2;

      // Reconstruct psi = C_R psi_prime
      A->cRightLinOp(psi, psi_prime, PLUS,0);
      A->cRightLinOp(psi, psi_prime, PLUS,1);

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
    EO3DPrecCentralPrecTimeFermActQprop() {}

    Handle< L<T,P,Q> > A;
    Handle< LinOpSystemSolver<T> > invA;
  };

  template<>
  SystemSolver<LF>* 
  EO3DPrecSpaceCentralPrecTimeConstDetWilsonTypeFermAct<LF,LCM,LCM>::qprop(Handle< FermState<LF,LCM,LCM> > state,
									   const GroupXML_t& invParam) const
  {
    return new EO3DPrecCentralPrecTimeFermActQprop<LF,LCM,LCM,EO3DPrecSpaceCentralPrecTimeLinearOperator>(Handle< EO3DPrecSpaceCentralPrecTimeLinearOperator<LF,LCM,LCM> >((*this).linOp(state)), 
												      Handle< LinOpSystemSolver<LF> >((*this).invLinOp(state,invParam)));
  }


} // namespace Chroma 


#endif
#endif
#endif

