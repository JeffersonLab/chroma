/*! \file
 *  \brief Propagator solver for a generic symmetric even-odd preconditioned fermion operator
 *
 *  Solve for the propagator of a symmetric even-odd non-preconditioned fermion operator
 */

#include "seoprec_wilstype_fermact_w.h"

namespace Chroma 
{ 
  //! Propagator of a generic symmetric even-odd preconditioned fermion linear operator
  /*! \ingroup qprop
   *
   * This routine is actually generic to all symmetric even-odd preconditioned fermions
   */
  template<typename T, typename P, typename Q>
  class SymEvenOddPrecActQprop : public SystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param A_         Linear operator ( Read )
     * \param invParam_  inverter parameters ( Read )
     */
    SymEvenOddPrecActQprop(Handle< SymEvenOddPrecLinearOperator<T,P,Q> > A_,
			   Handle< LinOpSystemSolver<T> > invA_) : A(A_), invA(invA_) 
      {}

    //! Destructor is automatic
    ~SymEvenOddPrecActQprop() {}

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
      /* chi' = L^{-1} x M_diag^{-1} chi */
      T Mdiag_inv_chi;
      QDPIO::cout << "Diag Inv" << std::endl;
      A->scaleEvenEvenInvLinOp(Mdiag_inv_chi,chi, PLUS);
      A->scaleOddOddInvLinOp(Mdiag_inv_chi,chi,PLUS);

      /* Now apply L^{-1} = [ 1          0  ]  = [ 1      0 ]
       *                    [ -A^-1_oo D_oe  1 ]    [ -M_oe   1]
       */
      {
    	  T tmp;
    	  A->oddEvenLinOp(tmp, Mdiag_inv_chi,PLUS);
    	  Mdiag_inv_chi[rb[1]] -= tmp;
      }

      // Now solve
      // [   1        0   ] [x'_e] = L^{-1} M_diag^{-1} [chi_e]
      // [   0        S   ] [x'_o]                      [chi_o]
      //
      //                           = [ Mdiag_inv_chi_e ]
      //                             [ Mdiag_inv_chi_o ]
      //
      // Call inverter -- for rb[1]
      SystemSolverResults_t res = (*invA)(psi, Mdiag_inv_chi);

      // rb[0] is trivial:
      psi[rb[0]]=Mdiag_inv_chi;

      /* Step (ii) */
      /* Reconstruct solution */
      /* psi = R^{-1} psi = [ 1  - A_ee^{-1} D_eo ] [ psi_e]
       *                    [ 0              1    ] [ psi_o]
       *
       *              = [ 1  -M_eo ] [ psi_e ]
       *                [  0    1  ] [ psi_o ]
       */
      {
    	  T tmp;

    	  A->evenOddLinOp(tmp, psi, PLUS);
    	  psi[rb[0]] -= tmp;
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
    SymEvenOddPrecActQprop() {}

    Handle< SymEvenOddPrecLinearOperator<T,P,Q> > A;
    Handle< LinOpSystemSolver<T> > invA;
  };


  typedef LatticeFermion LF;
  typedef multi1d<LatticeColorMatrix> LCM;


  template<>
  SystemSolver<LF>* 
  SymEvenOddPrecWilsonTypeFermAct<LF,LCM,LCM>::qprop(Handle< FermState<LF,LCM,LCM> > state,
						     const GroupXML_t& invParam) const
  {
    StopWatch swatch2;
    
    QDPIO::cout << "  ... constructing linop " ;
    swatch2.reset(); swatch2.start();
    Handle< SymEvenOddPrecLinearOperator<LF,LCM,LCM> > lh( (*this).linOp(state) );
    swatch2.stop();
    QDPIO::cout << " ..." << swatch2.getTimeInSeconds() << " sec" << std::endl;
  
  
    QDPIO::cout << "  ... constructing invLinOp " ;
    swatch2.reset(); swatch2.start();
    Handle< LinOpSystemSolver<LF> > ilh((*this).invLinOp(state,invParam));
    swatch2.stop(); 
    QDPIO::cout << " ..." << swatch2.getTimeInSeconds() << " sec" << std::endl;
  
    QDPIO::cout << "  ... constructing SymEvenOddPrecActQprop " ;
    swatch2.reset(); swatch2.start();
    SymEvenOddPrecActQprop<LF,LCM,LCM>* ret_val = new SymEvenOddPrecActQprop<LF,LCM,LCM>(lh , ilh);
     swatch2.stop();
    QDPIO::cout << " ..." << swatch2.getTimeInSeconds() << " sec" << std::endl;

    return ret_val;
  
}
  
} // namespace Chroma 
