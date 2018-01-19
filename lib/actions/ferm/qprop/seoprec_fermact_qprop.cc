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
