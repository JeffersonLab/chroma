// $Id: dwf_fermact_qprop_array_w.cc,v 3.4 2007-02-22 21:11:47 bjoo Exp $
/*! \file
 *  \brief Base class for unprec and even-odd preconditioned DWF qprop
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/eoprec_dwf_fermact_base_array_w.h"
#include "actions/ferm/fermacts/unprec_dwf_fermact_base_array_w.h"
#include "actions/ferm/linop/dwffld_w.h"



namespace Chroma 
{
  //! Propagator DWF linear operator
  /*! \ingroup qprop
   *
   * Propagator solver for DWF-like fermions
   */
  template<typename T>
  class DWFQprop : public SystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param qpropT_    5D solver ( Read )
     * \param PV_        Pauli-Villars linear operator ( Read )
     * \param m_q_       quark mass ( Read )
     */
    DWFQprop(Handle< SystemSolverArray<T> > qpropT_,
	     Handle< LinearOperatorArray<T> > A_,
	     Handle< LinearOperatorArray<T> > PV_,
	     const Real& m_q_) : 
      qpropT(qpropT_), A(A_), PV(PV_), m_q(m_q_) {}

    //! Destructor is automatic
    ~DWFQprop() {}

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

      SystemSolverResults_t res;
      const int N5 = PV->size();
  
      // Initialize the 5D fields
      multi1d<T> chi5(N5);
      {
	//  chi5 = (chi,0,0,0,..,0)^T
	chi5 = zero;
	chi5[0] = chi;

	// tmp5 = P . chi5
	multi1d<T> tmp5(N5);
	DwfFld(tmp5, chi5, PLUS);

	// chi5 = D5(1) . tmp5 =  D5(1) . P . (chi,0,0,..,0)^T 
	(*PV)(chi5, tmp5, PLUS);
      }

      // psi5 = (psi,0,0,0,...,0)^T
      multi1d<T> psi5(N5);
      psi5 = zero;
      psi5[0] = psi;

      // Solve  D5(m_q) . psi5 = chi5
      res = (*qpropT)(psi5, chi5);

      // Compute residual
      {
	multi1d<T>  r(N5);
	(*A)(r, psi5, PLUS);
	r -= chi5;
	res.resid = sqrt(norm2(r));
      }

      // Overall normalization
      Real ftmp1 = Real(1) / Real(1 - m_q);

      // Project out first slice after  chi5 <- P^(-1) . psi5
      DwfFld(chi5, psi5, MINUS);

      // Normalize and remove contact term
      psi = ftmp1*(chi5[0] - chi);

      END_CODE();

      return res;
    }

  private:
    // Hide default constructor
    DWFQprop() {}

    Handle< SystemSolverArray<T> > qpropT;
    Handle< LinearOperatorArray<T> > A;
    Handle< LinearOperatorArray<T> > PV;
    const Real m_q;
  };


  typedef LatticeFermion LF;
  typedef multi1d<LatticeColorMatrix> LCM;


  template<>
  SystemSolver<LF>* 
  EvenOddPrecDWFermActBaseArray<LF,LCM,LCM>::qprop(Handle< FermState<LF,LCM,LCM> > state,
						   const GroupXML_t& invParam) const
  {
    return new DWFQprop<LF>(Handle< SystemSolverArray<LF> >(this->qpropT(state,invParam)), 
			    Handle< LinearOperatorArray<LF> >(this->unprecLinOp(state,getQuarkMass())),
			    Handle< LinearOperatorArray<LF> >(this->unprecLinOp(state,Real(1))),
			    getQuarkMass());
  }
  



  template<>
  SystemSolver<LF>* 
  UnprecDWFermActBaseArray<LF,LCM,LCM>::qprop(Handle< FermState<LF,LCM,LCM> > state,
					      const GroupXML_t& invParam) const
  { 
    return new DWFQprop<LF>(Handle< SystemSolverArray<LF> >(this->qpropT(state,invParam)), 
			    Handle< LinearOperatorArray<LF> >(this->unprecLinOp(state,getQuarkMass())),
			    Handle< LinearOperatorArray<LF> >(this->unprecLinOp(state,Real(1))),
			    getQuarkMass());
  }
  

} // namespace Chroma
