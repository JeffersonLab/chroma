// -*- C++ -*-
// $Id: lg5eps_w.h,v 3.2 2007-02-22 21:11:46 bjoo Exp $
/*! \file
 *  \brief Internal pole epsilon operator. Just the unitary part
 *
 *   ie: gamma_5 eps(H_w)
 */

#ifndef __lg5eps_w_h__
#define __lg5eps_w_h__

#include "linearop.h"
#include "unprec_wilstype_fermact_w.h" 



namespace Chroma 
{ 
  //! Internal Overlap-pole operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   *
   *   Chi  =   gamma_5 * B . Psi 
   *  where  B  is the pole approx. to eps(H(m)) 
   *
   * Internally, it computes
   *   Chi  =    gamma_5 * B . Psi 
   * and then rescales at the end to the correct normalization
   *
   *  NOTE: B is hermitian, so       
   *     (1 + gamma_5 * B)^dag = (1 + B * gamma_5) 
   *                           = gamma_5 * (1 + gamma_5 * B) * gamma_5 
   */

  class lg5eps : public LinearOperator<LatticeFermion>
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Creation routine
    /*!
     * \param _MdagM          M^dag.M of underlying linop M      (Read)
     * \param _M              Underlying linop M	               (Read)
     * \param _numroot 	    number of poles in expansion       (Read)
     * \param _constP         constant coeff                     (Read)
     * \param _resP           numerator                          (Read)
     * \param _rootQ          denom                              (Read)
     * \param _OperEigVec     eigenvectors      	               (Read)
     * \param _EigValFunc     eigenvalues      	               (Read)
     * \param _NEig           number of eigenvalues              (Read)
     * \param _MaxCG          MaxCG inner CG                     (Read)
     * \param _RsdCG          residual for inner CG              (Read)
     */
    lg5eps(const UnprecWilsonTypeFermAct<T,P,Q>& S_aux,
	   Handle< FermState<T,P,Q> > state,
	   int _numroot, 
	   const Real& _constP, 
	   const multi1d<Real>& _resP,
	   const multi1d<Real>& _rootQ, 
	   int _NEig,
	   const multi1d<Real>& _EigValFunc,
	   const multi1d<LatticeFermion>& _EigVec,
	   int _MaxCG,
	   const Real& _RsdCG,
	   const int _ReorthFreq ) :
      M(S_aux.linOp(state)), MdagM(S_aux.lMdagM(state)), 
      numroot(_numroot), constP(_constP),
      resP(_resP), rootQ(_rootQ), EigVec(_EigVec), EigValFunc(_EigValFunc),
      NEig(_NEig), MaxCG(_MaxCG), RsdCG(_RsdCG),  ReorthFreq(_ReorthFreq) {}

    //! Destructor is automatic
    ~lg5eps() {}
 
    //! Only defined on the entire lattice
    const Subset& subset() const {return all;}

    //! Apply the operator onto a source vector
    void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;

    //! Apply the operator onto a source vector. Here epsilon is the
    //  RsdCG for the multi shift solve in evaluating the sign function.
    void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign, Real epsilon) const;

  private:
    Handle<const LinearOperator<LatticeFermion> > M;
    Handle<const LinearOperator<LatticeFermion> > MdagM;

    // Copy all of these rather than reference them.
    int numroot;
    const Real constP;
    const multi1d<Real> resP;
    const multi1d<Real> rootQ;
    const multi1d<LatticeFermion> EigVec;
    const multi1d<Real> EigValFunc;
    int NEig;
    int MaxCG;
    const Real RsdCG;
    const int   ReorthFreq;
  };


} // End Namespace Chroma


#endif
