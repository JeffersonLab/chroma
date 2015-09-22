// -*- C++ -*-
/*! \file
 *  \brief Even-odd preconditioned Clover fermion action
 */

#ifndef __prec_clover_back_fermact_w_h__
#define __prec_clover_back_fermact_w_h__

#include "eoprec_logdet_wilstype_fermact_w.h"
#include "actions/ferm/linop/lgherm_w.h"
#include "actions/ferm/fermacts/clover_back_fermact_params_w.h"
#include "actions/ferm/fermacts/eoprec_clover_fermact_w.h"

namespace Chroma 
{ 
  //! Name and registration
  /*! \ingroup fermacts */
  namespace EvenOddPrecCloverBackFermActEnv
  {
    extern const std::string name;
    bool registerAll();
  }
  

  //! Even-odd preconditioned Clover fermion action
  /*! \ingroup fermacts
   *
   * Even-odd preconditioned clover fermion action. 
   * Only defined on odd subset.
   */

  class EvenOddPrecCloverBackFermAct : public  EvenOddPrecCloverFermAct
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Partial constructor
    EvenOddPrecCloverBackFermAct():EvenOddPrecCloverFermAct() {}

    //! General FermState
    EvenOddPrecCloverBackFermAct(Handle< CreateFermState<T,P,Q> > cfs_, 
			     const CloverBackFermActParams& param_) : 
      cfs(cfs_), param(param_) {}

    //! Copy constructor
    EvenOddPrecCloverBackFermAct(const EvenOddPrecCloverBackFermAct& a) : 
      cfs(a.cfs), param(a.param) {}

    //! Produce a linear operator for this action
    EvenOddPrecLogDetLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const;

    //! Produce the gamma_5 hermitian operator H_w
    LinearOperator<LatticeFermion>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const 
      { 
	return new lgherm<LatticeFermion>(linOp(state));
      }

    //! Destructor is automatic
    ~EvenOddPrecCloverBackFermAct() {}

  protected:
    //! Return the fermion BC object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

    //! Assignment
    void operator=(const EvenOddPrecCloverBackFermAct& a) {}

  private:
    Handle< CreateFermState<T,P,Q> >  cfs;
    CloverBackFermActParams param;
  };

}; // End Namespace Chroma


#endif
