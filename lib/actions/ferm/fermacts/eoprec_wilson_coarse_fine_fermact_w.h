// -*- C++ -*-
// $Id: eoprec_wilson_coarse_fine_fermact_w.h,v 3.1 2008-01-21 20:18:50 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson fermion action supporting 2+2 anisotropy
 */

#ifndef __prec_wilson_coarse_fine_fermact_w_h__1
#define __prec_wilson_coarse_fine_fermact_w_h__

#include "eoprec_constdet_wilstype_fermact_w.h"
#include "actions/ferm/fermacts/wilson_coarse_fine_fermact_params_w.h"
#include "actions/ferm/linop/lgherm_w.h"

namespace Chroma
{
  //! Name and registration
  /*! \ingroup fermacts */
  namespace EvenOddPrecWilsonCoarseFineFermActEnv
  {
    extern const std::string name;
    bool registerAll();
  }
  

  //! Even-odd preconditioned WilsonCoarseFine fermion action
  /*! \ingroup fermacts
   *
   * Even-odd preconditioned wilson fermion action. 
   * Only defined on odd subset.
   */
  class EvenOddPrecWilsonCoarseFineFermAct : 
    public EvenOddPrecConstDetWilsonTypeFermAct<LatticeFermion, 
    multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General FermBC with Anisotropy
    EvenOddPrecWilsonCoarseFineFermAct(Handle< CreateFermState<T,P,Q> > cfs_, 
				       const WilsonCoarseFineFermActParams& p) :
      cfs(cfs_) {init(p);}

    //! Copy constructor
    EvenOddPrecWilsonCoarseFineFermAct(const EvenOddPrecWilsonCoarseFineFermAct& a) : 
      cfs(a.cfs), param(a.param) {}

    //! Produce a linear operator for this action
    EvenOddPrecConstDetLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const;

    //! Produce the gamma_5 hermitian operator H_w
    LinearOperator<T>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const 
      { 
	return new lgherm<LatticeFermion>(linOp(state));
      }

    //! Return a linear operator solver for this action to solve M*psi=chi 
    LinOpSystemSolver<T>* invLinOp(Handle< FermState<T,P,Q> > state,
				   const GroupXML_t& invParam) const;

    //! Return a linear operator solver for this action to solve MdagM*psi=chi 
    MdagMSystemSolver<T>* invMdagM(Handle< FermState<T,P,Q> > state,
				   const GroupXML_t& invParam) const;

    //! Destructor is automatic
    ~EvenOddPrecWilsonCoarseFineFermAct() {}

    Double getQuarkMass(void) const { 
      return param.Mass;
    }

  protected:
    //! Init params
    void init(const WilsonCoarseFineFermActParams& param_);

    //! Return the fermion BC object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

    //! Partial constructor
    EvenOddPrecWilsonCoarseFineFermAct() {}
    //! Assignment
    void operator=(const EvenOddPrecWilsonCoarseFineFermAct& a) {}

  private:
    struct Param_t
    {
      Real          Mass;
      multi1d<Real> coeffs;
    };

    Handle< CreateFermState<T,P,Q> >  cfs;
    Param_t param;
  };

}

#endif
