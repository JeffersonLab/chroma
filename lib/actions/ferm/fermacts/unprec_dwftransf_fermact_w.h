// -*- C++ -*-
// $Id: unprec_dwftransf_fermact_w.h,v 3.3 2006-10-19 16:01:29 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion action
 */

#ifndef __unprec_dwftransf_fermact_w_h__
#define __unprec_dwftransf_fermact_w_h__

#include "unprec_wilstype_fermact_w.h"
#include "actions/ferm/linop/unprec_dwftransf_linop_w.h"
#include "io/aniso_io.h"
#include "actions/ferm/invert/syssolver_cg_params.h"

namespace Chroma
{
  //! Name and registration
  namespace UnprecDWFTransfFermActEnv
  {
    extern const std::string name;
    bool registerAll();
  }


  //! Params for wilson ferm acts
  struct UnprecDWFTransfFermActParams
  {
    UnprecDWFTransfFermActParams() {}
    UnprecDWFTransfFermActParams(XMLReader& in, const std::string& path);
    
    Real Mass;
    Real b5;
    Real c5;
    SysSolverCGParams  invParam;
  };


  // Reader/writers
  void read(XMLReader& xml, const string& path, UnprecDWFTransfFermActParams& param);
  void write(XMLWriter& xml, const string& path, const UnprecDWFTransfFermActParams& param);


  //! Unpreconditioned DWFTransf fermion action
  /*! \ingroup fermacts
   *
   * Supports creation and application for fermion actions
   */
  class UnprecDWFTransfFermAct : public UnprecWilsonTypeFermAct<LatticeFermion, 
				 multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General FermBC
    UnprecDWFTransfFermAct(Handle< CreateFermState<T,P,Q> > cfs_, 
			const UnprecDWFTransfFermActParams& param_) : 
      cfs(cfs_), param(param_) {}

    //! Copy constructor
    UnprecDWFTransfFermAct(const UnprecDWFTransfFermAct& a) : 
      cfs(a.cfs), param(a.param) {}

    //! Produce a linear operator for this action
    UnprecLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const;

    //! Produce the gamma_5 hermitian operator H_w
    //  Actually, this operator is already Hermitian, so just return
    //  linop here... It is just a beastly hack... Maybe gamma5Herm
    //  should just be renamed Herm...
    LinearOperator<T>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const 
      { 
	return linOp(state);
      }

    //! Destructor is automatic
    ~UnprecDWFTransfFermAct() {}

  protected:
    //! Return the fermion BC object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

  private:
    UnprecDWFTransfFermAct() {} //hide default constructor
    //! Hide =
    void operator=(const UnprecDWFTransfFermAct& a) {}
   
  private:
    Handle< CreateFermState<T,P,Q> >  cfs;
    UnprecDWFTransfFermActParams param;
  };

}

#endif
