// -*- C++ -*-
// $Id: unprec_dwftransf_fermact_w.h,v 3.0 2006-04-03 04:58:47 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion action
 */

#ifndef __unprec_dwftransf_fermact_w_h__
#define __unprec_dwftransf_fermact_w_h__

#include "fermact.h"
#include "invtype.h"
#include "actions/ferm/linop/unprec_dwftransf_linop_w.h"
#include "io/aniso_io.h"


namespace Chroma
{
  //! Name and registration
  namespace UnprecDWFTransfFermActEnv
  {
    extern const std::string name;
    extern const bool registered;
  }


  //! Params for wilson ferm acts
  struct UnprecDWFTransfFermActParams
  {
    UnprecDWFTransfFermActParams() {}
    UnprecDWFTransfFermActParams(XMLReader& in, const std::string& path);
    
    Real Mass;
    Real b5;
    Real c5;
    Chroma::InvertParam_t invParam;
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
			   const Real& Mass_, 
			   const Real& b5_,
			   const Real& c5_,
			   const InvertParam_t& invParam_) : 
      cfs(cfs_) 
      {
	param.Mass=Mass_; 
	param.b5=b5_;
	param.c5=c5_;
	param.invParam=invParam_;
      }

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
