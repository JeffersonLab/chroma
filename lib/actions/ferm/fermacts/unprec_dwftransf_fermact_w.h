// -*- C++ -*-
// $Id: unprec_dwftransf_fermact_w.h,v 1.1 2004-11-02 10:33:49 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion action
 */

#ifndef __unprec_dwftransf_fermact_w_h__
#define __unprec_dwftransf_fermact_w_h__

#include "fermact.h"
#include "invtype.h"
#include "actions/ferm/linop/lgherm_w.h"
#include "actions/ferm/linop/unprec_dwftransf_linop_w.h"
#include "io/param_io.h"       // to get AnisoParam_t

using namespace QDP;
using namespace Chroma;

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
  void write(XMLReader& xml, const string& path, const UnprecDWFTransfFermActParams& param);


  //! Unpreconditioned DWFTransf fermion action
  /*! \ingroup fermact
   *
   * Supports creation and application for fermion actions
   */
  class UnprecDWFTransfFermAct : public UnprecWilsonTypeFermAct<LatticeFermion>
  {
  public:
    //! General FermBC
    UnprecDWFTransfFermAct(Handle< FermBC<LatticeFermion> > fbc_, 
			   const Real& Mass_, 
			   const Real& b5_,
			   const Real& c5_,
			   const InvertParam_t& invParam_) : 
      fbc(fbc_) {
      param.Mass=Mass_; 
      param.b5=b5_;
      param.c5=c5_;
      param.invParam=invParam_;
    }

    //! General FermBC
    UnprecDWFTransfFermAct(Handle< FermBC<LatticeFermion> > fbc_, 
			const UnprecDWFTransfFermActParams& param_) : 
      fbc(fbc_), param(param_) {}

    //! Copy constructor
    UnprecDWFTransfFermAct(const UnprecDWFTransfFermAct& a) : 
      fbc(a.fbc), param(a.param) {}

    //! Assignment
    UnprecDWFTransfFermAct& operator=(const UnprecDWFTransfFermAct& a)
      {fbc=a.fbc; param=a.param; return *this;}

    //! Return the fermion BC object for this action
    const FermBC<LatticeFermion>& getFermBC() const {return *fbc;}

    //! Produce a linear operator for this action
    const LinearOperator<LatticeFermion>* linOp(Handle<const ConnectState> state) const;

    //! Produce a linear operator M^dag.M for this action
    const LinearOperator<LatticeFermion>* lMdagM(Handle<const ConnectState> state) const;

    //! Produce the gamma_5 hermitian operator H_w

    const LinearOperator<LatticeFermion>* gamma5HermLinOp(Handle< const ConnectState> state) const { 
      return new lgherm<LatticeFermion>(linOp(state));
    }

    //! Compute dS_f/dU
    void dsdu(multi1d<LatticeColorMatrix>& result,
	      Handle<const ConnectState> state,
	      const LatticeFermion& psi) const;

    //! Destructor is automatic
    ~UnprecDWFTransfFermAct() {}

  private:
    UnprecDWFTransfFermAct() {} //hide default constructor
   
  private:
    Handle< FermBC<LatticeFermion> >  fbc;
    UnprecDWFTransfFermActParams param;
  };

}

#endif
