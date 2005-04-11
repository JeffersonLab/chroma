// -*- C++ -*-
// $Id: unprec_wilson_fermact_w.h,v 1.25 2005-04-11 01:59:59 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion action
 */

#ifndef __unprec_wilson_fermact_w_h__
#define __unprec_wilson_fermact_w_h__

#include "fermact.h"
#include "actions/ferm/linop/lgherm_w.h"
#include "io/aniso_io.h"


namespace Chroma
{
  //! Name and registration
  namespace UnprecWilsonFermActEnv
  {
    extern const std::string name;
    extern const bool registered;
  }


  //! Params for wilson ferm acts
  struct UnprecWilsonFermActParams
  {
    UnprecWilsonFermActParams() {}
    UnprecWilsonFermActParams(XMLReader& in, const std::string& path);
    
    Real Mass;
    AnisoParam_t anisoParam;
  };


  // Reader/writers
  void read(XMLReader& xml, const string& path, UnprecWilsonFermActParams& param);
  void write(XMLWriter& xml, const string& path, const UnprecWilsonFermActParams& param);


  //! Unpreconditioned Wilson fermion action
  /*! \ingroup fermacts
   *
   * Supports creation and application for fermion actions
   */
  class UnprecWilsonFermAct : public UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! General FermBC
    UnprecWilsonFermAct(Handle< FermBC<LatticeFermion> > fbc_, 
			const Real& Mass_) : 
      fbc(fbc_) {param.Mass=Mass_;}

    //! General FermBC
    UnprecWilsonFermAct(Handle< FermBC<LatticeFermion> > fbc_, 
			const UnprecWilsonFermActParams& param_) : 
      fbc(fbc_), param(param_) {}

    //! Copy constructor
    UnprecWilsonFermAct(const UnprecWilsonFermAct& a) : 
      fbc(a.fbc), param(a.param) {}

    //! Assignment
    UnprecWilsonFermAct& operator=(const UnprecWilsonFermAct& a)
      {fbc=a.fbc; param=a.param; return *this;}

    //! Return the fermion BC object for this action
    const FermBC<LatticeFermion>& getFermBC() const {return *fbc;}

    //! Produce a linear operator for this action
    const UnprecLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >* linOp(Handle<const ConnectState> state) const;

    //! Produce a linear operator M^dag.M for this action
    const LinearOperator<LatticeFermion>* lMdagM(Handle<const ConnectState> state) const;

    //! Produce the gamma_5 hermitian operator H_w

    const LinearOperator<LatticeFermion>* gamma5HermLinOp(Handle< const ConnectState> state) const { 
      return new lgherm<LatticeFermion>(linOp(state));
    }

    //! Destructor is automatic
    ~UnprecWilsonFermAct() {}

  private:
    UnprecWilsonFermAct() {} //hide default constructor
   
  private:
    Handle< FermBC<LatticeFermion> >  fbc;
    UnprecWilsonFermActParams param;
  };

}

#endif
