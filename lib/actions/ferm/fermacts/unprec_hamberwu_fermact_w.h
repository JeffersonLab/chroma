// -*- C++ -*-
// $Id: unprec_hamberwu_fermact_w.h,v 2.1 2005-12-03 04:20:20 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Hamber-Wu fermion action
 */

#ifndef __unprec_hamberwu_fermact_w_h__
#define __unprec_hamberwu_fermact_w_h__

#include "fermact.h"
#include "actions/ferm/linop/lgherm_w.h"
#include "io/aniso_io.h"


namespace Chroma
{
  //! Name and registration
  namespace UnprecHamberWuFermActEnv
  {
    extern const std::string name;
    extern const bool registered;
  }


  //! Params for wilson ferm acts
  struct UnprecHamberWuFermActParams
  {
    UnprecHamberWuFermActParams() {}
    UnprecHamberWuFermActParams(XMLReader& in, const std::string& path);
    
    Real Mass;
    Real u0;
    AnisoParam_t anisoParam;
  };


  // Reader/writers
  void read(XMLReader& xml, const string& path, UnprecHamberWuFermActParams& param);
  void write(XMLWriter& xml, const string& path, const UnprecHamberWuFermActParams& param);


  //! Unpreconditioned HamberWu fermion action
  /*! \ingroup fermacts
   *
   * Supports creation and application for fermion actions
   */
  class UnprecHamberWuFermAct : public UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! General FermBC
    UnprecHamberWuFermAct(Handle< FermBC<LatticeFermion> > fbc_, 
			  const UnprecHamberWuFermActParams& param_) : 
      fbc(fbc_), param(param_) {}

    //! Copy constructor
    UnprecHamberWuFermAct(const UnprecHamberWuFermAct& a) : 
      fbc(a.fbc), param(a.param) {}

    //! Assignment
    UnprecHamberWuFermAct& operator=(const UnprecHamberWuFermAct& a)
      {fbc=a.fbc; param=a.param; return *this;}

    //! Return the fermion BC object for this action
    const FermBC<LatticeFermion>& getFermBC() const {return *fbc;}

    //! Produce a linear operator for this action
    const UnprecLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >* linOp(Handle<const ConnectState> state) const;

    //! Produce a linear operator M^dag.M for this action
    const LinearOperator<LatticeFermion>* lMdagM(Handle<const ConnectState> state) const;

    //! Produce the gamma_5 hermitian operator H_w
    const LinearOperator<LatticeFermion>* hermitianLinOp(Handle< const ConnectState> state) const { 
      return new lgherm<LatticeFermion>(linOp(state));
    }

    //! Destructor is automatic
    ~UnprecHamberWuFermAct() {}

  private:
    UnprecHamberWuFermAct() {} //hide default constructor
   
  private:
    Handle< FermBC<LatticeFermion> >  fbc;
    UnprecHamberWuFermActParams param;
  };

}

#endif
