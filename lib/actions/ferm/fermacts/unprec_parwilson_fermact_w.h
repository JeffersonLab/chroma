// -*- C++ -*-
// $Id: unprec_parwilson_fermact_w.h,v 1.4 2004-12-12 21:22:15 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion action with parity breaking term
 */

#ifndef __unprec_parwilson_fermact_w_h__
#define __unprec_parwilson_fermact_w_h__

#include "fermact.h"

using namespace QDP;

namespace Chroma
{
  //! Name and registration
  namespace UnprecParWilsonFermActEnv
  {
    extern const std::string name;
    extern const bool registered;
  }


  //! Params for wilson ferm acts
  struct UnprecParWilsonFermActParams
  {
    UnprecParWilsonFermActParams() {}
    UnprecParWilsonFermActParams(XMLReader& in, const std::string& path);
    
    Real Mass;
    Real H;
  };


  // Reader/writers
  void read(XMLReader& xml, const string& path, UnprecParWilsonFermActParams& param);
  void write(XMLReader& xml, const string& path, const UnprecParWilsonFermActParams& param);


  //! Unpreconditioned Wilson fermion action with parity breaking term
  /*! \ingroup fermact
   *
   * Supports creation and application for fermion actions
   *
   * The kernel for Wilson fermions with a parity breaking term is
   *
   *      M  =  (d+M) + i*H*gamma_5  - (1/2) D'
   */
  class UnprecParWilsonFermAct : public UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! General FermBC
    UnprecParWilsonFermAct(Handle< FermBC<LatticeFermion> > fbc_, 
			   const UnprecParWilsonFermActParams& param_) : 
      fbc(fbc_), param(param_) {}

    //! Copy constructor
    UnprecParWilsonFermAct(const UnprecParWilsonFermAct& a) : 
      fbc(a.fbc), param(a.param) {}

    //! Assignment
    UnprecParWilsonFermAct& operator=(const UnprecParWilsonFermAct& a)
      {fbc=a.fbc; param=a.param; return *this;}

    //! Return the fermion BC object for this action
    const FermBC<LatticeFermion>& getFermBC() const {return *fbc;}

    //! Produce a linear operator for this action
    const UnprecLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >* linOp(Handle<const ConnectState> state) const;

    //! Produce a linear operator M^dag.M for this action
    const LinearOperator<LatticeFermion>* lMdagM(Handle<const ConnectState> state) const;

    //! Produce the gamma_5 hermitin op gamma_5 M
    const LinearOperator<LatticeFermion>* gamma5HermLinOp(Handle<const ConnectState> state) const {

      QDP_error_exit("gamma5HermLinOp not implemented yet for this action\n");
      return 0;
    }

    //! Destructor is automatic
    ~UnprecParWilsonFermAct() {}

  private:
    UnprecParWilsonFermAct() {} //hide default constructor
  
  private:
    Handle< FermBC<LatticeFermion> >  fbc;
    UnprecParWilsonFermActParams param;
  };

}

using namespace Chroma;

#endif
