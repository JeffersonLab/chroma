// -*- C++ -*-
// $Id: unprec_clover_fermact_w.h,v 1.8 2004-12-12 21:22:15 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Clover fermion action
 */

#ifndef __unprec_clover_fermact_w_h__
#define __unprec_clover_fermact_w_h__

#include "fermact.h"
#include "actions/ferm/linop/lgherm_w.h"

using namespace QDP;

namespace Chroma
{
  //! Name and registration
  namespace UnprecCloverFermActEnv
  {
    extern const std::string name;
    extern const bool registered;
  }


  //! Params for wilson ferm acts
  struct UnprecCloverFermActParams
  {
    UnprecCloverFermActParams(XMLReader& in, const std::string& path);
    
    Real Mass;
    Real ClovCoeff;
    Real u0;
  };


  // Reader/writers
  void read(XMLReader& xml, const string& path, UnprecCloverFermActParams& param);
  void write(XMLReader& xml, const string& path, const UnprecCloverFermActParams& param);


  //! Unpreconditioned Clover fermion action
  /*! \ingroup fermact
   *
   * Unpreconditioned clover fermion action
   */
  class UnprecCloverFermAct : public UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! General FermBC
    /*! Isotropic action */
    UnprecCloverFermAct(Handle< FermBC<LatticeFermion> > fbc_, 
			const UnprecCloverFermActParams& param_) : 
      fbc(fbc_), param(param_) {}

    //! Copy constructor
    UnprecCloverFermAct(const UnprecCloverFermAct& a) : 
      fbc(a.fbc), param(a.param) {}

    //! Assignment
    UnprecCloverFermAct& operator=(const UnprecCloverFermAct& a)
      {fbc=a.fbc; param=a.param; return *this;}

    //! Return the fermion BC object for this action
    const FermBC<LatticeFermion>& getFermBC() const {return *fbc;}

    //! Produce a linear operator for this action
    const UnprecLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >* linOp(Handle<const ConnectState> state) const;

    //! Produce a linear operator M^dag.M for this action
    const LinearOperator<LatticeFermion>* lMdagM(Handle<const ConnectState> state) const;

    const LinearOperator<LatticeFermion>* gamma5HermLinOp(Handle< const ConnectState> state) const { 
      return new lgherm<LatticeFermion>(linOp(state));
    }

    //! Destructor is automatic
    ~UnprecCloverFermAct() {}

  private:
    // Hide partial constructor
    UnprecCloverFermAct() {}

  private:
    Handle< FermBC<LatticeFermion> >  fbc;
    UnprecCloverFermActParams param;
  };

}

#endif
