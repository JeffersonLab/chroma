// -*- C++ -*-
// $Id: prec_wilson_fermact_w.h,v 1.8 2004-09-08 02:48:25 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson fermion action
 */

#ifndef __prec_wilson_fermact_w_h__
#define __prec_wilson_fermact_w_h__

#include "fermact.h"
#include "actions/ferm/linop/lgherm_w.h"
#include "io/param_io.h"       // to get AnisoParam_t

using namespace QDP;

namespace Chroma
{
  //! Name and registration
  namespace EvenOddPrecWilsonFermActEnv
  {
    extern const std::string name;
    extern const bool registered;
    //! Name to be used
  }
  

  //! Params for wilson ferm acts
  struct EvenOddPrecWilsonFermActParams
  {
    EvenOddPrecWilsonFermActParams() {}
    EvenOddPrecWilsonFermActParams(XMLReader& in, const std::string& path);
    
    Real Mass;
    AnisoParam_t anisoParam;
  };


  // Reader/writers
  void read(XMLReader& xml, const string& path, EvenOddPrecWilsonFermActParams& param);
  void write(XMLReader& xml, const string& path, const EvenOddPrecWilsonFermActParams& param);


  //! Even-odd preconditioned Wilson fermion action
  /*! \ingroup fermact
   *
   * Even-odd preconditioned wilson fermion action. 
   * Only defined on odd subset.
   */
  class EvenOddPrecWilsonFermAct : public EvenOddPrecWilsonTypeFermAct<LatticeFermion>
  {
  public:
    //! General FermBC
    EvenOddPrecWilsonFermAct(Handle< FermBC<LatticeFermion> > fbc_, 
			     const Real& Mass_) : 
      fbc(fbc_) {param.Mass=Mass_; initHeader(param.anisoParam);}

    //! General FermBC with Anisotropy
    EvenOddPrecWilsonFermAct(Handle< FermBC<LatticeFermion> > fbc_, 
			     const EvenOddPrecWilsonFermActParams& param_) :
      fbc(fbc_), param(param_) {}

    //! Copy constructor
    EvenOddPrecWilsonFermAct(const EvenOddPrecWilsonFermAct& a) : 
      fbc(a.fbc), param(a.param) {}

    //! Assignment
    EvenOddPrecWilsonFermAct& operator=(const EvenOddPrecWilsonFermAct& a)
      {fbc=a.fbc; param=a.param; return *this;}

    //! Return the fermion BC object for this action
    const FermBC<LatticeFermion>& getFermBC() const {return *fbc;}

    //! Produce a linear operator for this action
    const EvenOddPrecLinearOperator<LatticeFermion>* linOp(Handle<const ConnectState> state) const;

    //! Produce a linear operator M^dag.M for this action
    const LinearOperator<LatticeFermion>* lMdagM(Handle<const ConnectState> state) const;

    //! Produce the gamma_5 hermitian operator H_w
    const LinearOperator<LatticeFermion>* gamma5HermLinOp(Handle< const ConnectState> state) const { 
      return new lgherm<LatticeFermion>(linOp(state));
    }

    //! Override - compute dS_f/dU
    void dsdu(multi1d<LatticeColorMatrix>& result,
	      Handle<const ConnectState> state,
	      const LatticeFermion& psi) const;

    //! Destructor is automatic
    ~EvenOddPrecWilsonFermAct() {}

  private:
    Handle< FermBC<LatticeFermion> >  fbc;
    EvenOddPrecWilsonFermActParams param;
  };

}

#endif
