// -*- C++ -*-
// $Id: prec_parwilson_fermact_w.h,v 2.2 2006-01-17 16:01:46 bjoo Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson fermion action with parity breaking term
 */

#ifndef __prec_parwilson_fermact_w_h__
#define __prec_parwilson_fermact_w_h__

#include "fermact.h"


namespace Chroma
{
  //! Name and registration
  namespace EvenOddPrecParWilsonFermActEnv
  {
    extern const std::string name;
    extern const bool registered;
  }


  //! Params for wilson ferm acts
  struct EvenOddPrecParWilsonFermActParams
  {
    EvenOddPrecParWilsonFermActParams() {}
    EvenOddPrecParWilsonFermActParams(XMLReader& in, const std::string& path);
    
    Real Mass;
    Real H;
  };


  // Reader/writers
  void read(XMLReader& xml, const string& path, EvenOddPrecParWilsonFermActParams& param);
  void write(XMLWriter& xml, const string& path, const EvenOddPrecParWilsonFermActParams& param);


  //! Even-odd preconditioned Wilson fermion action with parity breaking term
  /*! \ingroup fermacts
   *
   * Even-odd preconditioned wilson fermion action with parity breaking term
   * Only defined on odd subset.
   *
   * The kernel for Wilson fermions with a parity breaking term is
   *
   *      M  =  (d+M) + i*H*gamma_5  - (1/2) D'
   */
  class EvenOddPrecParWilsonFermAct : public EvenOddPrecConstDetWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! General FermBC
    EvenOddPrecParWilsonFermAct(Handle< FermBC<LatticeFermion> > fbc_, 
				const EvenOddPrecParWilsonFermActParams& param_) : 
      fbc(fbc_), param(param_) {}

    //! Copy constructor
    EvenOddPrecParWilsonFermAct(const EvenOddPrecParWilsonFermAct& a) : 
      fbc(a.fbc), param(a.param) {}

    //! Assignment
    EvenOddPrecParWilsonFermAct& operator=(const EvenOddPrecParWilsonFermAct& a)
      {fbc=a.fbc; param=a.param; return *this;}

    //! Return the fermion BC object for this action
    const FermBC<LatticeFermion>& getFermBC() const {return *fbc;}

    //! Produce a linear operator for this action
    const EvenOddPrecConstDetLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >* linOp(Handle<const ConnectState> state) const;

    //! Produce the gamma_5 hermitin op gamma_5 M
    const LinearOperator<LatticeFermion>* hermitianLinOp(Handle<const ConnectState> state) const {
      QDP_error_exit("gamma5HermLinOp not implemented yet for this action\n");
      return 0;
    }

    //! Destructor is automatic
    ~EvenOddPrecParWilsonFermAct() {}

  private:
    Handle< FermBC<LatticeFermion> >  fbc;
    EvenOddPrecParWilsonFermActParams param;
  };

}


#endif
