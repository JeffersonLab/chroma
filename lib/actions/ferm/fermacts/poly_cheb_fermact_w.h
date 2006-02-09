// -*- C++ -*-
// $Id: poly_cheb_fermact_w.h,v 2.1 2006-02-09 22:27:01 edwards Exp $
/*! \file
 *  \brief Chebyshev polynomial fermion action
 */

#ifndef __poly_cheb_fermact_w_h__
#define __poly_cheb_fermact_w_h__

#include "polyfermact.h"
#include "actions/ferm/linop/lgherm_w.h"

namespace Chroma
{
  //! Name and registration
  namespace PolyChebFermActEnv
  {
    extern const std::string name;
    extern const bool registered;
  }


  //! Params for Chebyshev polynomial preconditioner
  struct PolyChebFermActParams
  {
    PolyChebFermActParams() {}
    PolyChebFermActParams(XMLReader& in, const std::string& path);
    
    struct PolyParams
    {
      int degree;
      Real UpperBound ;
      Real LowerBound ;
      int order ;
    } polyParams;

    std::string  AuxFermAct;      /*!< string holding fermact xml */
  };

  // Reader/writers
  void read(XMLReader& xml, const string& path, PolyChebFermActParams& param);
  void write(XMLWriter& xml, const string& path, const PolyChebFermActParams& param);


  //! Chebyshev Polynomial fermion action
  /*! \ingroup fermacts
   *
   */
  class PolyChebFermAct : public PolyWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! General FermBC
    PolyChebFermAct(Handle< FermBC<LatticeFermion> > fbc_, 
		    const PolyChebFermActParams& param_) : 
      fbc(fbc_), param(param_) {init();}

    //! Copy constructor
    PolyChebFermAct(const PolyChebFermAct& a) : 
      fbc(a.fbc), param(a.param) {}

    //! Assignment
    PolyChebFermAct& operator=(const PolyChebFermAct& a)
      {fbc=a.fbc; param=a.param; return *this;}

    //! Return the fermion BC object for this action
    const FermBC<LatticeFermion>& getFermBC() const {return *fbc;}

    //! Produce a linear operator for this action
    const DiffLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >* linOp(Handle<const ConnectState> state) const;

    //! Produce a linear operator M^dag.M for this action
    const DiffLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >* lMdagM(Handle<const ConnectState> state) const;

    //! Produce the gamma_5 hermitian operator H_w
    const LinearOperator<LatticeFermion>* hermitianLinOp(Handle< const ConnectState> state) const;

    //! Produce a linear operator for this action
    const DiffLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >* polyPrecLinOp(Handle<const ConnectState> state) const;

    //! Produce a linear operator M^dag.M for this action
    const DiffLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >* polyLinOp(Handle<const ConnectState> state) const;

    //! Destructor is automatic
    ~PolyChebFermAct() {}

  private:
    void init();
    PolyChebFermAct() {} //hide default constructor
   
  private:
    Handle< FermBC<LatticeFermion> >  fbc;
    PolyChebFermActParams param;

    // A handle for the PrecWilsonFermAct
    Handle<const WilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> > > fermact;
  };

}

#endif
