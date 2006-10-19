// -*- C++ -*-
// $Id: eoprec_parwilson_fermact_w.h,v 3.1 2006-10-19 16:01:27 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson fermion action with parity breaking term
 */

#ifndef __prec_parwilson_fermact_w_h__
#define __prec_parwilson_fermact_w_h__

#include "eoprec_constdet_wilstype_fermact_w.h"

namespace Chroma
{
  //! Name and registration
  namespace EvenOddPrecParWilsonFermActEnv
  {
    extern const std::string name;
    bool registerAll();
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
  class EvenOddPrecParWilsonFermAct : public EvenOddPrecConstDetWilsonTypeFermAct<LatticeFermion, 
				      multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General FermBC
    EvenOddPrecParWilsonFermAct(Handle< CreateFermState<T,P,Q> > fs_, 
				const EvenOddPrecParWilsonFermActParams& param_) : 
      fs(fs_), param(param_) {}

    //! Copy constructor
    EvenOddPrecParWilsonFermAct(const EvenOddPrecParWilsonFermAct& a) : 
      fs(a.fs), param(a.param) {}

    //! Produce a linear operator for this action
    EvenOddPrecConstDetLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const;

    //! Produce the gamma_5 hermitin op gamma_5 M
    LinearOperator<T>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const 
      {
	QDP_error_exit("gamma5HermLinOp not implemented yet for this action\n");
	return 0;
      }

    //! Destructor is automatic
    ~EvenOddPrecParWilsonFermAct() {}

  protected:
    //! Return the fermion create state for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *fs;}

    //! Partial constructor
    EvenOddPrecParWilsonFermAct() {}
    //! Hide =
    void operator=(const EvenOddPrecParWilsonFermAct& a) {}

  private:
    Handle< CreateFermState<T,P,Q> >  fs;
    EvenOddPrecParWilsonFermActParams param;
  };

}


#endif
