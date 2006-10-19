// -*- C++ -*-
// $Id: unprec_parwilson_fermact_w.h,v 3.2 2006-10-19 16:01:29 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion action with parity breaking term
 */

#ifndef __unprec_parwilson_fermact_w_h__
#define __unprec_parwilson_fermact_w_h__

#include "unprec_wilstype_fermact_w.h"


namespace Chroma
{
  //! Name and registration
  namespace UnprecParWilsonFermActEnv
  {
    extern const std::string name;
    bool registerAll();
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
  void write(XMLWriter& xml, const string& path, const UnprecParWilsonFermActParams& param);


  //! Unpreconditioned Wilson fermion action with parity breaking term
  /*! \ingroup fermacts
   *
   * Supports creation and application for fermion actions
   *
   * The kernel for Wilson fermions with a parity breaking term is
   *
   *      M  =  (d+M) + i*H*gamma_5  - (1/2) D'
   */
  class UnprecParWilsonFermAct : public UnprecWilsonTypeFermAct<LatticeFermion, 
				 multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General FermBC
    UnprecParWilsonFermAct(Handle< CreateFermState<T,P,Q> > cfs_, 
			   const UnprecParWilsonFermActParams& param_) : 
      cfs(cfs_), param(param_) {}

    //! Copy constructor
    UnprecParWilsonFermAct(const UnprecParWilsonFermAct& a) : 
      cfs(a.cfs), param(a.param) {}

    //! Produce a linear operator for this action
    UnprecLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const;

    //! Produce the gamma_5 hermitin op gamma_5 M
    LinearOperator<T>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const 
      {
	QDP_error_exit("hermitianLinOp not implemented yet for this action\n");
	return 0;
      }

    //! Destructor is automatic
    ~UnprecParWilsonFermAct() {}

  protected:
    //! Return the fermion create state for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

  private:
    UnprecParWilsonFermAct() {} //hide default constructor
    void operator=(const UnprecParWilsonFermAct& a) {}   // Hide =
  
  private:
    Handle< CreateFermState<T,P,Q> >  cfs;
    UnprecParWilsonFermActParams param;
  };

}


#endif
