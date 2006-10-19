// -*- C++ -*-
// $Id: unprec_hamberwu_fermact_w.h,v 3.2 2006-10-19 16:01:29 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Hamber-Wu fermion action
 */

#ifndef __unprec_hamberwu_fermact_w_h__
#define __unprec_hamberwu_fermact_w_h__

#include "unprec_wilstype_fermact_w.h"
#include "actions/ferm/linop/lgherm_w.h"
#include "io/aniso_io.h"


namespace Chroma
{
  //! Name and registration
  namespace UnprecHamberWuFermActEnv
  {
    extern const std::string name;
    bool registerAll();
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
  class UnprecHamberWuFermAct : public UnprecWilsonTypeFermAct<LatticeFermion, 
				multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General FermBC
    UnprecHamberWuFermAct(Handle< CreateFermState<T,P,Q> > cfs_, 
			  const UnprecHamberWuFermActParams& param_) : 
      cfs(cfs_), param(param_) {}

    //! Copy constructor
    UnprecHamberWuFermAct(const UnprecHamberWuFermAct& a) : 
      cfs(a.cfs), param(a.param) {}

    //! Produce a linear operator for this action
    UnprecLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state) const;

    //! Produce the gamma_5 hermitian operator H_w
    LinearOperator<T>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const 
      { 
	return new lgherm<T>(linOp(state));
      }

    //! Destructor is automatic
    ~UnprecHamberWuFermAct() {}

  protected:
    //! Return the fermion BC object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

  private:
    UnprecHamberWuFermAct() {} //hide default constructor
    void operator=(const UnprecHamberWuFermAct& a) {} // Hide =
   
  private:
    Handle< CreateFermState<T,P,Q> >  cfs;
    UnprecHamberWuFermActParams param;
  };

}

#endif
