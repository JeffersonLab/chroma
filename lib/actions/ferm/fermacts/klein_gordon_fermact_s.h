// -*- C++ -*-
// $Id: klein_gordon_fermact_s.h,v 3.1 2006-12-07 18:26:18 edwards Exp $
/*! \file
 *  \brief Klein-Gordon boson action masquerading action as a staggered action
 */

#ifndef __klein_gordon_fermact_s_h__
#define __klein_gordon_fermact_s_h__

#include "stagtype_fermact_s.h"
#include "io/aniso_io.h"

namespace Chroma 
{ 
  //! Name and registration
  namespace KleinGordonFermActEnv
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Params for klein-gordon action
  /*! \ingroup fermacts */
  struct KleinGordonFermActParams
  {
    KleinGordonFermActParams();
    KleinGordonFermActParams(XMLReader& in, const std::string& path);
    
    Real Mass;
    AnisoParam_t anisoParam;
  };


  // Reader/writers
  /*! \ingroup fermacts */
  void read(XMLReader& xml, const string& path, KleinGordonFermActParams& param);

  /*! \ingroup fermacts */
  void write(XMLWriter& xml, const string& path, const KleinGordonFermActParams& param);


  //! Klein-Gordon boson action
  /*! \ingroup fermacts
   *
   */
  class KleinGordonFermAct : public UnprecStaggeredTypeFermAct<
    LatticeStaggeredFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeStaggeredFermion      T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General CreateFermState
    KleinGordonFermAct(Handle< CreateFermState<T,P,Q> > cfs_, 
		       const KleinGordonFermActParams& p) :
      cfs(cfs_), param(p) {}
  
    //! Copy constructor
    KleinGordonFermAct(const KleinGordonFermAct& a) : 
      cfs(a.cfs), param(a.param) {}

    //! Create state should apply the BC
    FermState<T,P,Q>* createState(const Q& u_) const;

    //! Return the fermion BC object for this action
    const FermBC<T,P,Q>& getFermBC() const {return cfs->getBC();}

    //! Produce a linear operator for this action
    UnprecLinearOperator<T,P,Q>* linOp(Handle< FermState<T,P,Q> > state_) const;

    //! Produce a linear operator M^dag.M for this action
    DiffLinearOperator<T,P,Q>* lMdagM(Handle< FermState<T,P,Q> > state_) const;

    //! Return a linear operator solver for this action to solve MdagM*psi=chi 
    /*! Default implementation provided */
//    MdagMSystemSolver<T>* invMdagM(Handle< FermState<T,P,Q> > state,
//				   const GroupXML_t& invParam) const;

    //! Return a multi-shift linear operator solver for this action to solve (MdagM+shift)*psi=chi 
//    MdagMMultiSystemSolver<T>* mInvMdagM(Handle< FermState<T,P,Q> > state,
//					 const GroupXML_t& invParam) const;

    //! Return quark prop solver, solution of unpreconditioned system
    /*! This is already supplied in  chroma/lib/actions/ferm/qprop/fermact_qprop.cc */
//    SystemSolver<T>* qprop(Handle< FermState<T,P,Q> > state,
//			   const GroupXML_t& invParam) const;

    //! accessors 
    const Real getQuarkMass() const {return param.Mass;}

    //! Destructor is automatic
    ~KleinGordonFermAct() {}

  protected:
    //! Return the fermion BC object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

  private:
    KleinGordonFermAct() {} //hide default constructor
    void operator=(const KleinGordonFermAct& a) {} // Assignment
  
  private:
    Handle< CreateFermState<T,P,Q> >  cfs;
    KleinGordonFermActParams  param;
  };


}; // End Namespace Chroma

#endif
