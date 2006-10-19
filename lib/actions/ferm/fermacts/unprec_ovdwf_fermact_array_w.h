// -*- C++ -*-
// $Id: unprec_ovdwf_fermact_array_w.h,v 3.2 2006-10-19 16:01:29 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Overlap-DWF (Borici) action
 */

#ifndef __unprec_ovdwf_fermact_array_w_h__
#define __unprec_ovdwf_fermact_array_w_h__

#include "unprec_wilstype_fermact_w.h"
#include "actions/ferm/fermacts/unprec_dwf_fermact_base_array_w.h"


namespace Chroma
{
  //! Name and registration
  namespace UnprecOvDWFermActArrayEnv
  {
    extern const std::string name;
    bool registerAll();
  }
  

  //! Params for OvDWF
  struct UnprecOvDWFermActArrayParams
  {
    UnprecOvDWFermActArrayParams() {}
    UnprecOvDWFermActArrayParams(XMLReader& in, const std::string& path);
    
    Real OverMass;
    Real Mass;
    Real a5;
    int  N5;
  };


  // Reader/writers
  void read(XMLReader& xml, const string& path, UnprecOvDWFermActArrayParams& param);
  void write(XMLWriter& xml, const string& path, const UnprecOvDWFermActArrayParams& param);



  //! Unpreconditioned Overlap-style (Borici) OvDWF fermion action
  /*! \ingroup fermacts
   *
   * Unprecondition domain-wall fermion action. The conventions used here
   * are specified in Phys.Rev.D63:094505,2001 (hep-lat/0005002).
   */
  class UnprecOvDWFermActArray : public UnprecDWFermActBaseArray<LatticeFermion, 
				   multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General FermBC
    UnprecOvDWFermActArray(Handle< CreateFermState<T,P,Q> > cfs_, 
			   const Real& OverMass_, const Real& Mass_, int N5_) : 
      cfs(cfs_), OverMass(OverMass_), Mass(Mass_), N5(N5_) {a5=1;}

    //! General FermBC
    UnprecOvDWFermActArray(Handle< CreateFermState<T,P,Q> > cfs_, 
			   const UnprecOvDWFermActArrayParams& param) :
      cfs(cfs_), OverMass(param.OverMass), Mass(param.Mass), a5(param.a5), N5(param.N5) {}

    //! Copy constructor
    UnprecOvDWFermActArray(const UnprecOvDWFermActArray& a) : 
      cfs(a.cfs), OverMass(a.OverMass), Mass(a.Mass), a5(a.a5), N5(a.N5) {}

    //! Length of DW flavor index/space
    int size() const {return N5;}

    //! Return the quark mass
    Real getQuarkMass() const {return Mass;}

    //! Produce an unpreconditioned linear operator for this action with arbitrary quark mass
    UnprecDWLikeLinOpBaseArray<T,P,Q>* unprecLinOp(Handle< FermState<T,P,Q> > state, 
						   const Real& m_q) const;

    //! Destructor is automatic
    ~UnprecOvDWFermActArray() {}

  protected:
    //! Return the fermion BC object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

  private:
    UnprecOvDWFermActArray() {} //! General FermBC
    void operator=(const UnprecOvDWFermActArray& a) {} //! Hide =

  private:
    Handle< CreateFermState<T,P,Q> >  cfs;
    Real OverMass;
    Real Mass;
    Real a5;
    int  N5;
  };

}

#endif
