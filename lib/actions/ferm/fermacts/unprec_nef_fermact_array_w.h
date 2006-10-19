// -*- C++ -*-
// $Id: unprec_nef_fermact_array_w.h,v 3.6 2006-10-19 16:01:29 edwards Exp $
/*! \file
 *  \brief Unpreconditioned NEF domain-wall fermion action
 */

#ifndef __unprec_nef_fermact_array_w_h__
#define __unprec_nef_fermact_array_w_h__

#include "unprec_wilstype_fermact_w.h"
#include "actions/ferm/fermacts/unprec_dwf_fermact_base_array_w.h"


namespace Chroma
{
  //! Name and registration
  namespace UnprecNEFFermActArrayEnv
  {
    extern const std::string name;
    bool registerAll();
  }
  

  //! Params for NEFF
  struct UnprecNEFFermActArrayParams
  {
    UnprecNEFFermActArrayParams() {}
    UnprecNEFFermActArrayParams(XMLReader& in, const std::string& path);
    
    Real OverMass;
    Real b5;
    Real c5;
    Real Mass;
    int  N5;
  };


  // Reader/writers
  void read(XMLReader& xml, const string& path, UnprecNEFFermActArrayParams& param);
  void write(XMLWriter& xml, const string& path, const UnprecNEFFermActArrayParams& param);


  //! Unpreconditioned NEF fermion action
  /*! \ingroup fermacts
   *
   * Unprecondition NEF fermion action. The conventions used here
   * are specified in Phys.Rev.D63:094505,2001 (hep-lat/0005002).
   * See also Brower et.al. LATTICE04
   */
  class UnprecNEFFermActArray : public UnprecDWFermActBaseArray<LatticeFermion, 
				   multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General FermBC
    UnprecNEFFermActArray(Handle< CreateFermState<T,P,Q> > cfs_, 
			  const Real& OverMass_, 
			  const Real& b5_, const Real& c5_, 
			  const Real& Mass_, int N5_) : 
      cfs(cfs_), OverMass(OverMass_), b5(b5_),c5(c5_), Mass(Mass_), N5(N5_) {}

    //! General FermBC
    UnprecNEFFermActArray(Handle< CreateFermState<T,P,Q> > cfs_, 
			  const UnprecNEFFermActArrayParams& param) :
      cfs(cfs_), OverMass(param.OverMass), b5(param.b5), c5(param.c5), Mass(param.Mass), N5(param.N5) {}

    //! Copy constructor
    UnprecNEFFermActArray(const UnprecNEFFermActArray& a) : 
      cfs(a.cfs), OverMass(a.OverMass), b5(a.b5), c5(a.c5), 
      Mass(a.Mass), N5(a.N5) {}

    //! Length of DW flavor index/space
    int size() const {return N5;}

    //! Return the quark mass
    Real getQuarkMass() const {return Mass;}

    //! Produce an unpreconditioned linear operator for this action with arbitrary quark mass
    UnprecDWLikeLinOpBaseArray<T,P,Q>* unprecLinOp(Handle< FermState<T,P,Q> > state, 
						   const Real& m_q) const;

    //! Destructor is automatic
    ~UnprecNEFFermActArray() {}

    //! Given a complete propagator as a source, this does all the inversions needed
    /*!
     * This routine is actually generic to Domain Wall fermions (Array) fermions
     *
     * \param q_sol        quark propagator ( Write )
     * \param q_src        source ( Read )
     * \param xml_out      diagnostic output ( Modify )
     * \param state        gauge connection state ( Read )
     * \param t_src        time slice of source ( Read )
     * \param j_decay      direction of decay ( Read )
     * \param invParam     inverter parameters ( Read )
     * \param ncg_had      number of CG iterations ( Write )
     */
    void quarkProp(LatticePropagator& q_sol,
		   XMLWriter& xml_out,
		   const LatticePropagator& q_src,
		   int t_src, int j_decay,
		   Handle< FermState<T,P,Q> > state,
		   const GroupXML_t& invParam,
		   QuarkSpinType quarkSpinType,
		   bool obsvP,
		   int& ncg_had) const;
      
  protected:
    //! Return the fermion BC object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

  private:
    //! Partial constructor
    UnprecNEFFermActArray() {}
    //! Hide =
    void operator=(const UnprecNEFFermActArray& a) {}

  private:
    Handle< CreateFermState<T,P,Q> >  cfs;
    Real OverMass;
    Real b5;
    Real c5;
    Real Mass;
    int  N5;
  };

}

#endif
