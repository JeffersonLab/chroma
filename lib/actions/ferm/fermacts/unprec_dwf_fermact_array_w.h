// -*- C++ -*-
// $Id: unprec_dwf_fermact_array_w.h,v 3.6 2006-10-19 16:01:29 edwards Exp $
/*! \file
 *  \brief Unpreconditioned domain-wall fermion action
 */

#ifndef __unprec_dwf_fermact_array_w_h__
#define __unprec_dwf_fermact_array_w_h__

#include "actions/ferm/fermacts/unprec_dwf_fermact_base_array_w.h"
#include "io/aniso_io.h"
#include "actions/ferm/linop/lg5Rherm_w.h"

namespace Chroma
{
  //! Name and registration
  namespace UnprecDWFermActArrayEnv
  {
    extern const std::string name;
    bool registerAll();
  }
  

  //! Params for DWF
  struct UnprecDWFermActArrayParams
  {
    UnprecDWFermActArrayParams() {}
    UnprecDWFermActArrayParams(XMLReader& in, const std::string& path);
    
    Real OverMass;
    Real Mass;
    Real a5;
    int  N5;
    AnisoParam_t anisoParam;
  };


  // Reader/writers
  void read(XMLReader& xml, const string& path, UnprecDWFermActArrayParams& param);
  void write(XMLWriter& xml, const string& path, const UnprecDWFermActArrayParams& param);



  //! Unpreconditioned domain-wall fermion action
  /*! \ingroup fermacts
   *
   * Unprecondition domain-wall fermion action. The conventions used here
   * are specified in Phys.Rev.D63:094505,2001 (hep-lat/0005002).
   */

  class UnprecDWFermActArray : public UnprecDWFermActBaseArray<LatticeFermion, 
				   multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General FermBC
    UnprecDWFermActArray(Handle< CreateFermState<T,P,Q> > cfs_, 
			 const Real& OverMass_, const Real& Mass_, int N5_) : 
      cfs(cfs_)
      {
	param.a5=1;
	param.OverMass = OverMass_;
	param.Mass = Mass_;
	param.N5 = N5_;
      }

    //! General FermBC
    UnprecDWFermActArray(Handle< CreateFermState<T,P,Q> > cfs_, 
			 const UnprecDWFermActArrayParams& p) :
      cfs(cfs_), param(p) {}

    //! Copy constructor
    UnprecDWFermActArray(const UnprecDWFermActArray& a) : 
      cfs(a.cfs), param(a.param) {}

    //! Length of DW flavor index/space
    int size() const {return param.N5;}

    //! Return the quark mass
    Real getQuarkMass() const {return param.Mass;}

    //! Produce an unpreconditioned linear operator for this action with arbitrary quark mass
    UnprecDWLikeLinOpBaseArray<T,P,Q>* unprecLinOp(Handle< FermState<T,P,Q> > state, 
						   const Real& m_q) const;

    //! Destructor is automatic
    ~UnprecDWFermActArray() {}

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
    void quarkProp(LatticePropagator& q_sol,   // Oops, need to make propagator type more general
		   XMLWriter& xml_out,
		   const LatticePropagator& q_src,
		   int t_src, int j_decay,
		   Handle< FermState<T,P,Q> > state,
		   const GroupXML_t& invParam,
		   QuarkSpinType quarkSpinType,
		   bool obsvP,
		   int& ncg_had) const;


    //! Produce a hermitian version of the linear operator
    /*! This code is generic */
    LinearOperatorArray<T>* hermitianLinOp(Handle< FermState<T,P,Q> > state) const
    {
      return new lg5RHermArray<T>(linOp(state));
    }

  protected:
    //! Return the fermion BC object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

    //! Partial constructor
    UnprecDWFermActArray() {}
    //! Hide =
    void operator=(const UnprecDWFermActArray& a) {}

  private:
    Handle< CreateFermState<T,P,Q> >  cfs;
    UnprecDWFermActArrayParams param;
  };

}

#endif
