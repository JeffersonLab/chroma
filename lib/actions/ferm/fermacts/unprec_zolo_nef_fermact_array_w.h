// -*- C++ -*-
// $Id: unprec_zolo_nef_fermact_array_w.h,v 1.8 2005-01-02 05:21:10 edwards Exp $
/*! \file
 *  \brief Unpreconditioned NEF domain-wall fermion action
 */

#ifndef __unprec_zolo_nef_fermact_array_w_h__
#define __unprec_zolo_nef_fermact_array_w_h__

#include "fermact.h"
#include "actions/ferm/fermacts/unprec_dwf_fermact_base_array_w.h"
#include "actions/ferm/fermacts/overlap_state.h"
#include "io/overlap_state_info.h"

using namespace QDP;

namespace Chroma
{
  //! Name and registration
  namespace UnprecZoloNEFFermActArrayEnv
  {
    extern const std::string name;
    extern const bool registered;
  }
  

  //! Params for NEFF
  struct UnprecZoloNEFFermActArrayParams
  {
    UnprecZoloNEFFermActArrayParams() {}
    UnprecZoloNEFFermActArrayParams(XMLReader& in, const std::string& path);
    
    Real OverMass;
    Real Mass;
    Real b5;
    Real c5;
    int  N5;
    CoeffType approximation_type;
  };


  // Reader/writers
  void read(XMLReader& xml, const string& path, UnprecZoloNEFFermActArrayParams& param);
  void write(XMLWriter& xml, const string& path, const UnprecZoloNEFFermActArrayParams& param);


  //! Unpreconditioned NEF fermion action
  /*! \ingroup fermact
   *
   * Unprecondition NEF fermion action. The conventions used here
   * are specified in Phys.Rev.D63:094505,2001 (hep-lat/0005002).
   * See also Brower et.al. LATTICE04
   */
  class UnprecZoloNEFFermActArray : public UnprecDWFermActBaseArray< LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! General FermBC
    UnprecZoloNEFFermActArray(Handle< FermBC< multi1d<LatticeFermion> > > fbc_, 
			      const Real& OverMass_, 
			      const Real& Mass_, 
			      const Real& b5_,
			      const Real& c5_,
			      int N5_,
			      const CoeffType& approx_type_) : 
      fbc(fbc_), OverMass(OverMass_), Mass(Mass_), b5(b5_), c5(c5_), N5(N5_), approximation_type(approx_type_) {init();}

    //! General FermBC
    UnprecZoloNEFFermActArray(Handle< FermBC< multi1d<LatticeFermion> > > fbc_, 
			      const UnprecZoloNEFFermActArrayParams& param) :
      fbc(fbc_), OverMass(param.OverMass), Mass(param.Mass), b5(param.b5), c5(param.c5), N5(param.N5), approximation_type(param.approximation_type) {init();}

    //! Copy constructor
    UnprecZoloNEFFermActArray(const UnprecZoloNEFFermActArray& a) : 
      fbc(a.fbc), OverMass(a.OverMass), Mass(a.Mass), b5(a.b5), c5(a.c5), N5(a.N5), approximation_type(a.approximation_type)  {}

    //! Assignment
    UnprecZoloNEFFermActArray& operator=(const UnprecZoloNEFFermActArray& a)
      {
	fbc=a.fbc; 
	OverMass=a.OverMass; 
	Mass=a.Mass; 
	b5=a.b5;
	c5=a.c5;
	N5=a.N5; 
	approximation_type = a.approximation_type;
	return *this;
      }

    //! Return the fermion BC object for this action
    const FermBC< multi1d<LatticeFermion> >& getFermBC() const {return *fbc;}

    //! Length of DW flavor index/space
    int size() const {return N5;}

    //! Return the quark mass
    Real getQuarkMass() const {return Mass;}

    //! Produce an unpreconditioned linear operator for this action with arbitrary quark mass
    const UnprecDWLinOpBaseArray<LatticeFermion, multi1d<LatticeColorMatrix> >* unprecLinOp(Handle<const ConnectState> state, 
											    const Real& m_q) const;

    //! Destructor is automatic
    ~UnprecZoloNEFFermActArray() {}

    //! Given a complete propagator as a source, this does all the inversions needed
    /*! \ingroup qprop
     *
     * This routine is actually generic to Domain Wall fermions (Array) fermions
     *
     * \param q_sol    quark propagator ( Write )
     * \param q_src    source ( Read )
     * \param xml_out  diagnostic output ( Modify )
     * \param state    gauge connection state ( Read )
     * \param t_src    time slice of source ( Read )
     * \param j_decay  direction of decay ( Read )
     * \param invParam inverter parameters ( Read )
     * \param ncg_had  number of CG iterations ( Write )
     */
    void quarkProp(LatticePropagator& q_sol,   // Oops, need to make propagator type more general
		   XMLWriter& xml_out,
		   const LatticePropagator& q_src,
		   int t_src, int j_decay,
		   Handle<const ConnectState> state,
		   const InvertParam_t& invParam,
		   bool nonRelProp,
		   bool obsvP,
		   int& ncg_had);
      
    //! Create OverlapConnectState from XML
    const OverlapConnectState* 
    createState(const multi1d<LatticeColorMatrix>& u, 
		XMLReader& state_info_xml,
		const string& state_info_path) const;
    
    //! Given links, create the state needed for the linear operators
    /*! Override the parent */
    
    //! Create a ConnectState with just the gauge fields
    const OverlapConnectState*
    createState(const multi1d<LatticeColorMatrix>& u_) const ;

    //! Create a ConnectState with just the gauge fields, and a lower
    //!  approximation bound
    const OverlapConnectState*
    createState(const multi1d<LatticeColorMatrix>& u_, 
		const Real& approxMin_) const ;
 
    
    //! Create a connect State with just approximation range bounds
    const OverlapConnectState*
    createState(const multi1d<LatticeColorMatrix>& u_, 
		const Real& approxMin_, 
		const Real& approxMax_) const;


    //! Create OverlapConnectState with eigenvalues/vectors 
    const OverlapConnectState*
    createState(const multi1d<LatticeColorMatrix>& u_,
		const multi1d<Real>& lambda_lo_,
		const multi1d<LatticeFermion>& evecs_lo_, 
		const Real& lambda_hi_) const;


    //! Create from OverlapStateInfo Structure
    const OverlapConnectState*
    createState(const multi1d<LatticeColorMatrix>& u_,
		const OverlapStateInfo& state_info) const;



  private:
    void init();

    void initCoeffs(multi1d<Real>& b5_arr,
		    multi1d<Real>& c5_arr,
		    Handle<const ConnectState>& state) const ;

  private:
    Handle< FermBC< multi1d<LatticeFermion> > >  fbc;
    Real OverMass;
    Real Mass;
    Real b5;
    Real c5;
    int  N5;
    CoeffType approximation_type;
  };

}

#endif
