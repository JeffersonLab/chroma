// -*- C++ -*-
// $Id: unprec_nef_fermact_array_w.h,v 1.4 2004-09-16 15:20:54 kostas Exp $
/*! \file
 *  \brief Unpreconditioned NEF domain-wall fermion action
 */

#ifndef __unprec_nef_fermact_array_w_h__
#define __unprec_nef_fermact_array_w_h__

#include "fermact.h"
#include "actions/ferm/fermacts/unprec_dwf_fermact_base_array_w.h"

using namespace QDP;

namespace Chroma
{
  //! Name and registration
  namespace UnprecNEFFermActArrayEnv
  {
    extern const std::string name;
    extern const bool registered;
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
  void write(XMLReader& xml, const string& path, const UnprecNEFFermActArrayParams& param);


  //! Unpreconditioned NEF fermion action
  /*! \ingroup fermact
   *
   * Unprecondition NEF fermion action. The conventions used here
   * are specified in Phys.Rev.D63:094505,2001 (hep-lat/0005002).
   * See also Brower et.al. LATTICE04
   */
  class UnprecNEFFermActArray : public UnprecDWFermActBaseArray<LatticeFermion>
  {
  public:
    //! General FermBC
    UnprecNEFFermActArray(Handle< FermBC< multi1d<LatticeFermion> > > fbc_, 
			  const Real& OverMass_, 
			  const Real& b5_, const Real& c5_, 
			  const Real& Mass_, int N5_) : 
      fbc(fbc_), OverMass(OverMass_), b5(b5_),c5(c5_), Mass(Mass_), N5(N5_) {}

    //! General FermBC
    UnprecNEFFermActArray(Handle< FermBC< multi1d<LatticeFermion> > > fbc_, 
			  const UnprecNEFFermActArrayParams& param) :
      fbc(fbc_), OverMass(param.OverMass), b5(param.b5), c5(param.c5), Mass(param.Mass), N5(param.N5) {}

    //! Copy constructor
    UnprecNEFFermActArray(const UnprecNEFFermActArray& a) : 
      fbc(a.fbc), OverMass(a.OverMass), b5(a.b5), c5(a.c5), 
      Mass(a.Mass), N5(a.N5) {}

    //! Assignment
    UnprecNEFFermActArray& operator=(const UnprecNEFFermActArray& a)
      {
	fbc=a.fbc; 
	OverMass=a.OverMass; 
	Mass=a.Mass; 
	b5=a.b5; 
	c5=a.c5; 
	N5=a.N5; 
	return *this;
      }

    //! Return the fermion BC object for this action
    const FermBC< multi1d<LatticeFermion> >& getFermBC() const {return *fbc;}

    //! Length of DW flavor index/space
    int size() const {return N5;}

    //! Return the quark mass
    Real quark_mass() const {return Mass;}

    //! Produce a linear operator for this action
    const LinearOperator< multi1d<LatticeFermion> >* linOp(Handle<const ConnectState> state) const;

    //! Produce a linear operator M^dag.M for this action
    const LinearOperator< multi1d<LatticeFermion> >* lMdagM(Handle<const ConnectState> state) const;

    //! Produce a hermitian version of the linear operator
    /*! This code is generic */
    const LinearOperator< multi1d<LatticeFermion> >* gamma5HermLinOp(Handle<const ConnectState> state) const
      {
	// Have not implemented this yet, but it is generic
	QDPIO::cerr << "UnprecNEFFermActBase::gamma5HermLinOp not implemented" << endl;
	QDP_abort(1);
	return 0;
      }

    //! Produce a linear operator for this action but with quark mass 1
    const LinearOperator< multi1d<LatticeFermion> >* linOpPV(Handle<const ConnectState> state) const;

    //! Destructor is automatic
    ~UnprecNEFFermActArray() {}

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
    void dwf_quarkProp4(LatticePropagator& q_sol, 
			XMLWriter& xml_out,
			const LatticePropagator& q_src,
			int t_src, int j_decay,
			Handle<const ConnectState> state,
			const InvertParam_t& invParam,
			int& ncg_had);

//! Apply the Dminus operator on a vector in Ls. See my notes ;-)
    void Dminus(multi1d<LatticeFermion>& chi,
		const multi1d<LatticeFermion>& psi,
		Handle<const ConnectState> state,
		enum PlusMinus isign);
    
    void Dminus(LatticeFermion chi,
		const LatticeFermion& psi,
		Handle<const ConnectState> state,
		enum PlusMinus isign) ;
      
  private:
    Handle< FermBC< multi1d<LatticeFermion> > >  fbc;
    Real OverMass;
    Real b5;
    Real c5;
    Real Mass;
    int  N5;
  };

}

#endif
