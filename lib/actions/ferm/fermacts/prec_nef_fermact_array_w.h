// -*- C++ -*-
// $Id: prec_nef_fermact_array_w.h,v 1.11 2005-01-02 05:21:09 edwards Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned NEF fermion action
 */

#ifndef __prec_nef_fermact_array_w_h__
#define __prec_nef_fermact_array_w_h__

#include "fermact.h"
#include "actions/ferm/fermacts/prec_dwf_fermact_base_array_w.h"

using namespace QDP;

namespace Chroma
{
  //! Name and registration
  namespace EvenOddPrecNEFFermActArrayEnv
  {
    extern const std::string name;
    extern const bool registered;
    //! Name to be used
  }
  

  //! Params for NEFF
  struct EvenOddPrecNEFFermActArrayParams
  {
    EvenOddPrecNEFFermActArrayParams() {}
    EvenOddPrecNEFFermActArrayParams(XMLReader& in, const std::string& path);
    
    Real OverMass;
    Real b5;
    Real c5;
    Real Mass;
    int  N5;
  };


  // Reader/writers
  void read(XMLReader& xml, const string& path, EvenOddPrecNEFFermActArrayParams& param);
  void write(XMLWriter& xml, const string& path, const EvenOddPrecNEFFermActArrayParams& param);


  //! 4D style even-odd preconditioned domain-wall fermion action
  /*! \ingroup fermact
   *
   * 4D style even-odd preconditioned domain-wall fermion action. 
   * Follows notes of Orginos (10/2003)
   *
   * Hopefully, the conventions used here
   * are specified in Phys.Rev.D63:094505,2001 (hep-lat/0005002).
   */
  class EvenOddPrecNEFFermActArray : public EvenOddPrecDWFermActBaseArray< LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! General FermBC
    EvenOddPrecNEFFermActArray(Handle< FermBC< multi1d<LatticeFermion> > > fbc_, 
			       const Real& OverMass_, 
			       const Real& b5_, const Real& c5_,
			       const Real& Mass_, int N5_) : 
      fbc(fbc_), OverMass(OverMass_),  b5(b5_),c5(c5_),  Mass(Mass_), N5(N5_) 
      {
	QDPIO::cout << "Construct EvenOddPrecNEFFermActArray: OverMass = " << OverMass 
		    << "  Mass = " << Mass 
		    << "  N5 = " << N5 
		    << "  b5 = " << b5
		    << "  c5 = " << c5
		    << endl;
      }

    //! General FermBC
    EvenOddPrecNEFFermActArray(Handle< FermBC< multi1d<LatticeFermion> > > fbc_, 
			       const EvenOddPrecNEFFermActArrayParams& param) :
      fbc(fbc_), OverMass(param.OverMass), b5(param.b5), c5(param.c5), Mass(param.Mass), N5(param.N5) {}

    //! Copy constructor
    EvenOddPrecNEFFermActArray(const EvenOddPrecNEFFermActArray& a) : 
      fbc(a.fbc), OverMass(a.OverMass),  b5(a.b5), c5(a.c5), 
      Mass(a.Mass), N5(a.N5) {}

    //! Assignment
    EvenOddPrecNEFFermActArray& operator=(const EvenOddPrecNEFFermActArray& a)
      {
	fbc=a.fbc; 
	OverMass=a.OverMass; 
	b5=a.b5; 
	c5=a.c5; 
	Mass=a.Mass; 
	N5=a.N5; 
	return *this;
      }

    //! Return the fermion BC object for this action
    const FermBC< multi1d<LatticeFermion> >& getFermBC() const {return *fbc;}

    //! Length of DW flavor index/space
    int size() const {return N5;}

    //! Return the quark mass
    Real getQuarkMass() const {return Mass;}

    //! Produce an unpreconditioned linear operator for this action with arbitrary quark mass
    const UnprecDWLinOpBaseArray< LatticeFermion, multi1d<LatticeColorMatrix> >* unprecLinOp(Handle<const ConnectState> state, 
											     const Real& m_q) const;

    //! Produce an even-odd preconditioned linear operator for this action with arbitrary quark mass
    const EvenOddPrecDWLinOpBaseArray< LatticeFermion, multi1d<LatticeColorMatrix> >* precLinOp(Handle<const ConnectState> state, 
												const Real& m_q) const;

    //! Destructor is automatic
    ~EvenOddPrecNEFFermActArray() {}

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
