// -*- C++ -*-
// $Id: prec_dwf_fermact_array_w.h,v 1.9 2004-09-19 02:39:45 edwards Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned domain-wall fermion action
 */

#ifndef __prec_dwf_fermact_array_w_h__
#define __prec_dwf_fermact_array_w_h__

#include "fermact.h"
#include "actions/ferm/fermacts/prec_dwf_fermact_base_array_w.h"

using namespace QDP;

namespace Chroma
{
  //! Name and registration
  namespace EvenOddPrecDWFermActArrayEnv
  {
    extern const std::string name;
    extern const bool registered;
  }
  

  //! Params for DWF
  struct EvenOddPrecDWFermActArrayParams
  {
    EvenOddPrecDWFermActArrayParams(XMLReader& in, const std::string& path);
    
    Real OverMass;
    Real Mass;
    Real a5;
    int  N5;
  };


  // Reader/writers
  void read(XMLReader& xml, const string& path, EvenOddPrecDWFermActArrayParams& param);
  void write(XMLReader& xml, const string& path, const EvenOddPrecDWFermActArrayParams& param);


  //! 4D style even-odd preconditioned domain-wall fermion action
  /*! \ingroup fermact
   *
   * 4D style even-odd preconditioned domain-wall fermion action. 
   * Follows notes of Orginos (10/2003)
   *
   * Hopefully, the conventions used here
   * are specified in Phys.Rev.D63:094505,2001 (hep-lat/0005002).
   */
  
  class EvenOddPrecDWFermActArray : public EvenOddPrecDWFermActBaseArray<LatticeFermion>
  {
  public:
    //! General FermBC
    EvenOddPrecDWFermActArray(Handle< FermBC< multi1d<LatticeFermion> > > fbc_, 
			      const Real& OverMass_, const Real& Mass_, int N5_) : 
      fbc(fbc_), OverMass(OverMass_), Mass(Mass_), N5(N5_) 
      {
	a5=1;
	QDPIO::cout << "Construct EvenOddPrecDWFermActArray: OverMass = " << OverMass 
		    << "  Mass = " << Mass 
		    << "  N5 = " << N5 
		    << "  a5 = " << a5 
		    << endl;
      }

    //! General FermBC
    EvenOddPrecDWFermActArray(Handle< FermBC< multi1d<LatticeFermion> > > fbc_, 
			      const EvenOddPrecDWFermActArrayParams& param) :
      fbc(fbc_), OverMass(param.OverMass), Mass(param.Mass), a5(param.a5), N5(param.N5) {}

    //! Copy constructor
    EvenOddPrecDWFermActArray(const EvenOddPrecDWFermActArray& a) : 
      fbc(a.fbc), OverMass(a.OverMass), Mass(a.Mass), a5(a.a5), N5(a.N5) {}

    //! Assignment
    EvenOddPrecDWFermActArray& operator=(const EvenOddPrecDWFermActArray& a)
      {fbc=a.fbc; OverMass=a.OverMass; Mass=a.Mass; a5=a.a5; N5=a.N5; return *this;}

    //! Return the fermion BC object for this action
    const FermBC< multi1d<LatticeFermion> >& getFermBC() const {return *fbc;}

    //! Length of DW flavor index/space
    int size() const {return N5;}

    //! Return the quark mass
    Real quark_mass() const {return Mass;}

    //! Produce a linear operator for this action
    const EvenOddPrecDWLinOpBaseArray<LatticeFermion>* linOp(Handle<const ConnectState> state) const;

    //! Produce a linear operator M^dag.M for this action
    const LinearOperator< multi1d<LatticeFermion> >* lMdagM(Handle<const ConnectState> state) const;

    //! Produce a linear operator for this action but with quark mass 1
    const UnprecDWLinOpBaseArray<LatticeFermion>* linOpPV(Handle<const ConnectState> state) const;

    //! Destructor is automatic
    ~EvenOddPrecDWFermActArray() {}

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


  private:
    Handle< FermBC< multi1d<LatticeFermion> > >  fbc;
    
    Real OverMass;
    Real Mass;
    Real a5;
    int  N5;
  };

}

#endif
