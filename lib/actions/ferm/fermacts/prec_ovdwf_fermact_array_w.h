// -*- C++ -*-
// $Id: prec_ovdwf_fermact_array_w.h,v 1.9 2005-01-14 20:13:04 edwards Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned Overlap-DWF (Borici) action
 */

#ifndef __prec_ovdwf_fermact_array_w_h__
#define __prec_ovdwf_fermact_array_w_h__

#include "fermact.h"
#include "actions/ferm/fermacts/prec_dwf_fermact_base_array_w.h"


namespace Chroma
{
  //! Name and registration
  namespace EvenOddPrecOvDWFermActArrayEnv
  {
    extern const std::string name;
    extern const bool registered;
  }
  

  //! Params for DWF
  struct EvenOddPrecOvDWFermActArrayParams
  {
    EvenOddPrecOvDWFermActArrayParams(XMLReader& in, const std::string& path);
    
    Real OverMass;
    Real Mass;
    Real a5;
    int  N5;
  };


  // Reader/writers
  void read(XMLReader& xml, const string& path, EvenOddPrecOvDWFermActArrayParams& param);
  void write(XMLWriter& xml, const string& path, const EvenOddPrecOvDWFermActArrayParams& param);


  //! 4D style even-odd preconditioned Overlap-DWF (Borici) action
  /*! \ingroup fermact
   *
   * 4D style even-odd preconditioned Overlap-DWF (Borici) action
   * Follows notes of Brower (10/2003)
   *
   * Hopefully, the conventions used here
   * are specified in Phys.Rev.D63:094505,2001 (hep-lat/0005002).
   */
  class EvenOddPrecOvDWFermActArray : public EvenOddPrecDWFermActBaseArray< LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! General FermBC
    EvenOddPrecOvDWFermActArray(Handle< FermBC< multi1d<LatticeFermion> > > fbc_, 
				const Real& OverMass_, const Real& Mass_, int N5_) : 
      fbc(fbc_), OverMass(OverMass_), Mass(Mass_), N5(N5_) 
      {
	a5=1;
	QDPIO::cout << "Construct EvenOddPrecOvDWFermActArray: OverMass = " << OverMass 
		    << "  Mass = " << Mass 
		    << "  N5 = " << N5 
		    << "  a5 = " << a5 
		    << endl;
      }

    //! General FermBC
    EvenOddPrecOvDWFermActArray(Handle< FermBC< multi1d<LatticeFermion> > > fbc_, 
				const EvenOddPrecOvDWFermActArrayParams& param) :
      fbc(fbc_), OverMass(param.OverMass), Mass(param.Mass), a5(param.a5), N5(param.N5) {}

    //! Copy constructor
    EvenOddPrecOvDWFermActArray(const EvenOddPrecOvDWFermActArray& a) : 
      fbc(a.fbc), OverMass(a.OverMass), Mass(a.Mass), a5(a.a5), N5(a.N5) {}

    //! Assignment
    EvenOddPrecOvDWFermActArray& operator=(const EvenOddPrecOvDWFermActArray& a)
      {fbc=a.fbc; OverMass=a.OverMass; Mass=a.Mass; a5=a.a5; N5=a.N5; return *this;}

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
    ~EvenOddPrecOvDWFermActArray() {}

  private:
    Handle< FermBC< multi1d<LatticeFermion> > >  fbc;
    Real OverMass;
    Real Mass;
    Real a5;
    int  N5;
  };

}

#endif
