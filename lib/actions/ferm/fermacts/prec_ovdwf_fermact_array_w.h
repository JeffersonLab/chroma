// -*- C++ -*-
// $Id: prec_ovdwf_fermact_array_w.h,v 1.2 2004-09-08 02:48:25 edwards Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned Overlap-DWF (Borici) action
 */

#ifndef __prec_ovdwf_fermact_array_w_h__
#define __prec_ovdwf_fermact_array_w_h__

#include "fermact.h"
#include "actions/ferm/fermacts/prec_dwf_fermact_base_array_w.h"

using namespace QDP;

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
    
    Real WilsonMass;
    Real m_q;
    Real a5;
    int  N5;
  };


  // Reader/writers
  void read(XMLReader& xml, const string& path, EvenOddPrecOvDWFermActArrayParams& param);
  void write(XMLReader& xml, const string& path, const EvenOddPrecOvDWFermActArrayParams& param);


  //! 4D style even-odd preconditioned Overlap-DWF (Borici) action
  /*! \ingroup fermact
   *
   * 4D style even-odd preconditioned Overlap-DWF (Borici) action
   * Follows notes of Brower (10/2003)
   *
   * Hopefully, the conventions used here
   * are specified in Phys.Rev.D63:094505,2001 (hep-lat/0005002).
   */
  class EvenOddPrecOvDWFermActArray : public EvenOddPrecDWFermActBaseArray<LatticeFermion>
  {
  public:
    //! General FermBC
    EvenOddPrecOvDWFermActArray(Handle< FermBC< multi1d<LatticeFermion> > > fbc_, 
				const Real& WilsonMass_, const Real& m_q_, int N5_) : 
      fbc(fbc_), WilsonMass(WilsonMass_), m_q(m_q_), N5(N5_) 
      {
	a5=1;
	QDPIO::cout << "Construct EvenOddPrecOvDWFermActArray: WilsonMass = " << WilsonMass 
		    << "  m_q = " << m_q 
		    << "  N5 = " << N5 
		    << "  a5 = " << a5 
		    << endl;
      }

    //! General FermBC
    EvenOddPrecOvDWFermActArray(Handle< FermBC< multi1d<LatticeFermion> > > fbc_, 
				const EvenOddPrecOvDWFermActArrayParams& param) :
      fbc(fbc_), WilsonMass(param.WilsonMass), m_q(param.m_q), a5(param.a5), N5(param.N5) {}

    //! Copy constructor
    EvenOddPrecOvDWFermActArray(const EvenOddPrecOvDWFermActArray& a) : 
      fbc(a.fbc), WilsonMass(a.WilsonMass), m_q(a.m_q), a5(a.a5), N5(a.N5) {}

    //! Assignment
    EvenOddPrecOvDWFermActArray& operator=(const EvenOddPrecOvDWFermActArray& a)
      {fbc=a.fbc; WilsonMass=a.WilsonMass; m_q=a.m_q; a5=a.a5; N5=a.N5; return *this;}

    //! Return the fermion BC object for this action
    const FermBC< multi1d<LatticeFermion> >& getFermBC() const {return *fbc;}

    //! Length of DW flavor index/space
    int size() const {return N5;}

    //! Return the quark mass
    Real quark_mass() const {return m_q;}

    //! Produce a linear operator for this action
    const EvenOddPrecLinearOperator< multi1d<LatticeFermion> >* linOp(Handle<const ConnectState> state) const;

    //! Produce a linear operator M^dag.M for this action
    const LinearOperator< multi1d<LatticeFermion> >* lMdagM(Handle<const ConnectState> state) const;

    //! Produce a linear operator for this action but with quark mass 1
    const LinearOperator< multi1d<LatticeFermion> >* linOpPV(Handle<const ConnectState> state) const;

    //! Destructor is automatic
    ~EvenOddPrecOvDWFermActArray() {}

  private:
    Handle< FermBC< multi1d<LatticeFermion> > >  fbc;
    Real WilsonMass;
    Real m_q;
    Real a5;
    int  N5;
  };

}

#endif
