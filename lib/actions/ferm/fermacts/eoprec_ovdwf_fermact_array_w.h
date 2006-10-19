// -*- C++ -*-
// $Id: eoprec_ovdwf_fermact_array_w.h,v 3.1 2006-10-19 16:01:27 edwards Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned Overlap-DWF (Borici) action
 */

#ifndef __prec_ovdwf_fermact_array_w_h__
#define __prec_ovdwf_fermact_array_w_h__

#include "actions/ferm/fermacts/eoprec_dwf_fermact_base_array_w.h"


namespace Chroma
{
  //! Name and registration
  namespace EvenOddPrecOvDWFermActArrayEnv
  {
    extern const std::string name;
    bool registerAll();
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
  /*! \ingroup fermacts
   *
   * 4D style even-odd preconditioned Overlap-DWF (Borici) action
   * Follows notes of Brower (10/2003)
   *
   * Hopefully, the conventions used here
   * are specified in Phys.Rev.D63:094505,2001 (hep-lat/0005002).
   */
  class EvenOddPrecOvDWFermActArray : public EvenOddPrecDWFermActBaseArray<LatticeFermion, 
				   multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General FermBC
    EvenOddPrecOvDWFermActArray(Handle< CreateFermState<T,P,Q> > cfs_, 
				const Real& OverMass_, const Real& Mass_, int N5_) : 
      cfs(cfs_), OverMass(OverMass_), Mass(Mass_), N5(N5_) 
      {
	a5=1;
	QDPIO::cout << "Construct EvenOddPrecOvDWFermActArray: OverMass = " << OverMass 
		    << "  Mass = " << Mass 
		    << "  N5 = " << N5 
		    << "  a5 = " << a5 
		    << endl;
      }

    //! General FermBC
    EvenOddPrecOvDWFermActArray(Handle< CreateFermState<T,P,Q> > cfs_, 
				const EvenOddPrecOvDWFermActArrayParams& param) :
      cfs(cfs_), OverMass(param.OverMass), Mass(param.Mass), a5(param.a5), N5(param.N5) {}

    //! Copy constructor
    EvenOddPrecOvDWFermActArray(const EvenOddPrecOvDWFermActArray& a) : 
      cfs(a.cfs), OverMass(a.OverMass), Mass(a.Mass), a5(a.a5), N5(a.N5) {}

    //! Length of DW flavor index/space
    int size() const {return N5;}

    //! Return the quark mass
    Real getQuarkMass() const {return Mass;}

    //! Produce an unpreconditioned linear operator for this action with arbitrary quark mass
    UnprecDWLikeLinOpBaseArray<T,P,Q>* unprecLinOp(Handle< FermState<T,P,Q> > state, 
						   const Real& m_q) const;

    //! Produce an even-odd preconditioned linear operator for this action with arbitrary quark mass
    EvenOddPrecDWLikeLinOpBaseArray<T,P,Q>* precLinOp(Handle< FermState<T,P,Q> > state, 
						      const Real& m_q) const;

    //! Destructor is automatic
    ~EvenOddPrecOvDWFermActArray() {}

  protected:
    //! Return the fermion BC object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

    //! Partial constructor
    EvenOddPrecOvDWFermActArray() {}
    //! Hide =
    void operator=(const EvenOddPrecOvDWFermActArray& a) {}

  private:
    Handle< CreateFermState<T,P,Q> >  cfs;
    Real OverMass;
    Real Mass;
    Real a5;
    int  N5;
  };

}

#endif
