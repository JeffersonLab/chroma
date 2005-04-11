// -*- C++ -*-
// $Id: prec_dwf_fermact_array_w.h,v 1.18 2005-04-11 01:15:07 edwards Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned domain-wall fermion action
 */

#ifndef __prec_dwf_fermact_array_w_h__
#define __prec_dwf_fermact_array_w_h__

#include "fermact.h"
#include "actions/ferm/fermacts/prec_dwf_fermact_base_array_w.h"
#include "io/aniso_io.h"

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
    EvenOddPrecDWFermActArrayParams() {}
    EvenOddPrecDWFermActArrayParams(XMLReader& in, const std::string& path);
    
    Real OverMass;
    Real Mass;
    Real a5;
    int  N5;
    AnisoParam_t anisoParam;
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
  
  class EvenOddPrecDWFermActArray : public EvenOddPrecDWFermActBaseArray< LatticeFermion, multi1d<LatticeColorMatrix> >
  {
  public:
    //! General FermBC
    EvenOddPrecDWFermActArray(Handle< FermBC< multi1d<LatticeFermion> > > fbc_, 
			      const Real& OverMass_, const Real& Mass_, int N5_) : 
      fbc(fbc_)
      {
	param.a5=1;
	param.OverMass = OverMass_;
	param.Mass = Mass_;
	param.N5 = N5_;
	QDPIO::cout << "Construct EvenOddPrecDWFermActArray: OverMass = " << param.OverMass 
		    << "  Mass = " << param.Mass 
		    << "  N5 = " << param.N5 
		    << "  a5 = " << param.a5 
		    << endl;
      }

    //! General FermBC
    EvenOddPrecDWFermActArray(Handle< FermBC< multi1d<LatticeFermion> > > fbc_, 
			      const EvenOddPrecDWFermActArrayParams& p) :
      fbc(fbc_), param(p) {}

    //! Copy constructor
    EvenOddPrecDWFermActArray(const EvenOddPrecDWFermActArray& a) : 
      fbc(a.fbc), param(a.param) {}

    //! Assignment
    EvenOddPrecDWFermActArray& operator=(const EvenOddPrecDWFermActArray& a)
      {fbc=a.fbc; param=a.param; return *this;}

    //! Return the fermion BC object for this action
    const FermBC< multi1d<LatticeFermion> >& getFermBC() const {return *fbc;}

    //! Length of DW flavor index/space
    int size() const {return param.N5;}

    //! Return the quark mass
    Real getQuarkMass() const {return param.Mass;}

    //! Produce an unpreconditioned linear operator for this action with arbitrary quark mass
    const UnprecDWLinOpBaseArray< LatticeFermion, multi1d<LatticeColorMatrix> >* unprecLinOp(Handle<const ConnectState> state, 
											     const Real& m_q) const;

    //! Produce an even-odd preconditioned linear operator for this action with arbitrary quark mass
    const EvenOddPrecDWLinOpBaseArray< LatticeFermion, multi1d<LatticeColorMatrix> >* precLinOp(Handle<const ConnectState> state, 
												const Real& m_q) const;

    //! Return possibly optimized quark prop solver, solution of preconditioned system
    const SystemSolver< multi1d<LatticeFermion> >* qpropT(Handle<const ConnectState> state,
							  const InvertParam_t& invParam) const;

    //! Destructor is automatic
    ~EvenOddPrecDWFermActArray() {}

    //! Given a complete propagator as a source, this does all the inversions needed
    /*!
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
    
    EvenOddPrecDWFermActArrayParams param;
  };

}

#endif
