// -*- C++ -*-
// $Id: eoprec_dwf_fermact_array_w.h,v 3.2 2009-07-17 19:14:46 bjoo Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned domain-wall fermion action
 */

#ifndef __prec_dwf_fermact_array_w_h__
#define __prec_dwf_fermact_array_w_h__

#include "actions/ferm/fermacts/eoprec_dwf_fermact_base_array_w.h"
#include "io/aniso_io.h"

namespace Chroma
{
  //! Name and registration
  namespace EvenOddPrecDWFermActArrayEnv
  {
    extern const std::string name;
    bool registerAll();
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
  /*! \ingroup fermacts
   *
   * 4D style even-odd preconditioned domain-wall fermion action. 
   * Follows notes of Orginos (10/2003)
   *
   * Hopefully, the conventions used here
   * are specified in Phys.Rev.D63:094505,2001 (hep-lat/0005002).
   */
  
  class EvenOddPrecDWFermActArray : public EvenOddPrecDWFermActBaseArray<LatticeFermion, 
				   multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General FermBC
    EvenOddPrecDWFermActArray(Handle< CreateFermState<T,P,Q> > cfs_, 
			      const Real& OverMass_, const Real& Mass_, int N5_) : 
      cfs(cfs_)
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
    EvenOddPrecDWFermActArray(Handle< CreateFermState<T,P,Q> > cfs_, 
			      const EvenOddPrecDWFermActArrayParams& p) :
      cfs(cfs_), param(p) {}

    //! Copy constructor
    EvenOddPrecDWFermActArray(const EvenOddPrecDWFermActArray& a) : 
      cfs(a.cfs), param(a.param) {}

    //! Length of DW flavor index/space
    int size() const {return param.N5;}

    //! Return the quark mass
    Real getQuarkMass() const {return param.Mass;}

    //! Produce an unpreconditioned linear operator for this action with arbitrary quark mass
    UnprecDWLikeLinOpBaseArray<T,P,Q>* unprecLinOp(Handle< FermState<T,P,Q> > state, 
						   const Real& m_q) const;

    //! Produce an even-odd preconditioned linear operator for this action with arbitrary quark mass
    EvenOddPrecDWLikeLinOpBaseArray<T,P,Q>* precLinOp(Handle< FermState<T,P,Q> > state, 
						      const Real& m_q) const;

    //! Return possibly optimized quark prop solver, solution of preconditioned system
    SystemSolverArray<T>* qpropT(Handle< FermState<T,P,Q> > state,
				 const GroupXML_t& invParam) const;

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
    void quarkProp(PropTypeTraits<T>::Type_t& q_sol,
		   XMLWriter& xml_out,
		   const PropTypeTraits<T>::Type_t& q_src,
		   int t_src, int j_decay,
		   Handle< FermState<T,P,Q> > state,
		   const GroupXML_t& invParam,
		   QuarkSpinType quarkSpinType,
		   bool obsvP,
		   int& ncg_had) const;

  protected:
    //! Return the fermion create state for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

  private:
      EvenOddPrecDWFermActArray() {} // Partial constructor
      EvenOddPrecDWFermActArray& operator=(const EvenOddPrecDWFermActArray& a) { return *this; // satisfy return type
      } // Hide =

  private:
    Handle< CreateFermState<T,P,Q> >  cfs;
    EvenOddPrecDWFermActArrayParams param;
  };

}

#endif
