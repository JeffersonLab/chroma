// -*- C++ -*-
// $Id: eoprec_nef_fermact_array_w.h,v 3.1 2006-10-19 16:01:27 edwards Exp $
/*! \file
 *  \brief 4D style even-odd preconditioned NEF fermion action
 */

#ifndef __prec_nef_fermact_array_w_h__
#define __prec_nef_fermact_array_w_h__

#include "actions/ferm/fermacts/eoprec_dwf_fermact_base_array_w.h"


namespace Chroma
{
  //! Name and registration
  namespace EvenOddPrecNEFFermActArrayEnv
  {
    extern const std::string name;
    bool registerAll();
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
  /*! \ingroup fermacts
   *
   * 4D style even-odd preconditioned domain-wall fermion action. 
   * Follows notes of Orginos (10/2003)
   *
   * Hopefully, the conventions used here
   * are specified in Phys.Rev.D63:094505,2001 (hep-lat/0005002).
   */
  class EvenOddPrecNEFFermActArray : public EvenOddPrecDWFermActBaseArray<LatticeFermion, 
				   multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General constructor
    EvenOddPrecNEFFermActArray(Handle< CreateFermState<T,P,Q> > cfs_, 
			       const EvenOddPrecNEFFermActArrayParams& p) : 
      cfs(cfs_), params(p)
      {
	QDPIO::cout << "Construct EvenOddPrecNEFFermActArray: OverMass = " << params.OverMass 
		    << "  Mass = " << params.Mass 
		    << "  N5 = " << params.N5 
		    << "  b5 = " << params.b5
		    << "  c5 = " << params.c5
		    << endl;
      }

    //! Copy constructor
    EvenOddPrecNEFFermActArray(const EvenOddPrecNEFFermActArray& a) : 
      cfs(a.cfs), params(a.params) {}

    //! Length of DW flavor index/space
    int size() const {return params.N5;}

    //! Return the quark mass
    Real getQuarkMass() const {return params.Mass;}

    //! Produce an unpreconditioned linear operator for this action with arbitrary quark mass
    UnprecDWLikeLinOpBaseArray<T,P,Q>* unprecLinOp(Handle< FermState<T,P,Q> > state, 
						   const Real& m_q) const;

    //! Produce an even-odd preconditioned linear operator for this action with arbitrary quark mass
    EvenOddPrecDWLikeLinOpBaseArray<T,P,Q>* precLinOp(Handle< FermState<T,P,Q> > state, 
						      const Real& m_q) const;

    //! Destructor is automatic
    ~EvenOddPrecNEFFermActArray() {}

    //! Given a complete propagator as a source, this does all the inversions needed
    /*!
     * This routine is actually generic to Domain Wall fermions (Array) fermions
     *
     * \param q_sol         quark propagator ( Write )
     * \param q_src         source ( Read )
     * \param xml_out       diagnostic output ( Modify )
     * \param state         gauge connection state ( Read )
     * \param t_src         time slice of source ( Read )
     * \param j_decay       direction of decay ( Read )
     * \param invParam      inverter parameters ( Read )
     * \param obsvP         compute currents and residual mass ( Read )
     * \param ncg_had       number of CG iterations ( Write )
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

  protected:
    //! Return the fermion BC object for this action
    const CreateFermState<T,P,Q>& getCreateState() const {return *cfs;}

    //! Partial constructor
    EvenOddPrecNEFFermActArray() {}
    //! Hide =
    void operator=(const EvenOddPrecNEFFermActArray& a) {}

  private:
    Handle< CreateFermState<T,P,Q> >  cfs;
    EvenOddPrecNEFFermActArrayParams params;
  };

}

#endif
