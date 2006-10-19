// -*- C++ -*-
// $Id: eoprec_zolo_nef_fermact_array_w.h,v 3.1 2006-10-19 16:01:27 edwards Exp $
/*! \file
 *  \brief Unpreconditioned NEF domain-wall fermion action
 */

#ifndef __prec_zolo_nef_fermact_array_w_h__
#define __prec_zolo_nef_fermact_array_w_h__

#include "actions/ferm/fermacts/eoprec_dwf_fermact_base_array_w.h"
#include "io/enum_io/enum_coeffs_io.h"


namespace Chroma
{
  //! Name and registration
  namespace EvenOddPrecZoloNEFFermActArrayEnv
  {
    extern const std::string name;
    bool registerAll();
  }
  

  //! Params for NEFF
  struct EvenOddPrecZoloNEFFermActArrayParams
  {
    EvenOddPrecZoloNEFFermActArrayParams() {}
    EvenOddPrecZoloNEFFermActArrayParams(XMLReader& in, const std::string& path);
    
    Real OverMass;  //!< Mass of auxiliary Wilson action
    Real Mass;      //!< Fermion Mass
    Real b5;        //!< b5 in H_T expression
    Real c5;        //!< c5 in H_T expression
    int N5;         //!< Size of 5D extent
    CoeffType approximation_type;  //!< ZOLOTAREV | TANH | Other approximation coeffs
    Real ApproxMin; //!< Approximate min eigenvalue of H_T
    Real ApproxMax; //!< Approximate max eigenvalue of H_T
  };


  // Reader/writers
  void read(XMLReader& xml, const string& path, EvenOddPrecZoloNEFFermActArrayParams& param);
  void write(XMLWriter& xml, const string& path, const EvenOddPrecZoloNEFFermActArrayParams& param);


  //! EvenOddPreconditioned NEF fermion action
  /*! \ingroup fermacts
   *
   * EvenOddPrecondition NEF fermion action. The conventions used here
   * are specified in Phys.Rev.D63:094505,2001 (hep-lat/0005002).
   * See also Brower et.al. LATTICE04
   */
  class EvenOddPrecZoloNEFFermActArray : public EvenOddPrecDWFermActBaseArray<LatticeFermion, 
				   multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! General FermBC
    EvenOddPrecZoloNEFFermActArray(Handle< CreateFermState<T,P,Q> > cfs_, 
				   const EvenOddPrecZoloNEFFermActArrayParams& param_) : 
      cfs(cfs_), params(param_) {}

    //! Copy constructor
    EvenOddPrecZoloNEFFermActArray(const EvenOddPrecZoloNEFFermActArray& a) : 
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
    ~EvenOddPrecZoloNEFFermActArray() {}

    //! Given a complete propagator as a source, this does all the inversions needed
    /*!
     * This routine is actually generic to Domain Wall fermions (Array) fermions
     *
     * \param q_sol        quark propagator ( Write )
     * \param q_src        source ( Read )
     * \param xml_out      diagnostic output ( Modify )
     * \param state        gauge connection state ( Read )
     * \param t_src        time slice of source ( Read )
     * \param j_decay      direction of decay ( Read )
     * \param invParam     inverter parameters ( Read )
     * \param ncg_had      number of CG iterations ( Write )
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
    EvenOddPrecZoloNEFFermActArray() {}
    //! Assignment
    void operator=(const EvenOddPrecZoloNEFFermActArray& a) {}

  private:
    void initCoeffs(multi1d<Real>& b5_arr,
		    multi1d<Real>& c5_arr) const;

  private:
    Handle< CreateFermState<T,P,Q> >  cfs;
    EvenOddPrecZoloNEFFermActArrayParams params;
  };

}

#endif
