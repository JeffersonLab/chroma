// -*- C++ -*-
// $Id: unprec_ovext_fermact_array_w.h,v 1.10 2004-12-12 21:22:15 edwards Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) action
 */

#ifndef __unprec_ovext_fermact_array_w_h__
#define __unprec_ovext_fermact_array_w_h__

#include "fermact.h"

using namespace QDP;

namespace Chroma
{
  //! Name and registration
  namespace UnprecOvExtFermActArrayEnv
  {
    extern const std::string name;
    extern const bool registered;
  }
  

  //! Params for NEFF
  struct UnprecOvExtFermActArrayParams
  {
    UnprecOvExtFermActArrayParams() {}
    UnprecOvExtFermActArrayParams(XMLReader& in, const std::string& path);
    
    Real OverMass;
    Real a5;
    Real Mass;
    int  N5;
  };


  // Reader/writers
  void read(XMLReader& xml, const string& path, UnprecOvExtFermActArrayParams& param);
  void write(XMLReader& xml, const string& path, const UnprecOvExtFermActArrayParams& param);


  //! Unpreconditioned Extended-Overlap (N&N) linear operator
  /*!
   * \ingroup fermact
   *
   * This operator applies the extended version of the hermitian overlap operator
   *   Chi  =   ((1+Mass)/(1-Mass)*gamma_5 + B) . Psi
   *  where  B  is the continued fraction of the pole approx. to eps(H(m))
   */
  class UnprecOvExtFermActArray : public UnprecWilsonTypeFermAct< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> >
  {
  public:
    //! General FermBC
    UnprecOvExtFermActArray(Handle< FermBC< multi1d<LatticeFermion> > > fbc_, 
			    const Real& OverMass_, const Real& Mass_, int N5_) : 
      fbc(fbc_), OverMass(OverMass_), Mass(Mass_), N5(N5_) {a5=1;}

    //! General FermBC
    UnprecOvExtFermActArray(Handle< FermBC< multi1d<LatticeFermion> > > fbc_, 
			    const UnprecOvExtFermActArrayParams& param) :
      fbc(fbc_), OverMass(param.OverMass), Mass(param.Mass), a5(param.a5), N5(param.N5) {}

    //! Copy constructor
    UnprecOvExtFermActArray(const UnprecOvExtFermActArray& a) : 
      fbc(a.fbc), OverMass(a.OverMass), Mass(a.Mass), a5(a.a5), N5(a.N5) {}

    //! Assignment
    UnprecOvExtFermActArray& operator=(const UnprecOvExtFermActArray& a)
      {fbc=a.fbc; OverMass=a.OverMass; Mass=a.Mass; a5=a.a5; N5=a.N5; return *this;}

    //! Return the fermion BC object for this action
    const FermBC< multi1d<LatticeFermion> >& getFermBC() const {return *fbc;}

    //! Length of DW flavor index/space
    int size() const {return N5;}

    //! Return the quark mass
    Real quark_mass() const {return Mass;}

    //! Produce a linear operator for this action
    const UnprecLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> >* linOp(Handle<const ConnectState> state) const;

    //! Produce a linear operator M^dag.M for this action
    const LinearOperator< multi1d<LatticeFermion> >* lMdagM(Handle<const ConnectState> state) const;

    //! Produce a hermitian version of the linear operator
    const LinearOperator< multi1d<LatticeFermion> >* gamma5HermLinOp(Handle<const ConnectState> state) const
      {
	QDPIO::cerr << "UnprecOvExtFermActArray::gamma5HermLinOp not implemented" << endl;
	QDP_abort(1);
	return 0;
      }

    //! Compute quark propagator over base type
    /*! 
     * Solves  M.psi = chi
     *
     * \param psi      quark propagator ( Modify )
     * \param u        gauge field ( Read )
     * \param chi      source ( Modify )
     * \param invParam inverter parameters ( Read (
     * \param ncg_had  number of CG iterations ( Write )
     *
     * NOTE: maybe this should produce a quark prop foundry class object 
     */
    void qprop(LatticeFermion& psi, 
	       Handle<const ConnectState> state, 
	       const LatticeFermion& chi, 
	       const InvertParam_t& invParam,
	       int& ncg_had) const;

    //! Destructor is automatic
    ~UnprecOvExtFermActArray() {}

  private:
    // Hide partial constructor
    UnprecOvExtFermActArray() {}

  private:
    Handle< FermBC< multi1d<LatticeFermion> > >  fbc;
    Real OverMass;
    Real Mass;
    Real a5;
    int  N5;
  };

}

using namespace Chroma;


#endif
