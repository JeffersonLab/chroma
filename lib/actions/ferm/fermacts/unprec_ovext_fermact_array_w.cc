// $Id: unprec_ovext_fermact_array_w.cc,v 1.12 2004-12-24 04:23:20 edwards Exp $
/*! \file
 *  \brief Unpreconditioned extended-Overlap (5D) (Naryanan&Neuberger) action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/unprec_ovext_fermact_array_w.h"
#include "actions/ferm/linop/unprec_ovext_linop_array_w.h"
#include "actions/ferm/linop/lmdagm.h"

#include "actions/ferm/invert/invcg2_array.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermbcs/fermbcs_w.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace UnprecOvExtFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct< multi1d<LatticeFermion> >* createFermAct(XMLReader& xml_in,
								const std::string& path)
    {
      return new UnprecOvExtFermActArray(WilsonTypeFermBCArrayEnv::reader(xml_in, path), 
					 UnprecOvExtFermActArrayParams(xml_in, path));
    }

    //! Name to be used
    const std::string name = "UNPRECONDITIONED_OVEXT";

    //! Register the Wilson fermact
    const bool registered = Chroma::TheWilsonTypeFermActArrayFactory::Instance().registerObject(name, createFermAct); 
  }


  //! Read parameters
  UnprecOvExtFermActArrayParams::UnprecOvExtFermActArrayParams(XMLReader& xml, 
							       const std::string& path)
  {
    XMLReader paramtop(xml, path);

    // Read the stuff for the action
    read(paramtop, "OverMass", OverMass);
    read(paramtop, "Mass", Mass);
    read(paramtop, "N5", N5);

    if (paramtop.count("a5") != 0) 
      read(paramtop, "a5", a5);
    else
      a5 = 1.0;
  }


  //! Read parameters
  void read(XMLReader& xml, const string& path, UnprecOvExtFermActArrayParams& param)
  {
    UnprecOvExtFermActArrayParams tmp(xml, path);
    param = tmp;
  }



  //! Produce a linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  const UnprecLinearOperator< multi1d<LatticeFermion>, multi1d<LatticeColorMatrix> >* 
  UnprecOvExtFermActArray::linOp(Handle<const ConnectState> state) const
  {
    return new UnprecOvExtLinOpArray(state->getLinks(),OverMass,Mass,N5);
  }


  //! Produce a M^dag.M linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  const LinearOperator< multi1d<LatticeFermion> >* 
  UnprecOvExtFermActArray::lMdagM(Handle<const ConnectState> state) const
  {
    return new lmdagm<multi1d<LatticeFermion> >(linOp(state));
  }

  //! Propagator of an un-preconditioned Extended-Overlap linear operator
  /*!
   * \param psi      quark propagator ( Modify )
   * \param state    gauge field ( Read )
   * \param chi      source ( Read )
   * \param invParam inverter parameters ( Read (
   * \param ncg_had  number of CG iterations ( Write )
   */
  void 
  UnprecOvExtFermActArray::qprop(LatticeFermion& psi, 
				 Handle<const ConnectState> state, 
				 const LatticeFermion& chi, 
				 const InvertParam_t& invParam,
				 int& ncg_had) const
  {
    START_CODE();

    const int  N5 = size();   // array size better match
    const Real Mass = quark_mass();
    int n_count;
  
    int G5 = Ns*Ns - 1;

    // Initialize the 5D fields
    multi1d<LatticeFermion> chi5(N5);
    multi1d<LatticeFermion> psi5(N5);
    psi5 = zero;
    chi5 = zero;

    psi5[0] = psi;
    chi5[0] = Gamma(G5) * chi;

    // Construct the linear operator
    Handle<const LinearOperator< multi1d<LatticeFermion> > > A(linOp(state));

    switch(invParam.invType)
    {
    case CG_INVERTER: 
      // psi5 = (H_o)^(-2) chi5
      InvCG2(*A, chi5, psi5, invParam.RsdCG, invParam.MaxCG, n_count);

      // chi5 = H_o * (H_o)^(-2) * gamma_5 * chi
      (*A)(chi5, psi5, MINUS);
      break;
  
    case MR_INVERTER:
    case BICG_INVERTER:
      QDP_error_exit("Unsupported inverter type", invParam.invType);
      break;
  
    default:
      QDP_error_exit("Unknown inverter type", invParam.invType);
    }
  
    if ( n_count == invParam.MaxCG )
      QDP_error_exit("no convergence in the inverter", n_count);
  
    ncg_had = n_count;
  
    // Overall normalization
    Real ftmp1 = Real(1) / Real(1 - Mass);

    // Normalize and remove contact term
    psi = ftmp1*(chi5[0] - chi);

    END_CODE();
  }

}
