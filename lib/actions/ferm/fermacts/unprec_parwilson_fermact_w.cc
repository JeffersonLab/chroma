// $Id: unprec_parwilson_fermact_w.cc,v 1.3 2004-09-08 02:48:26 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion action with parity breaking term
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/unprec_parwilson_fermact_w.h"
#include "actions/ferm/linop/unprec_parwilson_linop_w.h"
#include "actions/ferm/linop/lmdagm.h"

#include "actions/ferm/fermacts/fermfactory_w.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace UnprecParWilsonFermActEnv
  {
    //! Callback function
    WilsonTypeFermAct<LatticeFermion>* createFermAct(Handle< FermBC<LatticeFermion> > fbc,
						     XMLReader& xml_in,
						     const std::string& path)
    {
      return new UnprecParWilsonFermAct(fbc, UnprecParWilsonFermActParams(xml_in, path));
    }

    //! Name to be used
    const std::string name = "PARWILSON";

    //! Register the ParWilson fermact
    const bool registered = Chroma::TheWilsonTypeFermActFactory::Instance().registerObject(name, createFermAct); 
  }


  //! Read parameters
  UnprecParWilsonFermActParams::UnprecParWilsonFermActParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "Mass", Mass);
    read(paramtop, "H", H);
  }

  //! Read parameters
  void read(XMLReader& xml, const string& path, UnprecParWilsonFermActParams& param)
  {
    UnprecParWilsonFermActParams tmp(xml, path);
    param = tmp;
  }



  //! Produce a linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  const LinearOperator<LatticeFermion>*
  UnprecParWilsonFermAct::linOp(Handle<const ConnectState> state) const
  {
    return new UnprecParWilsonLinOp(state->getLinks(),param.Mass,param.H); 
  }

  //! Produce a M^dag.M linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state    gauge field     	       (Read)
   */
  const LinearOperator<LatticeFermion>*
  UnprecParWilsonFermAct::lMdagM(Handle<const ConnectState> state) const
  {
  return new lmdagm<LatticeFermion>(linOp(state));
  }


  //! Computes the derivative of the fermionic action respect to the link field
  /*!
   *         |  dS      dS_f
   * ds_u -- | ----   + -----   ( Write )
   *         |  dU       dU
   *
   * psi -- [1./(M_dag*M)]*chi_  ( read ) 
   *
   * \param ds_u     result      ( Write )
   * \param state    gauge field ( Read )
   * \param psi      solution to linear system ( Read )
   */
  void
  UnprecParWilsonFermAct::dsdu(multi1d<LatticeColorMatrix> & ds_u,
			       Handle<const ConnectState> state,
			       const LatticeFermion& psi) const
  {
    START_CODE();
  
    QDPIO::cerr << "UnprecParWilsonFermAct::dsdu not implemented" << endl;
    QDP_abort(1);

    END_CODE();
  }

}
