// $Id: unprec_clover_fermact_w.cc,v 1.6 2004-09-08 02:48:25 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Clover fermion action
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_clover_linop_w.h"
#include "actions/ferm/fermacts/unprec_clover_fermact_w.h"
#include "actions/ferm/linop/lmdagm.h"

#include "actions/ferm/fermacts/fermfactory_w.h"

namespace Chroma
{

  //! Hooks to register the class with the fermact factory
  namespace UnprecCloverFermActEnv
  {
    //! Callback function
    WilsonTypeFermAct<LatticeFermion>* createFermAct(Handle< FermBC<LatticeFermion> > fbc,
						     XMLReader& xml_in,
						     const std::string& path)
    {
      return new UnprecCloverFermAct(fbc, UnprecCloverFermActParams(xml_in, path));
    }

    //! Name to be used
    const std::string name = "UNPRECONDITIONED_CLOVER";

    //! Register the Clover fermact
    const bool registered = TheWilsonTypeFermActFactory::Instance().registerObject(name, createFermAct);
  }


  //! Read parameters
  UnprecCloverFermActParams::UnprecCloverFermActParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    // Read the stuff for the action
    if (paramtop.count("Mass") != 0) 
    {
      read(paramtop, "Mass", Mass);
      if (paramtop.count("Kappa") != 0) 
      {
	QDPIO::cerr << "Error: found both a Kappa and a Mass tag" << endl;
	QDP_abort(1);
      }
    }
    else if (paramtop.count("Kappa") != 0)
    {
      Real Kappa;
      read(paramtop, "Kappa", Kappa);
      Mass = kappaToMass(Kappa);    // Convert Kappa to Mass
    }
    else
    {
      QDPIO::cerr << "Error: neither Mass or Kappa found" << endl;
      QDP_abort(1);
    }

    read(paramtop, "ClovCoeff", ClovCoeff);
    read(paramtop, "u0", u0);
  }

  //! Read parameters
  void read(XMLReader& xml, const string& path, UnprecCloverFermActParams& param)
  {
    UnprecCloverFermActParams tmp(xml, path);
    param = tmp;
  }



  //! Produce a linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  const LinearOperator<LatticeFermion>* 
  UnprecCloverFermAct::linOp(Handle<const ConnectState> state) const
  {
    return new UnprecCloverLinOp(state->getLinks(),Mass,ClovCoeff,u0);
  }

  //! Produce a M^dag.M linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  const LinearOperator<LatticeFermion>* 
  UnprecCloverFermAct::lMdagM(Handle<const ConnectState> state) const
  {
    return new lmdagm<LatticeFermion>(linOp(state));
  }

}

