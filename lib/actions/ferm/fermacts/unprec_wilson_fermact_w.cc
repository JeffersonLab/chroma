// $Id: unprec_wilson_fermact_w.cc,v 1.26 2005-01-02 05:21:10 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"
#include "actions/ferm/linop/unprec_wilson_linop_w.h"
#include "actions/ferm/linop/lmdagm.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermbcs/fermbcs_w.h"

namespace Chroma
{

  //! Hooks to register the class with the fermact factory
  namespace UnprecWilsonFermActEnv
  {
    //! Callback function
    WilsonTypeFermAct<LatticeFermion,multi1d<LatticeColorMatrix> >* createFermAct4D(XMLReader& xml_in,
										    const std::string& path)
    {
      return new UnprecWilsonFermAct(WilsonTypeFermBCEnv::reader(xml_in, path), 
				     UnprecWilsonFermActParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    FermionAction<LatticeFermion>* createFermAct(XMLReader& xml_in,
						 const std::string& path)
    {
      return createFermAct4D(xml_in, path);
    }

    //! Name to be used
    const std::string name = "UNPRECONDITIONED_WILSON";

    //! Register all the factories
    bool registerAll()
    {
      return Chroma::TheFermionActionFactory::Instance().registerObject(name, createFermAct)
	   & Chroma::TheWilsonTypeFermActFactory::Instance().registerObject(name, createFermAct4D);
    }

    //! Register the fermact
    const bool registered = registerAll();
  }


  //! Read parameters
  UnprecWilsonFermActParams::UnprecWilsonFermActParams(XMLReader& xml, const string& path)
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

    // There is always an aniso Param for wilson, so set it to default
    initHeader(anisoParam);

    //  Read optional aniso
    if (paramtop.count("AnisoParam") != 0) 
      read(paramtop, "AnisoParam", anisoParam);
  }


  //! Read parameters
  void read(XMLReader& xml, const string& path, UnprecWilsonFermActParams& param)
  {
    UnprecWilsonFermActParams tmp(xml, path);
    param = tmp;
  }



  //! Produce a linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  const UnprecLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >*
  UnprecWilsonFermAct::linOp(Handle<const ConnectState> state) const
  {
    if (param.anisoParam.anisoP)
    {
      QDPIO::cerr << "UnprecWilsonFermAct::linOp - no aniso support" << endl;
      QDP_abort(1);
    }

    return new UnprecWilsonLinOp(state->getLinks(), param.Mass); 
  }

  //! Produce a M^dag.M linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state    gauge field     	       (Read)
   */
  const LinearOperator<LatticeFermion>*
  UnprecWilsonFermAct::lMdagM(Handle<const ConnectState> state) const
  {
    return new lmdagm<LatticeFermion>(linOp(state));
  }

}
