// $Id: unprec_graphene_fermact_w.cc,v 3.1 2007-12-31 23:24:26 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Graphene fermion action.
 *
 * This formulation follows Borici's variant of Creutz's graphene
 * fermion construction. Borici's variant is described in
 * arXiv:0712.4401 and Cruetz's original construction is described
 * in arXiv:0712.1201
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/unprec_graphene_fermact_w.h"
#include "actions/ferm/linop/unprec_graphene_linop_w.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"

#include "io/param_io.h"

namespace Chroma
{

  //! Hooks to register the class with the fermact factory
  namespace UnprecGrapheneFermActEnv
  {
    //! Callback function
    WilsonTypeFermAct<LatticeFermion,
		      multi1d<LatticeColorMatrix>,
		      multi1d<LatticeColorMatrix> >* createFermAct4D(XMLReader& xml_in,
								     const std::string& path)
    {
      return new UnprecGrapheneFermAct(CreateFermStateEnv::reader(xml_in, path), 
				       WilsonFermActParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    FermionAction<LatticeFermion,
		  multi1d<LatticeColorMatrix>,
		  multi1d<LatticeColorMatrix> >* createFermAct(XMLReader& xml_in,
							       const std::string& path)
    {
      return createFermAct4D(xml_in, path);
    }

    //! Name to be used
    const std::string name = "UNPRECONDITIONED_GRAPHENE";

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheFermionActionFactory::Instance().registerObject(name, createFermAct);
	success &= Chroma::TheWilsonTypeFermActFactory::Instance().registerObject(name, createFermAct4D);
	registered = true;
      }
      return success;
    }
  }


  //! Produce a linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  UnprecLinearOperator<LatticeFermion,
		       multi1d<LatticeColorMatrix>,
		       multi1d<LatticeColorMatrix> >*
  UnprecGrapheneFermAct::linOp(Handle< FermState<T,P,Q> > state) const
  {
    return new UnprecGrapheneLinOp(state,param.Mass,param.anisoParam); 
  }

}
