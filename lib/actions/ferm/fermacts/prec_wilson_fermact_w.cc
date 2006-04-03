// $Id: prec_wilson_fermact_w.cc,v 3.0 2006-04-03 04:58:46 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Wilson fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/prec_wilson_fermact_w.h"
#include "actions/ferm/linop/prec_wilson_linop_w.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/ferm_createstate_reader_w.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace EvenOddPrecWilsonFermActEnv
  {
    //! Callback function
    WilsonTypeFermAct<LatticeFermion, 
		      multi1d<LatticeColorMatrix>,
		      multi1d<LatticeColorMatrix> >* createFermAct4D(XMLReader& xml_in,
								     const std::string& path)
    {
      return new EvenOddPrecWilsonFermAct(CreateFermStateEnv::reader(xml_in, path), 
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
    const std::string name("WILSON");

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= Chroma::TheFermionActionFactory::Instance().registerObject(name, createFermAct);
      foo &= Chroma::TheWilsonTypeFermActFactory::Instance().registerObject(name, createFermAct4D);
    }

    //! Register the fermact
    const bool registered = registerAll();
  }


  //! Produce a linear operator for this action
  /*!
   * The operator acts on the odd subset
   *
   * \param state 	    gauge field     	       (Read)
   */
  EvenOddPrecConstDetLinearOperator<LatticeFermion,
				    multi1d<LatticeColorMatrix>,
				    multi1d<LatticeColorMatrix> >* 
  EvenOddPrecWilsonFermAct::linOp(Handle< FermState<T,P,Q> > state) const
  {
    return new EvenOddPrecWilsonLinOp(state,param.Mass,param.anisoParam);
  }

}
