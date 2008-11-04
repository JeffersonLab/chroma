// $Id: eoprec_twm_fermact_array_w.cc,v 1.1 2008-11-04 18:42:58 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned Twisted-mass where each flavor is one of two array elements
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/eoprec_twm_fermact_array_w.h"
#include "actions/ferm/linop/eoprec_twm_linop_array_w.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"

namespace Chroma
{
  //! Hooks to register the class with the fermact factory
  namespace EvenOddPrecTwmFermActArrayEnv
  {
    //! Callback function
    WilsonTypeFermAct5D<LatticeFermion,
			multi1d<LatticeColorMatrix>,
			multi1d<LatticeColorMatrix> >* createFermAct5D(XMLReader& xml_in,
								       const std::string& path)
    {
      return new EvenOddPrecTwmFermActArray(CreateFermStateEnv::reader(xml_in, path), 
					    EvenOddPrecTwmFermActArrayParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    FermionAction<LatticeFermion,
		  multi1d<LatticeColorMatrix>,
		  multi1d<LatticeColorMatrix> >* createFermAct(XMLReader& xml_in,
							       const std::string& path)
    {
      return createFermAct5D(xml_in, path);
    }

    //! Name to be used
    const std::string name = "TWM";

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheFermionActionFactory::Instance().registerObject(name, createFermAct);
	success &= Chroma::TheWilsonTypeFermAct5DFactory::Instance().registerObject(name, createFermAct5D);
	registered = true;
      }
      return success;
    }
  }


  //! Read parameters
  EvenOddPrecTwmFermActArrayParams::EvenOddPrecTwmFermActArrayParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "mu_sigma", mu_sigma);
    read(paramtop, "mu_delta", mu_delta);
    read(paramtop, "Mass", Mass);
  }

  //! Read parameters
  void read(XMLReader& xml, const string& path, EvenOddPrecTwmFermActArrayParams& param)
  {
    EvenOddPrecTwmFermActArrayParams tmp(xml, path);
    param = tmp;
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
  EvenOddPrecTwmFermActArray::linOp(Handle< FermState<T,P,Q> > state) const
  {
    return new EvenOddPrecTwmLinOpArray(state,param.Mass,param.mu_sigma,param.mu_delta);
  }


  // Return quark prop solver, solution of unpreconditioned system
  SystemSolverArray<LatticeFermion>* 
  EvenOddPrecTwmFermActArray::qpropT(Handle< FermState<T,P,Q> > state,
				     const GroupXML_t& invParam) const
  {
    QDP_error_exit("EvenOddPrecTwmFermActArray::qpropT not implemented yet for this action yet\n");
    return 0;
  }

}
