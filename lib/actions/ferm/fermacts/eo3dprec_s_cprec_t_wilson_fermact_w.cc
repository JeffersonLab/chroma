// $Id: eo3dprec_s_cprec_t_wilson_fermact_w.cc,v 3.2 2007-12-04 16:25:02 bjoo Exp $
/*! \file
 *  \brief EO3DPreconditioned Wilson fermion action
 */

#include "qdp_config.h"
#if QDP_NS == 4
#if QDP_NC == 3
#if QDP_ND == 4

#include "chromabase.h"
#include "actions/ferm/fermacts/eo3dprec_s_cprec_t_wilson_fermact_w.h"
#include "actions/ferm/linop/eo3dprec_s_cprec_t_wilson_linop_w.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"

#include "io/param_io.h"

namespace Chroma
{

  //! Hooks to register the class with the fermact factory
  namespace EO3DPrecSpaceCentralPrecTimeWilsonFermActEnv
  {
    //! Callback function
    WilsonTypeFermAct<LatticeFermion,
		      multi1d<LatticeColorMatrix>,
		      multi1d<LatticeColorMatrix> >* createFermAct4D(XMLReader& xml_in,
								     const std::string& path)
    {
      return new EO3DPrecSpaceCentralPrecTimeWilsonFermAct(CreateFermStateEnv::reader(xml_in, path), 
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
    const std::string name = "EO3DPREC_SPACE_CPREC_TIME_WILSON";

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
  EO3DPrecSpaceCentralPrecTimeLinearOperator<LatticeFermion,
		       multi1d<LatticeColorMatrix>,
		       multi1d<LatticeColorMatrix> >*
  EO3DPrecSpaceCentralPrecTimeWilsonFermAct::linOp(Handle< FermState<T,P,Q> > state) const
  {
    return new EO3DPrecSCprecTWilsonLinOp(state,param.Mass,param.anisoParam); 
  }

}

#endif
#endif
#endif
