// $Id: unprec_dwftransf_fermact_w.cc,v 1.3 2004-11-16 06:09:09 bjoo Exp $
/*! \file
 *  \brief Unpreconditioned Wilson fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/unprec_dwftransf_fermact_w.h"
#include "actions/ferm/linop/unprec_dwftransf_linop_w.h"
#include "actions/ferm/linop/lmdagm.h"
#include "io/param_io.h"
#include "actions/ferm/fermacts/fermfactory_w.h"

#include <string>

using namespace std;
using namespace QDP;
using namespace Chroma;

namespace Chroma
{

  //! Hooks to register the class with the fermact factory
  namespace UnprecDWFTransfFermActEnv
  {
    //! Callback function
    WilsonTypeFermAct<LatticeFermion>* createFermAct(Handle< FermBC<LatticeFermion> > fbc,
						     XMLReader& xml_in,
						     const std::string& path)
    {
      return new UnprecDWFTransfFermAct(fbc, UnprecDWFTransfFermActParams(xml_in, path));
    }

    //! Name to be used
    const std::string name = "UNPRECONDITIONED_DWFTRANSF";

    //! Register the DWFTransf fermact
    const bool registered = TheWilsonTypeFermActFactory::Instance().registerObject(name, createFermAct);
  }


  //! Read parameters
  UnprecDWFTransfFermActParams::UnprecDWFTransfFermActParams(XMLReader& xml, const string& path)
  {
    try {
      XMLReader paramtop(xml, path);
    
      // Read the stuff for the action
      if (paramtop.count("Mass") != 0) {
	read(paramtop, "Mass", Mass);
	if (paramtop.count("Kappa") != 0) {
	  QDPIO::cerr << "Error: found both a Kappa and a Mass tag" << endl;
	  QDP_abort(1);
	}
      }
      else if (paramtop.count("Kappa") != 0) {
	Real Kappa;
	read(paramtop, "Kappa", Kappa);
	Mass = kappaToMass(Kappa);    // Convert Kappa to Mass
      }
      else {
	QDPIO::cerr << "Error: neither Mass or Kappa found" << endl;
	QDP_abort(1);
      }
      
      // Read b5 c5 and the solver params
      read(paramtop, "b5", b5);
      read(paramtop, "c5", c5);
      // read(paramtop, "invParam", invParam);
      read(paramtop, "RsdCG", invParam.RsdCG);
      read(paramtop, "MaxCG", invParam.MaxCG);
    }
    catch( const string& e) { 
      QDPIO::cerr << "Caught exception reading XML " << e << endl;
      QDP_abort(1);
    }
      
  }


  //! Read parameters
  void read(XMLReader& xml, const string& path, UnprecDWFTransfFermActParams& param)
  {
    UnprecDWFTransfFermActParams tmp(xml, path);
    param = tmp;
  }



  //! Produce a linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  const LinearOperator<LatticeFermion>*
  UnprecDWFTransfFermAct::linOp(Handle<const ConnectState> state) const
  {
    return new UnprecDWFTransfLinOp(state->getLinks(),
				    param.Mass,
				    param.b5,
				    param.c5,
				    param.invParam);
  }


  //! Produce a M^dag.M linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state    gauge field     	       (Read)
   */
  /*
  const LinearOperator<LatticeFermion>*
  UnprecDWFTransfFermAct::lMdagM(Handle<const ConnectState> state) const
  {
    return new UnprecDWFTransfMdagMLinOp(state->getLinks(),
					 param.Mass,
					 param.b5,
					 param.c5,
					 param.invParam);
  }
  */


  void
  UnprecDWFTransfFermAct::dsdu(multi1d<LatticeColorMatrix> & ds_u,
			    Handle<const ConnectState> state,
			    const LatticeFermion& psi) const
  {
    START_CODE();
    QDPIO::cerr << "Force term not implemented for DWFTRansf" << endl;
    QDP_abort(1);
    END_CODE();
  }

}
