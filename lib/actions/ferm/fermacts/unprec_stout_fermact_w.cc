#include "chromabase.h"
#include "actions/ferm/fermacts/unprec_stout_fermact_w.h"
#include "actions/ferm/fermacts/stout_fermact_params.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermbcs/fermbcs_reader_w.h"

#include "io/param_io.h"
#include "fermbc.h"

namespace Chroma {

  namespace UnprecStoutWilsonTypeFermActEnv {
       //! Callback function
    WilsonTypeFermAct<LatticeFermion,multi1d<LatticeColorMatrix> >* createFermAct4D(XMLReader& xml_in,
										    const std::string& path)
    {
      // Create a Dummy FBC. This breaks the mold a little but
      // is OK because this proxy will always use the BC's in the
      // internal action.
      Handle< FermBC<LatticeFermion> > fbc(new PeriodicFermBC<LatticeFermion>());

      return new UnprecStoutWilsonTypeFermAct(fbc, 
					      StoutFermActParams(xml_in, path));
    }

  

    //! Callback function
    /*! Differs in return type */
    FermionAction<LatticeFermion>* createFermAct(XMLReader& xml_in,
						 const std::string& path)
    {
      return createFermAct4D(xml_in, path);
    }
  
    //! Name to be used
    const std::string name = "UNPRECONDITIONED_STOUT";
    
    //! Register all the factories
    bool registerAll()
    {
      return Chroma::TheFermionActionFactory::Instance().registerObject(name, createFermAct)
	& Chroma::TheWilsonTypeFermActFactory::Instance().registerObject(name, createFermAct4D);
    }

    const bool registered = registerAll();

  }

};
