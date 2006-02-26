#include "chromabase.h"
#include "actions/ferm/fermacts/unprec_stout_fermact_array_w.h"
#include "actions/ferm/fermacts/stout_fermact_params.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermbcs/fermbcs_reader_w.h"

#include "io/param_io.h"
#include "fermbc.h"

namespace Chroma {

  namespace UnprecStoutWilsonTypeFermAct5DEnv {
       //! Callback function
    WilsonTypeFermAct5D<LatticeFermion,multi1d<LatticeColorMatrix> >* createFermAct5D(XMLReader& xml_in,
										    const std::string& path)
    {
      // Create a Dummy FBC. This breaks the mold a little but
      // is OK because this proxy will always use the BC's in the
      // internal action.
      Handle< FermBC<multi1d<LatticeFermion > > > fbc(new PeriodicFermBC<multi1d< LatticeFermion > >());

      return new UnprecStoutWilsonTypeFermAct5D(fbc, 
					      StoutFermActParams(xml_in, path));
    }

  

    //! Callback function
    /*! Differs in return type */
    FermionAction<LatticeFermion>* createFermAct(XMLReader& xml_in,
						 const std::string& path)
    {
      return createFermAct5D(xml_in, path);
    }
  
    //! Name to be used
    const std::string name = "UNPRECONDITIONED_STOUT_5D";
    
    //! Register all the factories
    bool registerAll()
    {
      return Chroma::TheFermionActionFactory::Instance().registerObject(name, createFermAct)
	& Chroma::TheWilsonTypeFermAct5DFactory::Instance().registerObject(name, createFermAct5D);
    }

    const bool registered = registerAll();

  }

};
