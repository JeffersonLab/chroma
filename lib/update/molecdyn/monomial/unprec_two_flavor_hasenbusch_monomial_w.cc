// $Id: unprec_two_flavor_hasenbusch_monomial_w.cc,v 2.3 2006-01-14 05:22:32 edwards Exp $
/*! @file
 * @brief Two-flavor collection of unpreconditioned 4D ferm monomials
 */

#include "chromabase.h"
#include "update/molecdyn/monomial/unprec_two_flavor_hasenbusch_monomial_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/invert/invcg2.h"

#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"
#include "actions/ferm/fermacts/unprec_parwilson_fermact_w.h"

#if 0
#include "actions/ferm/fermacts/unprec_stout_fermact_w.h"
#endif

#include "update/molecdyn/predictor/chrono_predictor_factory.h"
#include "update/molecdyn/predictor/zero_guess_predictor.h"


namespace Chroma 
{ 
 
  namespace UnprecTwoFlavorHasenbuschWilsonTypeFermMonomialEnv 
  {
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialWilson(XMLReader& xml, const string& path) 
    {
      return new UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial(
	UnprecWilsonFermActEnv::name,
	TwoFlavorHasenbuschWilsonTypeFermMonomialParams(xml, path));
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialParWilson(XMLReader& xml, const string& path) 
    {
      return new UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial(
	UnprecParWilsonFermActEnv::name,
	TwoFlavorHasenbuschWilsonTypeFermMonomialParams(xml, path));
    }

#if 0
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialStout(XMLReader& xml, const string& path) 
    {
      return new UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial(
	UnprecStoutWilsonTypeFermActEnv::name,
	TwoFlavorHasenbuschWilsonTypeFermMonomialParams(xml, path));
    }
#endif   
 
    //! Register all the objects
    bool registerAll()
    {
      bool foo = true;
      const std::string prefix = "TWO_FLAVOR_";
      const std::string suffix = "_HASENBUSCH_FERM_MONOMIAL";

      // Use a pattern to register all the qualifying fermacts
      foo &= UnprecWilsonFermActEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecWilsonFermActEnv::name+suffix, 
							   createMonomialWilson);

      foo &= UnprecParWilsonFermActEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecParWilsonFermActEnv::name+suffix, 
							   createMonomialParWilson);

#if 0
      foo &= UnprecStoutWilsonTypeFermActEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecStoutWilsonTypeFermActEnv::name+suffix, 
							   createMonomialStout);
#endif
      return foo;
    }

    //! Register the fermact
    const bool registered = registerAll();
  } //end namespace Unprec TwoFlavorHasenbuschWilsonFermMonomialEnv


  // Constructor
  UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial::UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial(
    const string& name_,
    const TwoFlavorHasenbuschWilsonTypeFermMonomialParams& param_) 
  {
    inv_param = param_.inv_param;

    std::istringstream is(param_.ferm_act);
    XMLReader fermact_reader(is);

    // Get the name of the ferm act
    std::string fermact_string;
    try { 
      read(fermact_reader, "/FermionAction/FermAct", fermact_string);
      if ( fermact_string != name_ ) { 
	QDPIO::cerr << "Fermion action is not " << name_
		    << " but is: " << fermact_string << endl;
	QDP_abort(1);
      }
    }
    catch( const std::string& e) { 
      QDPIO::cerr << "Error grepping the fermact name: " << e<<  endl;
      QDP_abort(1);
    }


    QDPIO::cout << "UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial: construct " << fermact_string << endl;

   
    const FermionAction<LatticeFermion>* tmp_act = TheFermionActionFactory::Instance().createObject(fermact_string, fermact_reader, "/FermionAction");

    const UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >* downcast=dynamic_cast<const UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >*>(tmp_act);


    // Check success of the downcast 
    if( downcast == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to UnprecWilsonTypeFermAct in UnprecTwoFlavorWilsonTypeFermMonomial()" << endl;
      QDP_abort(1);
    }

    fermact = downcast;    

    // -------------------
    std::istringstream is_prec(param_.ferm_act_prec);
    XMLReader fermact_prec_reader(is_prec);

    // Get the name of the ferm act
    std::string fermact_prec_string;
    try { 
      read(fermact_prec_reader, "/FermionActionPrec/FermAct", fermact_prec_string);
      if ( fermact_prec_string != name_ ) { 
	QDPIO::cerr << "Fermion action is not " << name_
		    << " but is: " << fermact_prec_string << endl;
	QDP_abort(1);
      }
    }
    catch( const std::string& e) { 
      QDPIO::cerr << "Error grepping the fermact name: " << e<<  endl;
      QDP_abort(1);
    }


    QDPIO::cout << "UnprecTwoFlavorHasenbuschWilsonTypeFermMonomial: construct " << fermact_prec_string << endl;

   
    const FermionAction<LatticeFermion>* tmp_act_prec = TheFermionActionFactory::Instance().createObject(fermact_prec_string, fermact_prec_reader, "/FermionActionPrec");

    const UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >* downcast_prec=dynamic_cast<const UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >*>(tmp_act_prec);


    // Check success of the downcast 
    if( downcast_prec == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to UnprecWilsonTypeFermAct in UnprecTwoFlavorWilsonTypeFermMonomial()" << endl;
      QDP_abort(1);
    }

    fermact_prec = downcast_prec;    

    //------------------------------------

    // Get Chronological predictor
    AbsChronologicalPredictor4D<LatticeFermion>* tmp=0x0;
    if( param_.predictor_xml == "" ) {
      // No predictor specified use zero guess
       tmp = new ZeroGuess4DChronoPredictor();
    }
    else {

      
      try { 
	std::string chrono_name;
	std::istringstream chrono_is(param_.predictor_xml);
	XMLReader chrono_xml(chrono_is);
	read(chrono_xml, "/ChronologicalPredictor/Name", chrono_name);
	tmp = The4DChronologicalPredictorFactory::Instance().createObject(chrono_name, 
								 chrono_xml, 
								 "/ChronologicalPredictor");
      }
      catch(const std::string& e ) { 
	QDPIO::cerr << "Caught Exception Reading XML: " << e << endl;
	QDP_abort(1);
      }


    }
     
    if( tmp == 0x0 ) { 
      QDPIO::cerr << "Failed to create the 4D ChronoPredictor" << endl;
      QDP_abort(1);
    }
    chrono_predictor = tmp;

    QDPIO::cout << "UnprecTwoFlavorWilsonTypeFermMonomial: finished " << fermact_string << endl;
  }
}; //end namespace Chroma


