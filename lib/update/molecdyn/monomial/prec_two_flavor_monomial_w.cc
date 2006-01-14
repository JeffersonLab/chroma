// $Id: prec_two_flavor_monomial_w.cc,v 2.4 2006-01-14 05:56:59 edwards Exp $
/*! @file
 * @brief Two-flavor collection of even-odd preconditioned 4D ferm monomials
 */

#include "chromabase.h"
#include "update/molecdyn/monomial/prec_two_flavor_monomial_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"

#include "actions/ferm/fermacts/prec_wilson_fermact_w.h"
#include "actions/ferm/fermacts/prec_parwilson_fermact_w.h"

#if 0
#include "actions/ferm/fermacts/prec_stout_fermact_w.h"
#endif

#include "update/molecdyn/predictor/chrono_predictor.h"
#include "update/molecdyn/predictor/chrono_predictor_factory.h"

#include "update/molecdyn/predictor/zero_guess_predictor.h"


namespace Chroma 
{ 
 
  namespace EvenOddPrecTwoFlavorWilsonTypeFermMonomialEnv 
  {
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialWilson(XMLReader& xml, const string& path) 
    {
      QDPIO::cout << "Create Monomial: " << EvenOddPrecWilsonFermActEnv::name << endl;

      return new EvenOddPrecTwoFlavorWilsonTypeFermMonomial(
	EvenOddPrecWilsonFermActEnv::name,
	TwoFlavorWilsonTypeFermMonomialParams(xml, path));
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialParWilson(XMLReader& xml, const string& path) 
    {
      QDPIO::cout << "Create Monomial: " << EvenOddPrecParWilsonFermActEnv::name << endl;

      return new EvenOddPrecTwoFlavorWilsonTypeFermMonomial(
	EvenOddPrecParWilsonFermActEnv::name,
	TwoFlavorWilsonTypeFermMonomialParams(xml, path));
    }


#if 0 
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialStout(XMLReader& xml, const string& path) 
    {
      QDPIO::cout << "Create Monomial: " << EvenOddPrecStoutWilsonTypeFermActEnv::name << endl;

      return new EvenOddPrecTwoFlavorWilsonTypeFermMonomial(
	EvenOddPrecStoutWilsonTypeFermActEnv::name,
	TwoFlavorWilsonTypeFermMonomialParams(xml, path));
    }
#endif
    
    //! Register all the objects
    bool registerAll()
    {
      bool foo = true;
      const std::string prefix = "TWO_FLAVOR_";
      const std::string suffix = "_FERM_MONOMIAL";

      // Use a pattern to register all the qualifying fermacts
      foo &= EvenOddPrecWilsonFermActEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecWilsonFermActEnv::name+suffix, 
							   createMonomialWilson);

      foo &= EvenOddPrecParWilsonFermActEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecParWilsonFermActEnv::name+suffix, 
							   createMonomialParWilson);
#if 0
      foo &= EvenOddPrecStoutWilsonTypeFermActEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecStoutWilsonTypeFermActEnv::name+suffix,
							   createMonomialStout);
#endif      
      return foo;
    }

    //! Register the fermact
    const bool registered = registerAll();
  }; //end namespace EvenOddPrec TwoFlavorWilsonFermMonomialEnv



  // Constructor
  EvenOddPrecTwoFlavorWilsonTypeFermMonomial::EvenOddPrecTwoFlavorWilsonTypeFermMonomial(
    const string& name_,
    const TwoFlavorWilsonTypeFermMonomialParams& param_) 
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

    QDPIO::cout << "EvanOddPrecTwoFlavorWilsonTypeFermMonomial: construct " << fermact_string << endl;


    const FermionAction<LatticeFermion>* tmp_act = TheFermionActionFactory::Instance().createObject(fermact_string, fermact_reader, "/FermionAction");
  

    const EvenOddPrecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >* downcast=dynamic_cast<const EvenOddPrecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >*>(tmp_act);

    // Check success of the downcast 
    if( downcast == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to EvenOddPrecWilsonTypeFermAct in EvenOddPrecTwoFlavorWilsonTypeFermMonomial()" << endl;
      QDP_abort(1);
    }

    fermact = downcast;    

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
      QDPIO::cerr << "Failed to create ZeroGuess4DChronoPredictor" << endl;
      QDP_abort(1);
    }
    chrono_predictor = tmp;

    
  }
  
}; //end namespace Chroma


