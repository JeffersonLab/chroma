// $Id: prec_two_flavor_wilson_monomial_w.cc,v 1.2 2005-01-13 16:10:30 bjoo Exp $
/*! @file
 * @brief Two-flavor collection of even-odd preconditioned 4D ferm monomials
 */

#include "chromabase.h"
#include "update/molecdyn/monomial/prec_two_flavor_wilson_monomial_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"

#include "io/param_io.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/invert/invcg2.h"

#include "actions/ferm/fermacts/prec_wilson_fermact_w.h"
#include "actions/ferm/fermacts/prec_parwilson_fermact_w.h"

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
	EvenOddPrecTwoFlavorWilsonTypeFermMonomialParams(xml, path));
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialParWilson(XMLReader& xml, const string& path) 
    {
      QDPIO::cout << "Create Monomial: " << EvenOddPrecParWilsonFermActEnv::name << endl;

      return new EvenOddPrecTwoFlavorWilsonTypeFermMonomial(
	EvenOddPrecParWilsonFermActEnv::name,
	EvenOddPrecTwoFlavorWilsonTypeFermMonomialParams(xml, path));
    }
    
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
      return foo;
    }

    //! Register the fermact
    const bool registered = registerAll();
  }; //end namespace EvenOddPrec TwoFlavorWilsonFermMonomialEnv


  // Read the parameters
  EvenOddPrecTwoFlavorWilsonTypeFermMonomialParams::EvenOddPrecTwoFlavorWilsonTypeFermMonomialParams(XMLReader& xml_in, const string& path)
  {
    // Get the top of the parameter XML tree
    XMLReader paramtop(xml_in, path);
    
    try {
      // Read the inverter Parameters
      read(paramtop, "./InvertParam", inv_param);
      XMLReader xml_tmp(paramtop, "./FermionAction");
      std::ostringstream os;
      xml_tmp.print(os);
      ferm_act = os.str();

      if( paramtop.count("./ChronologicalPredictor") == 0 ) {
	predictor_xml="";
      }
      else {
	XMLReader chrono_xml_reader(paramtop, "./ChronologicalPredictor");
	std::ostringstream chrono_os;
	chrono_xml_reader.print(chrono_os);
	predictor_xml = chrono_os.str();
      }
      
    }
    catch(const string& s) {
      QDPIO::cerr << "Caught Exception while reading parameters: " << s <<endl;
      QDP_abort(1);
    }

    QDPIO::cout << "EvenOddPrecTwoFlavorWilsonTypeFermMonomialParams: read \n" << ferm_act << endl;
  }

  //! Read Parameters
  void read(XMLReader& xml, const std::string& path,
	    EvenOddPrecTwoFlavorWilsonTypeFermMonomialParams& params) {
    EvenOddPrecTwoFlavorWilsonTypeFermMonomialParams tmp(xml, path);
    params = tmp;
  }

  //! Write Parameters
  void write(XMLWriter& xml, const std::string& path,
	     const EvenOddPrecTwoFlavorWilsonTypeFermMonomialParams& params) {
    // Not implemented
  }


  // Constructor
  EvenOddPrecTwoFlavorWilsonTypeFermMonomial::EvenOddPrecTwoFlavorWilsonTypeFermMonomial(
    const string& name_,
    const EvenOddPrecTwoFlavorWilsonTypeFermMonomialParams& param_) 
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

    const FermionAction<LatticeFermion>* tmp_act = TheFermionActionFactory::Instance().createObject(fermact_string, fermact_reader, "./FermionAction");
  

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


  // Do inversion M^dag M X = phi
  int
  EvenOddPrecTwoFlavorWilsonTypeFermMonomial::getX(
    LatticeFermion& X, 
    const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const
  {
    // Upcast the fermact
    const FermAct4D<LatticeFermion>& FA = getFermAct();

    // Make the state
    Handle< const ConnectState > state(FA.createState(s.getQ()));

    // Guess is passed in
   
    // Get linop
    Handle< const LinearOperator<LatticeFermion> > M(FA.linOp(state));

    // Do the inversion...
    int n_count = invert(X, *M, getPhi());
    return n_count;
  }



  // Get X = (A^dag*A)^{-1} eta
  int
  EvenOddPrecTwoFlavorWilsonTypeFermMonomial::invert(
    LatticeFermion& X, 
    const LinearOperator<LatticeFermion>& M,
    const LatticeFermion& eta) const
  {
    int n_count =0;

    // Do the inversion...
    switch( inv_param.invType) {
    case CG_INVERTER:
    {
      // Solve MdagM X = eta
      InvCG2(M, eta, X, inv_param.RsdCG, inv_param.MaxCG, n_count);
      QDPIO::cout << "2Flav::invert,  n_count = " << n_count << endl;
    }
    break;
    default:
    {
      QDPIO::cerr << "Currently only CG Inverter is implemented" << endl;
      QDP_abort(1);
    }
    break;
    };

    return n_count;
  }

  
}; //end namespace Chroma


