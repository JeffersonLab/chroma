#include "chromabase.h"
#include "update/molecdyn/field_state.h"
#include "update/molecdyn/abs_monomial.h"
#include "update/molecdyn/monomial_factory.h"

#include "io/param_io.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"


#include "linearop.h"
#include "fermact.h"

#include "actions/ferm/invert/invcg2.h"

#include "update/molecdyn/unprec_two_flavor_wilson_ferm_monomial_w.h"

#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"
#include "update/molecdyn/zero_guess_predictor.h"
using namespace QDP;

#include <string>
using namespace std;

namespace Chroma { 
 
  namespace UnprecTwoFlavorWilsonTypeFermMonomialEnv {
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const string& path) {

      
      return new UnprecTwoFlavorWilsonTypeFermMonomial(UnprecTwoFlavorWilsonTypeFermMonomialParams(xml, path));
    }
    

    bool registerAll() 
    {
      bool foo = true;
      const std::string prefix = string("TWO_FLAVOR_");
      const std::string suffix = string("_FERM_MONOMIAL");
    
    
      foo &= UnprecWilsonFermActEnv::registered;

#if 0
      // This causes a segfault? Why?
      foo &= TheMonomialFactory::Instance().registerObject(string(prefix + UnprecWilsonFermActEnv::name + suffix), createMonomial);
#else
      // Ugly hack
     string monomial_name = prefix + "UNPRECONDITIONED_WILSON" + suffix;
      foo &= TheMonomialFactory::Instance().registerObject(monomial_name, createMonomial);
#endif
      return foo;
    }

    const bool registered=registerAll();
    
  }; //end namespace Unprec TwoFlavorWilsonTypeFermMonomialEnv

  // Read the parameters
  UnprecTwoFlavorWilsonTypeFermMonomialParams::UnprecTwoFlavorWilsonTypeFermMonomialParams(XMLReader& xml_in, const string& path)
  {
    // Get the top of the parameter XML tree
    XMLReader paramtop(xml_in, path);
    
    try {
      // Read the inverter Parameters
      read(paramtop, "./InvertParam", inv_param);
  
      // Read the fermact
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

  }

  //! Read Parameters
  void read(XMLReader& xml, const std::string& path,
	    UnprecTwoFlavorWilsonTypeFermMonomialParams& params) {
    UnprecTwoFlavorWilsonTypeFermMonomialParams tmp(xml, path);
    params = tmp;
  }

  //! Write Parameters
  void write(XMLWriter& xml, const std::string& path,
	     const UnprecTwoFlavorWilsonTypeFermMonomialParams& params) {
    // Not implemented
  }

  UnprecTwoFlavorWilsonTypeFermMonomial::UnprecTwoFlavorWilsonTypeFermMonomial(const UnprecTwoFlavorWilsonTypeFermMonomialParams& param_) {

    inv_param = param_.inv_param;

    std::istringstream is(param_.ferm_act);
    XMLReader fermact_reader(is);

    // Get the name of the ferm act
    std::string fermact_string;
    try { 

      read(fermact_reader, "/FermionAction/FermAct", fermact_string);
    }
    catch( const std::string& e) { 
      QDPIO::cerr << "Error grepping the fermact name: " << e<<  endl;
      QDP_abort(1);
    }

    const FermionAction<LatticeFermion>* tmp_act = TheFermionActionFactory::Instance().createObject(fermact_string, fermact_reader, "./FermionAction");
  

    const UnprecWilsonTypeFermAct<LatticeFermion, multi1d<LatticeColorMatrix> >* downcast=dynamic_cast<const UnprecWilsonTypeFermAct<LatticeFermion, multi1d<LatticeColorMatrix> > *>(tmp_act);

    // Check success of the downcast 
    if( downcast == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to UnprecWilsonTypeFermAct in UnprecTwoFlavorWilsonTypeFermMonomial()" << endl;
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
	read(chrono_xml, "./ChronologicalPredictor/Name", chrono_name);
	tmp = The4DChronologicalFactory::Instance().createObject(chrono_name, 
								 chrono_xml, 
								 "./ChronologicalPredictor");
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

  void
  UnprecTwoFlavorWilsonTypeFermMonomial::getX(LatticeFermion& X, 
					  const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const
  {
    // Do the MdagM game:
    const UnprecWilsonTypeFermAct<LatticeFermion, multi1d<LatticeColorMatrix> >& S_w = getFermAct();

    // Make the state
    Handle< const ConnectState > state(S_w.createState(s.getQ()));

   
    // Do the inversion...
    switch( inv_param.invType) {
    case CG_INVERTER:
      {
	// Get linop
	Handle< const LinearOperator<LatticeFermion> > M(S_w.linOp(state));
	int n_count =0;

	// Solve MdagM X = phi
	// Get initial guess from predictor
	*(getChronologicalPredictor())(X);

	// do solve
	InvCG2(*M, getPhi(), X, inv_param.RsdCG, inv_param.MaxCG, n_count);

	// Register nwe vector with predictor
	*(getChronologicalPredictor()).newVector(X);
      }
      break;
    default:
      {
	QDPIO::cerr << "Currently only CG Inverter is implemented" << endl;
	QDP_abort(1);
      }
      break;
    };
  }

}; //end namespace Chroma


