#include "chromabase.h"
#include "update/molecdyn/field_state.h"
#include "update/molecdyn/abs_monomial.h"
#include "update/molecdyn/monomial_factory.h"

#include "io/param_io.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/prec_wilson_fermact_w.h"
#include "linearop.h"
#include "actions/ferm/invert/invcg2.h"

#include "update/molecdyn/prec_two_flavor_wilson_ferm_monomial_w.h"
using namespace QDP;

#include <string>
using namespace std;

namespace Chroma { 
 
  namespace EvenOddPrecTwoFlavorWilsonFermMonomialEnv {
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,		   
	      multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const string& path) {
      return new EvenOddPrecTwoFlavorWilsonFermMonomial(EvenOddPrecTwoFlavorWilsonFermMonomialParams(xml, path));
    }
    
    const std::string name="TWO_FLAVOR_WILSON_FERM_MONOMIAL";
    const bool registered=TheMonomialFactory::Instance().registerObject(name, createMonomial);
    
  }; //end namespace EvenOddPrec TwoFlavorWilsonFermMonomialEnv

  // Read the parameters
  EvenOddPrecTwoFlavorWilsonFermMonomialParams::EvenOddPrecTwoFlavorWilsonFermMonomialParams(XMLReader& xml_in, const string& path)
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
    }
    catch(const string& s) {
      QDPIO::cerr << "Caught Exception while reading parameters: " << s <<endl;
      QDP_abort(1);
    }

    std::istringstream is(ferm_act);
    
    XMLReader sanity_check(is);
    try { 
      std::string ferm_name;
      read(sanity_check, "/FermionAction/FermAct", ferm_name);
      if ( ferm_name != EvenOddPrecWilsonFermActEnv::name ) { 
	QDPIO::cerr << "Fermion action is not " << EvenOddPrecWilsonFermActEnv::name
		    << " but is: " << ferm_name << endl;
	QDP_abort(1);
      }
    }
    catch(const std::string& s) { 
      QDPIO::cerr << "Unable to sanity check fermion action" << endl;
      QDP_abort(1);
    }
  }

  //! Read Parameters
  void read(XMLReader& xml, const std::string& path,
	    EvenOddPrecTwoFlavorWilsonFermMonomialParams& params) {
    EvenOddPrecTwoFlavorWilsonFermMonomialParams tmp(xml, path);
    params = tmp;
  }

  //! Write Parameters
  void write(XMLWriter& xml, const std::string& path,
	     const EvenOddPrecTwoFlavorWilsonFermMonomialParams& params) {
    // Not implemented
  }

  EvenOddPrecTwoFlavorWilsonFermMonomial::EvenOddPrecTwoFlavorWilsonFermMonomial(const EvenOddPrecTwoFlavorWilsonFermMonomialParams& param_) {
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
  

    const EvenOddPrecWilsonFermAct* downcast=dynamic_cast<const EvenOddPrecWilsonFermAct*>(tmp_act);

    // Check success of the downcast 
    if( downcast == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to EvenOddPrecWilsonFermAct in EvenOddPrecTwoFlavorWilsonFermMonomial()" << endl;
      QDP_abort(1);
    }

    fermact = downcast;    
  }

  void
  EvenOddPrecTwoFlavorWilsonFermMonomial::getX(LatticeFermion& X, 
					  const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const
  {
    // Do the MdagM game:
    const EvenOddPrecWilsonFermAct& S_w = getFermAct();

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
	X=zero;
	InvCG2(*M, getPhi(), X, inv_param.RsdCG, inv_param.MaxCG, n_count);
	QDPIO::cout << "getX: n_count = " << n_count << endl;
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


