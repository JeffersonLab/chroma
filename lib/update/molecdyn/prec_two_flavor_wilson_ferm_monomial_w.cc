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

using namespace Chroma;

namespace Chroma { 

 
  namespace EvenOddPrecTwoFlavorWilsonTypeFermMonomialEnv {
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,		   
	      multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const string& path) {
      return new EvenOddPrecTwoFlavorWilsonTypeFermMonomial(EvenOddPrecTwoFlavorWilsonTypeFermMonomialParams(xml, path));
    }

    bool registerAll() 
    {
      bool foo = true;
      std::string prefix = string("TWO_FLAVOR_");
      std::string suffix = string("_FERM_MONOMIAL");
    
    
      foo &= EvenOddPrecWilsonFermActEnv::registered;

#if 0
      // This causes a segfault? Why?
      foo &= TheMonomialFactory::Instance().registerObject(string(prefix + EvenOddPrecWilsonFermActEnv::name + suffix), createMonomial);
#else
      // Ugly hack
     string monomial_name = prefix + "WILSON" + suffix;
      foo &= TheMonomialFactory::Instance().registerObject(monomial_name, createMonomial);
#endif
      return foo;
    }

    const bool registered=registerAll();

    
  }; //end namespace EvenOddPrec TwoFlavorWilsonTypeFermMonomialEnv

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
    }
    catch(const string& s) {
      QDPIO::cerr << "Caught Exception while reading parameters: " << s <<endl;
      QDP_abort(1);
    }

    
#if 0
    // This sanity checking is not needed anymore, since this 
    // monomial may also take other wilsonesque actions
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
#endif
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

  EvenOddPrecTwoFlavorWilsonTypeFermMonomial::EvenOddPrecTwoFlavorWilsonTypeFermMonomial(const EvenOddPrecTwoFlavorWilsonTypeFermMonomialParams& param_) {
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
  

    const EvenOddPrecWilsonTypeFermAct<LatticeFermion, multi1d<LatticeColorMatrix> >* downcast=dynamic_cast<const EvenOddPrecWilsonTypeFermAct<LatticeFermion, multi1d<LatticeColorMatrix> >*>(tmp_act);

    // Check success of the downcast 
    if( downcast == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to EvenOddPrecWilsonTypeFermAct in EvenOddPrecTwoFlavorWilsonTypeFermMonomial()" << endl;
      QDP_abort(1);
    }

    fermact = downcast;    
  }

  void
  EvenOddPrecTwoFlavorWilsonTypeFermMonomial::getX(LatticeFermion& X, 
					  const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const
  {
    // Do the MdagM game:
    const EvenOddPrecWilsonTypeFermAct<LatticeFermion, multi1d<LatticeColorMatrix> >& S_w = getFermAct();

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


