// $Id: unprec_one_flavor_rat_monomial_w.cc,v 1.1 2005-01-28 02:15:33 edwards Exp $
/*! @file
 * @brief One-flavor collection of unpreconditioned 4D ferm monomials
 */

#include "chromabase.h"
#include "update/molecdyn/monomial/unprec_one_flavor_rat_monomial_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"

#include "io/param_io.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/invert/minvcg.h"

#include "actions/ferm/fermacts/unprec_wilson_fermact_w.h"
#include "actions/ferm/fermacts/unprec_parwilson_fermact_w.h"

namespace Chroma 
{ 
 
  namespace UnprecOneFlavorWilsonTypeFermRatMonomialEnv 
  {
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialWilson(XMLReader& xml, const string& path) 
    {
      return new UnprecOneFlavorWilsonTypeFermRatMonomial(
	UnprecWilsonFermActEnv::name,
	UnprecOneFlavorWilsonTypeFermRatMonomialParams(xml, path));
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialParWilson(XMLReader& xml, const string& path) 
    {
      return new UnprecOneFlavorWilsonTypeFermRatMonomial(
	UnprecParWilsonFermActEnv::name,
	UnprecOneFlavorWilsonTypeFermRatMonomialParams(xml, path));
    }
    
    //! Register all the objects
    bool registerAll()
    {
      bool foo = true;
      const std::string prefix = "ONE_FLAVOR_";
      const std::string suffix = "_FERM_RAT_MONOMIAL";

      // Use a pattern to register all the qualifying fermacts
      foo &= UnprecWilsonFermActEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecWilsonFermActEnv::name+suffix, 
							   createMonomialWilson);

      foo &= UnprecParWilsonFermActEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecParWilsonFermActEnv::name+suffix, 
							   createMonomialParWilson);
      return foo;
    }

    //! Register the fermact
    const bool registered = registerAll();
  } //end namespace Unprec OneFlavorWilsonFermRatMonomialEnv


  // Read the parameters
  UnprecOneFlavorWilsonTypeFermRatMonomialParams::UnprecOneFlavorWilsonTypeFermRatMonomialParams(XMLReader& xml_in, const string& path)
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

    QDPIO::cout << "UnprecOneFlavorWilsonTypeFermRatMonomialParams: read " << ferm_act << endl;
  }

  //! Read Parameters
  void read(XMLReader& xml, const std::string& path,
	    UnprecOneFlavorWilsonTypeFermRatMonomialParams& params) {
    UnprecOneFlavorWilsonTypeFermRatMonomialParams tmp(xml, path);
    params = tmp;
  }

  //! Write Parameters
  void write(XMLWriter& xml, const std::string& path,
	     const UnprecOneFlavorWilsonTypeFermRatMonomialParams& params) {
    // Not implemented
  }

  // Constructor
  UnprecOneFlavorWilsonTypeFermRatMonomial::UnprecOneFlavorWilsonTypeFermRatMonomial(
    const string& name_,
    const UnprecOneFlavorWilsonTypeFermRatMonomialParams& param_) 
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


    QDPIO::cout << "UnprecOneFlavorWilsonTypeFermRatMonomial: construct " << fermact_string << endl;

    const FermionAction<LatticeFermion>* tmp_act = TheFermionActionFactory::Instance().createObject(fermact_string, fermact_reader, "./FermionAction");
  

    const UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >* downcast=dynamic_cast<const UnprecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >*>(tmp_act);

    // Check success of the downcast 
    if( downcast == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to UnprecWilsonTypeFermAct in UnprecOneFlavorWilsonTypeFermRatMonomial()" << endl;
      QDP_abort(1);
    }

    fermact = downcast;    

    //*********************************************************************
    // HACK FOR NOW - arbitrarily set the coefficients
    int N = 4;
    FPartFracCoeff.resize(N);
    FPartFracRoot.resize(N);
    HBPartFracCoeff.resize(N);
    HBPartFracRoot.resize(N);

    FPartFracCoeff = 1;
    FPartFracRoot = 1;
    HBPartFracCoeff = 1;
    HBPartFracRoot = 1;
    //*********************************************************************

    QDPIO::cout << "UnprecOneFlavorWilsonTypeFermRatMonomial: finished " << fermact_string << endl;
  }


  // Do inversion M^dag M X = phi ?
  int
  UnprecOneFlavorWilsonTypeFermRatMonomial::getX(
    multi1d<LatticeFermion>& X, 
    const multi1d<Real>& shifts, 
    const LatticeFermion& chi, 
    const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const
  {
    // Upcast the fermact
    const FermAct4D<LatticeFermion>& FA = getFermAct();

    // Make the state
    Handle< const ConnectState > state(FA.createState(s.getQ()));

    // Initial guess for X passed in
    
    // Get linop
    Handle< const LinearOperator<LatticeFermion> > MdagM(FA.lMdagM(state));

    int n_count =0;
    multi1d<Real> RsdCG(shifts.size());
    RsdCG = inv_param.RsdCG;

    // Do the inversion...
    switch( inv_param.invType) {
    case CG_INVERTER:
    {
      // Solve MdagM X = eta
      MInvCG(*MdagM, chi, X, shifts, RsdCG, inv_param.MaxCG, n_count);
      QDPIO::cout << "2Flav::invert, n_count = " << n_count << endl;
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


