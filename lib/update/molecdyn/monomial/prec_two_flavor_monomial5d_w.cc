// $Id: prec_two_flavor_monomial5d_w.cc,v 1.1 2005-01-28 02:15:32 edwards Exp $
/*! @file
 * @brief Two-flavor collection of even-odd preconditioned 5D ferm monomials
 */

#include "update/molecdyn/monomial/prec_two_flavor_monomial5d_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"

#include "io/param_io.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/invert/invcg2_array.h"

#include "actions/ferm/fermacts/prec_dwf_fermact_array_w.h"
#include "actions/ferm/fermacts/prec_ovdwf_fermact_array_w.h"
#include "actions/ferm/fermacts/prec_nef_fermact_array_w.h"
#include "actions/ferm/fermacts/prec_zolo_nef_fermact_array_w.h"
#include "actions/ferm/fermacts/prec_ovlap_contfrac5d_fermact_array_w.h"

#include "update/molecdyn/predictor/chrono_predictor.h"
#include "update/molecdyn/predictor/chrono_predictor_factory.h"

#include "update/molecdyn/predictor/zero_guess_predictor.h"


namespace Chroma 
{ 
 
  namespace EvenOddPrecTwoFlavorWilsonTypeFermMonomial5DEnv 
  {
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialDWF(XMLReader& xml, const string& path) 
    {
      return new EvenOddPrecTwoFlavorWilsonTypeFermMonomial5D(
	EvenOddPrecDWFermActArrayEnv::name,
	EvenOddPrecTwoFlavorWilsonTypeFermMonomial5DParams(xml, path));
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialOvDWF(XMLReader& xml, const string& path) 
    {
      return new EvenOddPrecTwoFlavorWilsonTypeFermMonomial5D(
	EvenOddPrecOvDWFermActArrayEnv::name,
	EvenOddPrecTwoFlavorWilsonTypeFermMonomial5DParams(xml, path));
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialNEF(XMLReader& xml, const string& path) 
    {
      return new EvenOddPrecTwoFlavorWilsonTypeFermMonomial5D(
	EvenOddPrecNEFFermActArrayEnv::name,
	EvenOddPrecTwoFlavorWilsonTypeFermMonomial5DParams(xml, path));
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialZoloNEF(XMLReader& xml, const string& path) 
    {
      return new EvenOddPrecTwoFlavorWilsonTypeFermMonomial5D(
	EvenOddPrecZoloNEFFermActArrayEnv::name,
	EvenOddPrecTwoFlavorWilsonTypeFermMonomial5DParams(xml, path));
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialContFrac(XMLReader& xml, const string& path) 
    {
      return new EvenOddPrecTwoFlavorWilsonTypeFermMonomial5D(
	EvenOddPrecOvlapContFrac5DFermActArrayEnv::name,
	EvenOddPrecTwoFlavorWilsonTypeFermMonomial5DParams(xml, path));
    }
    
    //! Register all the objects
    bool registerAll()
    {
      bool foo = true;
      const std::string prefix = "TWO_FLAVOR_";
      const std::string suffix = "_FERM_MONOMIAL";

      // Use a pattern to register all the qualifying fermacts
      foo &= EvenOddPrecDWFermActArrayEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecDWFermActArrayEnv::name+suffix, 
							   createMonomialDWF);

      foo &= EvenOddPrecOvDWFermActArrayEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecOvDWFermActArrayEnv::name+suffix, 
							   createMonomialOvDWF);

      foo &= EvenOddPrecNEFFermActArrayEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecNEFFermActArrayEnv::name+suffix, 
							   createMonomialNEF);

      foo &= EvenOddPrecZoloNEFFermActArrayEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecZoloNEFFermActArrayEnv::name+suffix, 
							   createMonomialZoloNEF);

      foo &= EvenOddPrecOvlapContFrac5DFermActArrayEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecOvlapContFrac5DFermActArrayEnv::name+suffix, 
							   createMonomialContFrac);

      return foo;
    }

    //! Register the fermact
    const bool registered = registerAll();
  }; //end namespace EvenOddPrec TwoFlavorWilsonFermMonomialEnv


  // Read the parameters
  EvenOddPrecTwoFlavorWilsonTypeFermMonomial5DParams::EvenOddPrecTwoFlavorWilsonTypeFermMonomial5DParams(XMLReader& xml_in, const string& path)
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

    QDPIO::cout << "EvenOddPrecTwoFlavorWilsonTypeFermMonomial5DParams: read " << ferm_act << endl;
  }

  //! Read Parameters
  void read(XMLReader& xml, const std::string& path,
	    EvenOddPrecTwoFlavorWilsonTypeFermMonomial5DParams& params) {
    EvenOddPrecTwoFlavorWilsonTypeFermMonomial5DParams tmp(xml, path);
    params = tmp;
  }

  //! Write Parameters
  void write(XMLWriter& xml, const std::string& path,
	     const EvenOddPrecTwoFlavorWilsonTypeFermMonomial5DParams& params) {
    // Not implemented
  }


  // Constructor
  EvenOddPrecTwoFlavorWilsonTypeFermMonomial5D::EvenOddPrecTwoFlavorWilsonTypeFermMonomial5D(
    const string& name_,
    const EvenOddPrecTwoFlavorWilsonTypeFermMonomial5DParams& param_) 
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
  

    const EvenOddPrecWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >* downcast=dynamic_cast<const EvenOddPrecWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >*>(tmp_act);

    // Check success of the downcast 
    if( downcast == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to EvenOddPrecWilsonTypeFermAct5D in EvenOddPrecTwoFlavorWilsonTypeFermMonomial5D()" << endl;
      QDP_abort(1);
    }

    fermact = downcast;    

    // Get Chronological predictor
    AbsChronologicalPredictor5D<LatticeFermion>* tmp=0x0;
    if( param_.predictor_xml == "" ) {
      // No predictor specified use zero guess
      tmp = new ZeroGuess5DChronoPredictor(fermact->size());
    }
    else {
      try { 
	std::string chrono_name;
	std::istringstream chrono_is(param_.predictor_xml);
	XMLReader chrono_xml(chrono_is);
	read(chrono_xml, "/ChronologicalPredictor/Name", chrono_name);
	tmp = The5DChronologicalPredictorFactory::Instance().createObject(chrono_name, fermact->size(), chrono_xml, "/ChronologicalPredictor");
      }
      catch(const std::string& e ) { 
	QDPIO::cerr << "Caught Exception Reading XML: " << e << endl;
	QDP_abort(1);
      }
    }
    
    if( tmp == 0x0 ) { 
      QDPIO::cerr << "Failed to create the 5D ChronoPredictor" << endl;
      QDP_abort(1);
    }
    chrono_predictor = tmp;
  }

  // Do inversion M^dag M X = phi
  // X allocated and passed in 
  int
  EvenOddPrecTwoFlavorWilsonTypeFermMonomial5D::getX(
    multi1d<LatticeFermion>& X, 
    const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const
  {
    // Upcast the fermact
    const FermAct5D<LatticeFermion>& FA = getFermAct();

    // Make the state
    Handle< const ConnectState > state(FA.createState(s.getQ()));

    // Get linop
    Handle< const LinearOperator< multi1d<LatticeFermion> > > M(FA.linOp(state));
    // Get PV
    Handle< const LinearOperator< multi1d<LatticeFermion> > > PV(FA.linOpPV(state));

    multi1d<LatticeFermion> VdagPhi(FA.size());
    
    (*PV)(VdagPhi, getPhi(), MINUS);

    // Do the inversion...
    int n_count = invert(X, *M, VdagPhi);
    return n_count;
  }

  
  // Get X = (PV^dag*PV)^{-1} eta
  // X allocated and passed in
  int 
  EvenOddPrecTwoFlavorWilsonTypeFermMonomial5D::getXPV(
    multi1d<LatticeFermion>& X, 
    const multi1d<LatticeFermion>& eta, 
    const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const
  {
    // Upcast the fermact
    const FermAct5D<LatticeFermion>& FA = getFermAct();

    // Make the state
    Handle< const ConnectState > state(FA.createState(s.getQ()));
   
    // Get linop
    Handle< const LinearOperator< multi1d<LatticeFermion> > > M(FA.linOpPV(state));

    // Do the inversion...
    int n_count = invert(X, *M, eta);
    return n_count;
  }


  // Get X = (A^dag*A)^{-1} eta
  int
  EvenOddPrecTwoFlavorWilsonTypeFermMonomial5D::invert(
    multi1d<LatticeFermion>& X, 
    const LinearOperator< multi1d<LatticeFermion> >& M,
    const multi1d<LatticeFermion>& eta) const
  {
    int n_count =0;

    // Do the inversion...
    switch( inv_param.invType) {
    case CG_INVERTER:
    {
      // Solve MdagM X = eta
      InvCG2(M, eta, X, inv_param.RsdCG, inv_param.MaxCG, n_count);
      QDPIO::cout << "2Flav5D::invert,  n_count = " << n_count << endl;
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


