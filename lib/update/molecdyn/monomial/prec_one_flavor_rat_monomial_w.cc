// $Id: prec_one_flavor_rat_monomial_w.cc,v 2.1 2005-10-04 19:23:19 bjoo Exp $
/*! @file
 * @brief One-flavor collection of even-odd preconditioned 4D ferm monomials
 */

#include "update/molecdyn/monomial/prec_one_flavor_rat_monomial_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"
#include "update/molecdyn/monomial/genapprox.h"

#include "io/param_io.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/invert/minvcg.h"

#include "actions/ferm/fermacts/prec_wilson_fermact_w.h"
#include "actions/ferm/fermacts/prec_parwilson_fermact_w.h"
#include "actions/ferm/fermacts/prec_stout_fermact_w.h"

namespace Chroma 
{ 
 
  namespace EvenOddPrecOneFlavorWilsonTypeFermRatMonomialEnv 
  {
    //! Does the work
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomial(const string& name, 
							    XMLReader& xml, const string& path)
    {
      QDPIO::cout << "Create Monomial: " << name << endl;

      return new EvenOddPrecOneFlavorWilsonTypeFermRatMonomial(
	name, EvenOddPrecOneFlavorWilsonTypeFermRatMonomialParams(xml, path));
    }
    
    //! Does the work
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomial(const string& name, 
							    XMLReader& xml, const string& path,
							    int expNumPower, int expDenPower) 
    {
      return new EvenOddPrecOneFlavorWilsonTypeFermRatMonomial(
	name, EvenOddPrecOneFlavorWilsonTypeFermRatMonomialParams(xml, path, 
								  expNumPower, expDenPower));
    }
    

    //----------------------------------------------------------------------
    // One flavor
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialWilson1(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecWilsonFermActEnv::name, xml, path, 1, 1);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialParWilson1(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecParWilsonFermActEnv::name, xml, path, 1, 1);
    }

    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialStout1(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecStoutWilsonTypeFermActEnv::name, xml, path, 1, 1);
    }
    
    //----------------------------------------------------------------------
    // Three flavor
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialWilson3(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecWilsonFermActEnv::name, xml, path, 3, 1);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialParWilson3(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecParWilsonFermActEnv::name, xml, path, 3, 1);
    }

    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialStout3(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecStoutWilsonTypeFermActEnv::name, xml, path, 3, 1);
    }
    
    //----------------------------------------------------------------------
    // Generic fractional flavor
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialWilson(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecWilsonFermActEnv::name, xml, path);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialParWilson(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecParWilsonFermActEnv::name, xml, path);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialStout(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecStoutWilsonTypeFermActEnv::name, xml, path);
    }

    //------------------------------------------------------
    //! Register one flavor
    bool registerOne(const string& prefix, const string& suffix)
    {
      bool foo = true;

      // Use a pattern to register all the qualifying fermacts
      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecWilsonFermActEnv::name+suffix, 
							   createMonomialWilson1);

      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecParWilsonFermActEnv::name+suffix, 
							   createMonomialParWilson1);

      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecStoutWilsonTypeFermActEnv::name+suffix, 
							   createMonomialStout1);
      return foo;
    }

    //------------------------------------------------------
    //! Register three flavor objects
    bool registerThree(const string& prefix, const string& suffix)
    {
      bool foo = true;

      // Use a pattern to register all the qualifying fermacts
      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecWilsonFermActEnv::name+suffix, 
							   createMonomialWilson3);

      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecParWilsonFermActEnv::name+suffix, 
							   createMonomialParWilson3);

      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecStoutWilsonTypeFermActEnv::name+suffix, 
							   createMonomialStout3);
      return foo;
    }

    //------------------------------------------------------
    //! Register generic fractional flavor
    bool registerGeneric(const string& prefix, const string& suffix)
    {
      bool foo = true;

      // Use a pattern to register all the qualifying fermacts
      foo &= EvenOddPrecWilsonFermActEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecWilsonFermActEnv::name+suffix, 
							   createMonomialWilson);

      foo &= EvenOddPrecParWilsonFermActEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecParWilsonFermActEnv::name+suffix, 
							   createMonomialParWilson);

      foo &= EvenOddPrecStoutWilsonTypeFermActEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecStoutWilsonTypeFermActEnv::name+suffix, 
							   createMonomialStout);
      return foo;
    }

    //------------------------------------------------------
    //! Register all the objects
    bool registerAll()
    {
      bool foo = true;
      const std::string suffix = "_FERM_RAT_MONOMIAL";

      foo &= registerOne(string("ONE_FLAVOR_"), suffix);
      foo &= registerThree(string("THREE_FLAVOR_"), suffix);
      foo &= registerGeneric(string("FRACTIONAL_FLAVOR_"), suffix);
    }

    //! Register the fermact
    const bool registered = registerAll();
  }; //end namespace EvenOddPrec OneFlavorWilsonFermRatMonomialEnv


  //! Remez input
  void read(XMLReader& xml, const string& path, 
	    EvenOddPrecOneFlavorWilsonTypeFermRatMonomialParams::Remez_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "lowerMin", input.lowerMin);
    read(inputtop, "upperMax", input.upperMax);
    read(inputtop, "forceDegree", input.forceDegree);
    read(inputtop, "actionDegree", input.actionDegree);

    if (inputtop.count("digitPrecision") != 0)
      read(inputtop, "digitPrecision", input.digitPrecision);
    else
      input.digitPrecision = 50;
  }


  // Read the parameters
  EvenOddPrecOneFlavorWilsonTypeFermRatMonomialParams::EvenOddPrecOneFlavorWilsonTypeFermRatMonomialParams(XMLReader& xml_in, const string& path)
  {
    // Get the top of the parameter XML tree
    XMLReader paramtop(xml_in, path);
    
    try {
      // Read the inverter Parameters
      read(paramtop, "InvertParam", inv_param);
      read(paramtop, "Remez", remez);
      read(paramtop, "expNumPower", expNumPower);
      read(paramtop, "expDenPower", expDenPower);
      read(paramtop, "nthRoot", nthRoot);
      XMLReader xml_tmp(paramtop, "./FermionAction");
      std::ostringstream os;
      xml_tmp.print(os);
      ferm_act = os.str();
    }
    catch(const string& s) {
      QDPIO::cerr << "Caught Exception while reading parameters: " << s <<endl;
      QDP_abort(1);
    }

    QDPIO::cout << "EvenOddPrecOneFlavorWilsonTypeFermRatMonomialParams: read \n" << ferm_act << endl;
  }

  // Read the parameters
  EvenOddPrecOneFlavorWilsonTypeFermRatMonomialParams::EvenOddPrecOneFlavorWilsonTypeFermRatMonomialParams(
    XMLReader& xml_in, const string& path, int expNumPower_, int expDenPower_)
  {
    // Get the top of the parameter XML tree
    XMLReader paramtop(xml_in, path);
    
    expNumPower = expNumPower_;
    expDenPower = expDenPower_;

    try {
      // Read the inverter Parameters
      read(paramtop, "InvertParam", inv_param);
      read(paramtop, "Remez", remez);
      read(paramtop, "nthRoot", nthRoot);
      XMLReader xml_tmp(paramtop, "./FermionAction");
      std::ostringstream os;
      xml_tmp.print(os);
      ferm_act = os.str();
    }
    catch(const string& s) {
      QDPIO::cerr << "Caught Exception while reading parameters: " << s <<endl;
      QDP_abort(1);
    }

    QDPIO::cout << "EvenOddPrecOneFlavorWilsonTypeFermRatMonomialParams: read \n" << ferm_act << endl;
  }

  //! Read Parameters
  void read(XMLReader& xml, const std::string& path,
	    EvenOddPrecOneFlavorWilsonTypeFermRatMonomialParams& params) {
    EvenOddPrecOneFlavorWilsonTypeFermRatMonomialParams tmp(xml, path);
    params = tmp;
  }

  //! Write Parameters
  void write(XMLWriter& xml, const std::string& path,
	     const EvenOddPrecOneFlavorWilsonTypeFermRatMonomialParams& params) {
    // Not implemented
  }


  // Constructor
  EvenOddPrecOneFlavorWilsonTypeFermRatMonomial::EvenOddPrecOneFlavorWilsonTypeFermRatMonomial(
    const string& name_,
    const EvenOddPrecOneFlavorWilsonTypeFermRatMonomialParams& param) 
  {
    inv_param = param.inv_param;
    nthRoot   = param.nthRoot;

    cout << "Param.ferm_act is : "<< param.ferm_act << endl;

    std::istringstream is(param.ferm_act);
    XMLReader fermact_reader(is);

    cout << "Fermact reader holds: " << endl;
    fermact_reader.print(cout);
    cout << flush << endl;

    // Get the name of the ferm act
    std::string fermact_string;
    try { 
  
      read(fermact_reader, "/FermionAction/FermAct", fermact_string);
      if ( fermact_string != name_ ) { 
	QDPIO::cerr << "Fermion action is not " << name_
		    << " but is: " << fermact_string << endl;
	QDP_abort(1);
      }
      QDPIO::cout << "Fermact string is " << fermact_string << endl;
    }
    catch( const std::string& e) { 
      QDPIO::cerr << "Error grepping the fermact name: " << e<<  endl;
      QDP_abort(1);
    }


    const FermionAction<LatticeFermion>* tmp_act = TheFermionActionFactory::Instance().createObject(fermact_string, fermact_reader, "/FermionAction");
  

    const EvenOddPrecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >* downcast=dynamic_cast<const EvenOddPrecWilsonTypeFermAct< LatticeFermion, multi1d<LatticeColorMatrix> >*>(tmp_act);

    // Check success of the downcast 
    if( downcast == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to EvenOddPrecWilsonTypeFermAct in EvenOddPrecOneFlavorWilsonTypeFermRatMonomial()" << endl;
      QDP_abort(1);
    }

    fermact = downcast;    

    //*********************************************************************
    // Remez approx
    // M term
    QDPIO::cout << "Normal operator PFE" << endl;
    generateApprox(fpfe, spfe, sipfe,
		   param.remez.lowerMin, param.remez.upperMax, 
		   -param.expNumPower, 2*param.expDenPower*nthRoot, 
		   param.remez.forceDegree, param.remez.actionDegree,
		   param.remez.digitPrecision);
    //*********************************************************************

    QDPIO::cout << "DONECONSTRUCTIN " << endl;
  }


  // Do inversion M^dag M X = phi ?
  int
  EvenOddPrecOneFlavorWilsonTypeFermRatMonomial::getX(
    multi1d<LatticeFermion>& X, 
    const multi1d<Real>& shifts, 
    const LatticeFermion& chi, 
    const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const
  {
    // Upcast the fermact
    const FermAct4D<LatticeFermion>& FA = getFermAct();

    // Make the state
    Handle< const ConnectState > state(FA.createState(s.getQ()));

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
      QDPIO::cout << "1Flav::invert, n_count = " << n_count << endl;
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


