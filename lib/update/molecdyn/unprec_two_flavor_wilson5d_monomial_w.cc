// $Id: unprec_two_flavor_wilson5d_monomial_w.cc,v 1.4 2005-01-10 19:57:02 edwards Exp $
/*! @file
 * @brief Two-flavor collection of unpreconditioned 5D ferm monomials
 */

#include "chromabase.h"
#include "update/molecdyn/unprec_two_flavor_wilson5d_monomial_w.h"
#include "update/molecdyn/monomial_factory.h"

#include "io/param_io.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/invert/invcg2_array.h"

#include "actions/ferm/fermacts/unprec_dwf_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_ovdwf_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_ovext_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_ovlap_contfrac5d_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_nef_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_zolo_nef_fermact_array_w.h"

namespace Chroma 
{ 
 
  namespace UnprecTwoFlavorWilsonTypeFermMonomial5DEnv 
  {
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialDWF(XMLReader& xml, const string& path) 
    {
      return new UnprecTwoFlavorWilsonTypeFermMonomial5D(
	UnprecDWFermActArrayEnv::name,
	UnprecTwoFlavorWilsonTypeFermMonomial5DParams(xml, path));
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialZoloNEF(XMLReader& xml, const string& path) 
    {
      return new UnprecTwoFlavorWilsonTypeFermMonomial5D(
	UnprecZoloNEFFermActArrayEnv::name,
	UnprecTwoFlavorWilsonTypeFermMonomial5DParams(xml, path));
    }
    
    //! Register all the objects
    bool registerAll()
    {
      bool foo = true;
      const std::string prefix = "TWO_FLAVOR_";
      const std::string suffix = "_FERM_MONOMIAL";

      // Use a pattern to register all the qualifying fermacts
      foo &= UnprecDWFermActArrayEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecDWFermActArrayEnv::name+suffix, 
							   createMonomialDWF);

      foo &= UnprecZoloNEFFermActArrayEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecZoloNEFFermActArrayEnv::name+suffix, 
							   createMonomialZoloNEF);
      return foo;
    }

    //! Register the fermact
    const bool registered = registerAll();
  } //end namespace Unprec TwoFlavorWilsonFermMonomialEnv


  // Read the parameters
  UnprecTwoFlavorWilsonTypeFermMonomial5DParams::UnprecTwoFlavorWilsonTypeFermMonomial5DParams(XMLReader& xml_in, const string& path)
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

    QDPIO::cout << "UnprecTwoFlavorWilsonTypeFermMonomial5DParams: read " << ferm_act << endl;
  }

  //! Read Parameters
  void read(XMLReader& xml, const std::string& path,
	    UnprecTwoFlavorWilsonTypeFermMonomial5DParams& params) {
    UnprecTwoFlavorWilsonTypeFermMonomial5DParams tmp(xml, path);
    params = tmp;
  }

  //! Write Parameters
  void write(XMLWriter& xml, const std::string& path,
	     const UnprecTwoFlavorWilsonTypeFermMonomial5DParams& params) {
    // Not implemented
  }

  // Constructor
  UnprecTwoFlavorWilsonTypeFermMonomial5D::UnprecTwoFlavorWilsonTypeFermMonomial5D(
    const string& name_,
    const UnprecTwoFlavorWilsonTypeFermMonomial5DParams& param_) 
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


    QDPIO::cout << "UnprecTwoFlavorWilsonTypeFermMonomial5D: construct " << fermact_string << endl;

    const FermionAction<LatticeFermion>* tmp_act = TheFermionActionFactory::Instance().createObject(fermact_string, fermact_reader, "./FermionAction");
  

    const UnprecWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >* downcast=dynamic_cast<const UnprecWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >*>(tmp_act);

    // Check success of the downcast 
    if( downcast == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to UnprecWilsonTypeFermAct5D in UnprecTwoFlavorWilsonTypeFermMonomial5D()" << endl;
      QDP_abort(1);
    }

    fermact = downcast;    

    QDPIO::cout << "UnprecTwoFlavorWilsonTypeFermMonomial5D: finished " << fermact_string << endl;
  }


  // Do inversion M^dag M X = V^dag phi ?
  void
  UnprecTwoFlavorWilsonTypeFermMonomial5D::getX(
    multi1d<LatticeFermion>& X, 
    const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const
  {
    // Upcast the fermact
    const FermAct5D<LatticeFermion>& FA = getFermAct();

    // Make the state
    Handle< const ConnectState > state(FA.createState(s.getQ()));

    X.resize(FA.size());
    X = zero;
   
    // Get linop
    Handle< const LinearOperator< multi1d<LatticeFermion> > > M(FA.linOp(state));

    Handle< const LinearOperator< multi1d<LatticeFermion> > > PV(FA.linOpPV(state));

    //    // Get PV linop
    // Handle< const LinearOperator< multi1d<LatticeFermion> > > PV(FA.linOpPV
    // Do the inversion...
    multi1d<LatticeFermion> VdagPhi(FA.size());
    (*PV)(VdagPhi, getPhi(), MINUS);

    /*
    int n_count = invert(X, *M, getPhi());
    */
    int n_count = invert(X, *M, VdagPhi);
  }

  
  // Get X = (PV^dag*PV)^{-1} eta
  void
  UnprecTwoFlavorWilsonTypeFermMonomial5D::getXPV(
    multi1d<LatticeFermion>& X, 
    const multi1d<LatticeFermion>& eta, 
    const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const
  {
    // Upcast the fermact
    const FermAct5D<LatticeFermion>& FA = getFermAct();

    // Make the state
    Handle< const ConnectState > state(FA.createState(s.getQ()));

    X.resize(FA.size());
    X=zero;
   
    // Get linop
    Handle< const LinearOperator< multi1d<LatticeFermion> > > M(FA.linOpPV(state));
    

    // Do the inversion...
    int n_count = invert(X, *M, eta);
  }


  // Get X = (A^dag*A)^{-1} eta
  int
  UnprecTwoFlavorWilsonTypeFermMonomial5D::invert(
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
      QDPIO::cout << "2Flav5D::invert, n_count = " << n_count << endl;
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


