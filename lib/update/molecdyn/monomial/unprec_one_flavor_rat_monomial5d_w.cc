// $Id: unprec_one_flavor_rat_monomial5d_w.cc,v 2.1 2005-10-04 19:23:19 bjoo Exp $
/*! @file
 * @brief One-flavor collection of unpreconditioned 5D ferm monomials
 */

#include "update/molecdyn/monomial/unprec_one_flavor_rat_monomial5d_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"
#include "update/molecdyn/monomial/genapprox.h"

#include "io/param_io.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/invert/minvcg_array.h"
#include "actions/ferm/linop/lmdagm.h"

#include "actions/ferm/fermacts/unprec_dwf_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_ovdwf_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_ovext_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_ovlap_contfrac5d_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_nef_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_zolo_nef_fermact_array_w.h"
#include "actions/ferm/fermacts/unprec_stout_fermact_array_w.h"

namespace Chroma 
{ 
 
  namespace UnprecOneFlavorWilsonTypeFermRatMonomial5DEnv 
  {
    //! Does the work
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomial(const string& name, 
							    XMLReader& xml, const string& path)
    {
      QDPIO::cout << "Create Monomial: " << name << endl;
      return new UnprecOneFlavorWilsonTypeFermRatMonomial5D(
	name, UnprecOneFlavorWilsonTypeFermRatMonomial5DParams(xml, path));
    }
    
    //! Does the work
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomial(const string& name, 
							    XMLReader& xml, const string& path,
							    int expNumPower, int expDenPower) 
    {
      QDPIO::cout << "Create Monomial: " << name << endl;
      return new UnprecOneFlavorWilsonTypeFermRatMonomial5D(
	name, UnprecOneFlavorWilsonTypeFermRatMonomial5DParams(xml, path, 
								    expNumPower, expDenPower));
    }
    

    //----------------------------------------------------------------------
    // One flavor
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialDWF1(XMLReader& xml, const string& path) 
    {
      return createMonomial(UnprecDWFermActArrayEnv::name, xml, path, 1, 1);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialOvDWF1(XMLReader& xml, const string& path) 
    {
      return createMonomial(UnprecOvDWFermActArrayEnv::name, xml, path, 1, 1);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialNEF1(XMLReader& xml, const string& path) 
    {
      return createMonomial(UnprecNEFFermActArrayEnv::name, xml, path, 1, 1);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialZoloNEF1(XMLReader& xml, const string& path) 
    {
      return createMonomial(UnprecZoloNEFFermActArrayEnv::name, xml, path, 1, 1);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialContFrac1(XMLReader& xml, const string& path) 
    {
      return createMonomial(UnprecOvlapContFrac5DFermActArrayEnv::name, xml, path, 1, 1);
    }

    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialOvExt1(XMLReader& xml, const string& path) 
    {
      return createMonomial(UnprecOvExtFermActArrayEnv::name, xml, path, 1, 1);
    }

    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialStout1(XMLReader& xml, const string& path) 
    {
      return createMonomial(UnprecStoutWilsonTypeFermAct5DEnv::name, xml, path, 1, 1);
    }


    //----------------------------------------------------------------------
    // Three flavor
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialDWF3(XMLReader& xml, const string& path) 
    {
      return createMonomial(UnprecDWFermActArrayEnv::name, xml, path, 3, 1);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialOvDWF3(XMLReader& xml, const string& path) 
    {
      return createMonomial(UnprecOvDWFermActArrayEnv::name, xml, path, 3, 1);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialNEF3(XMLReader& xml, const string& path) 
    {
      return createMonomial(UnprecNEFFermActArrayEnv::name, xml, path, 3, 1);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialZoloNEF3(XMLReader& xml, const string& path) 
    {
      return createMonomial(UnprecZoloNEFFermActArrayEnv::name, xml, path, 3, 1);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialContFrac3(XMLReader& xml, const string& path) 
    {
      return createMonomial(UnprecOvlapContFrac5DFermActArrayEnv::name, xml, path, 3, 1);
    }

    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialOvExt3(XMLReader& xml, const string& path) 
    {
      return createMonomial(UnprecOvExtFermActArrayEnv::name, xml, path, 3, 1);
    }


    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialStout3(XMLReader& xml, const string& path) 
    {
      return createMonomial(UnprecStoutWilsonTypeFermAct5DEnv::name, xml, path, 3, 1);
    }

    //----------------------------------------------------------------------
    // Generic fractional flavor
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialDWF(XMLReader& xml, const string& path) 
    {
      return createMonomial(UnprecDWFermActArrayEnv::name, xml, path);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialOvDWF(XMLReader& xml, const string& path) 
    {
      return createMonomial(UnprecOvDWFermActArrayEnv::name, xml, path);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialNEF(XMLReader& xml, const string& path) 
    {
      return createMonomial(UnprecNEFFermActArrayEnv::name, xml, path);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialZoloNEF(XMLReader& xml, const string& path) 
    {
      return createMonomial(UnprecZoloNEFFermActArrayEnv::name, xml, path);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialContFrac(XMLReader& xml, const string& path) 
    {
      return createMonomial(UnprecOvlapContFrac5DFermActArrayEnv::name, xml, path);
    }

    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialOvExt(XMLReader& xml, const string& path) 
    {
      return createMonomial(UnprecOvExtFermActArrayEnv::name, xml, path);
    }

    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialStout(XMLReader& xml, const string& path) 
    {
      return createMonomial(UnprecStoutWilsonTypeFermAct5DEnv::name, xml, path);
    }

    //------------------------------------------------------
    //! Register one flavor
    bool registerOne(const string& prefix, const string& suffix)
    {
      bool foo = true;

      // Use a pattern to register all the qualifying fermacts
      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecDWFermActArrayEnv::name+suffix, 
							   createMonomialDWF1);

      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecOvDWFermActArrayEnv::name+suffix, 
							   createMonomialOvDWF1);

      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecNEFFermActArrayEnv::name+suffix, 
							   createMonomialNEF1);

      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecZoloNEFFermActArrayEnv::name+suffix, 
							   createMonomialZoloNEF1);

      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecOvlapContFrac5DFermActArrayEnv::name+suffix, 
							   createMonomialContFrac1);

      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecOvExtFermActArrayEnv::name+suffix, 
							   createMonomialOvExt1);
      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecStoutWilsonTypeFermAct5DEnv::name+suffix, 
							   createMonomialStout1);

      return foo;
    }

    //------------------------------------------------------
    //! Register three flavor objects
    bool registerThree(const string& prefix, const string& suffix)
    {
      bool foo = true;

      // Use a pattern to register all the qualifying fermacts
      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecDWFermActArrayEnv::name+suffix, 
							   createMonomialDWF3);

      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecOvDWFermActArrayEnv::name+suffix, 
							   createMonomialOvDWF3);

      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecNEFFermActArrayEnv::name+suffix, 
							   createMonomialNEF3);

      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecZoloNEFFermActArrayEnv::name+suffix, 
							   createMonomialZoloNEF3);

      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecOvlapContFrac5DFermActArrayEnv::name+suffix, 
							   createMonomialContFrac3);

      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecOvExtFermActArrayEnv::name+suffix, 
							   createMonomialOvExt3);

      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecStoutWilsonTypeFermAct5DEnv::name+suffix, 
							   createMonomialStout3);
      return foo;
    }

    //------------------------------------------------------
    //! Register generic fractional flavor
    bool registerGeneric(const string& prefix, const string& suffix)
    {
      bool foo = true;

      // Use a pattern to register all the qualifying fermacts
      foo &= UnprecDWFermActArrayEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecDWFermActArrayEnv::name+suffix, 
							   createMonomialDWF);

      foo &= UnprecOvDWFermActArrayEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecOvDWFermActArrayEnv::name+suffix, 
							   createMonomialOvDWF);

      foo &= UnprecNEFFermActArrayEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecNEFFermActArrayEnv::name+suffix, 
							   createMonomialNEF);

      foo &= UnprecZoloNEFFermActArrayEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecZoloNEFFermActArrayEnv::name+suffix, 
							   createMonomialZoloNEF);

      foo &= UnprecOvlapContFrac5DFermActArrayEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecOvlapContFrac5DFermActArrayEnv::name+suffix, 
							   createMonomialContFrac);

      foo &= UnprecOvExtFermActArrayEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecOvExtFermActArrayEnv::name+suffix, 
							   createMonomialOvExt);

      foo &= UnprecStoutWilsonTypeFermAct5DEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+UnprecStoutWilsonTypeFermAct5DEnv::name+suffix, 
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
  } //end namespace Unprec OneFlavorWilsonFermRatMonomialEnv


  //! Remez input
  void read(XMLReader& xml, const string& path, 
	    UnprecOneFlavorWilsonTypeFermRatMonomial5DParams::Remez_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "lowerMin", input.lowerMin);
    read(inputtop, "upperMax", input.upperMax);
    read(inputtop, "lowerMinPV", input.lowerMinPV);
    read(inputtop, "upperMaxPV", input.upperMaxPV);
    read(inputtop, "degree", input.degree);
    read(inputtop, "degreePV", input.degreePV);

    if (inputtop.count("digitPrecision") != 0)
      read(inputtop, "digitPrecision", input.digitPrecision);
    else
      input.digitPrecision = 50;
  }


  // Read the parameters
  UnprecOneFlavorWilsonTypeFermRatMonomial5DParams::UnprecOneFlavorWilsonTypeFermRatMonomial5DParams(XMLReader& xml_in, const string& path)
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
      read(paramtop, "nthRootPV", nthRootPV);
      XMLReader xml_tmp(paramtop, "./FermionAction");
      std::ostringstream os;
      xml_tmp.print(os);
      ferm_act = os.str();
    }
    catch(const string& s) {
      QDPIO::cerr << "Caught Exception while reading parameters: " << s <<endl;
      QDP_abort(1);
    }

    QDPIO::cout << "UnprecOneFlavorWilsonTypeFermRatMonomial5DParams: read " << ferm_act << endl;
  }

  // Read the parameters
  UnprecOneFlavorWilsonTypeFermRatMonomial5DParams::UnprecOneFlavorWilsonTypeFermRatMonomial5DParams(
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
      read(paramtop, "nthRootPV", nthRootPV);
      XMLReader xml_tmp(paramtop, "./FermionAction");
      std::ostringstream os;
      xml_tmp.print(os);
      ferm_act = os.str();
    }
    catch(const string& s) {
      QDPIO::cerr << "Caught Exception while reading parameters: " << s <<endl;
      QDP_abort(1);
    }

    QDPIO::cout << "UnprecOneFlavorWilsonTypeFermRatMonomial5DParams: read " << ferm_act << endl;
  }

  //! Read Parameters
  void read(XMLReader& xml, const std::string& path,
	    UnprecOneFlavorWilsonTypeFermRatMonomial5DParams& params) {
    UnprecOneFlavorWilsonTypeFermRatMonomial5DParams tmp(xml, path);
    params = tmp;
  }

  //! Write Parameters
  void write(XMLWriter& xml, const std::string& path,
	     const UnprecOneFlavorWilsonTypeFermRatMonomial5DParams& params) {
    // Not implemented
  }

  // Constructor
  UnprecOneFlavorWilsonTypeFermRatMonomial5D::UnprecOneFlavorWilsonTypeFermRatMonomial5D(
    const string& name_,
    const UnprecOneFlavorWilsonTypeFermRatMonomial5DParams& param) 
  {
    inv_param = param.inv_param;
    nthRoot   = param.nthRoot;
    nthRootPV = param.nthRoot;

    std::istringstream is(param.ferm_act);
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


    QDPIO::cout << "UnprecOneFlavorWilsonTypeFermRatMonomial5D: construct " << fermact_string << endl;

    const FermionAction<LatticeFermion>* tmp_act = TheFermionActionFactory::Instance().createObject(fermact_string, fermact_reader, "/FermionAction");
  

    const UnprecWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >* downcast=dynamic_cast<const UnprecWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >*>(tmp_act);

    // Check success of the downcast 
    if( downcast == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to UnprecWilsonTypeFermAct5D in UnprecOneFlavorWilsonTypeFermRatMonomial5D()" << endl;
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
		   param.remez.degree, param.remez.degree,
		   param.remez.digitPrecision);

    // PV term
    QDPIO::cout << "PV operator PFE" << endl;
    generateApprox(fpvpfe, spvpfe, sipvpfe,
		   param.remez.lowerMinPV, param.remez.upperMaxPV, 
		   param.expNumPower, 2*param.expDenPower*nthRootPV, 
		   param.remez.degreePV, param.remez.degreePV,
		   param.remez.digitPrecision);
    //*********************************************************************

    QDPIO::cout << "UnprecOneFlavorWilsonTypeFermRatMonomial5D: finished " << fermact_string << endl;
  }


  //! Multi-mass solver  (M^dag*M + q_i)^{-1} chi  using partfrac
  int
  UnprecOneFlavorWilsonTypeFermRatMonomial5D::getX(
    multi1d< multi1d<LatticeFermion> >& X, 
    const multi1d<Real>& shifts, 
    const multi1d<LatticeFermion>& chi, 
    const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const
  {
    // Upcast the fermact
    const FermAct5D<LatticeFermion>& FA = getFermAct();

    // Make the state
    Handle< const ConnectState > state(FA.createState(s.getQ()));

    // Get linop
    Handle< const LinearOperator< multi1d<LatticeFermion> > > MdagM(FA.lMdagM(state));

    int n_count = invert(X, shifts, *MdagM, chi);
    return n_count;
  }

  
  //! Multi-mass solver  (V^dag*V + q_i)^{-1} chi  using partfrac
  int
  UnprecOneFlavorWilsonTypeFermRatMonomial5D::getXPV(
    multi1d< multi1d<LatticeFermion> >& X, 
    const multi1d<Real>& shifts, 
    const multi1d<LatticeFermion>& chi, 
    const AbsFieldState<multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >& s) const
  {
    // Upcast the fermact
    const FermAct5D<LatticeFermion>& FA = getFermAct();

    // Make the state
    Handle< const ConnectState > state(FA.createState(s.getQ()));

    // Get linop
    Handle< const LinearOperator< multi1d<LatticeFermion> > > 
      MdagM(new lmdagm< multi1d<LatticeFermion> >(FA.linOpPV(state)));
    
    // Do the inversion...
    int n_count = invert(X, shifts, *MdagM, chi);
    return n_count;
  }


  //! Get X = (A^dag*A + q_i)^{-1} eta
  int
  UnprecOneFlavorWilsonTypeFermRatMonomial5D::invert(
    multi1d< multi1d<LatticeFermion> >& X, 
    const multi1d<Real>& shifts, 
    const LinearOperator< multi1d<LatticeFermion> >& A,
    const multi1d<LatticeFermion>& eta) const
  {
    int n_count = 0;
    multi1d<Real> RsdCG(shifts.size());
    RsdCG = inv_param.RsdCG;

    // X allocated and passed in
//    X=zero;
   
    // Do the inversion...
    switch( inv_param.invType) {
    case CG_INVERTER:
    {
      // Solve A^dag*M X = eta
      MInvCG(A, eta, X, shifts, RsdCG, inv_param.MaxCG, n_count);
      QDPIO::cout << "1Flav5D::invert, n_count = " << n_count << endl;
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

} //end namespace Chroma
