// $Id: prec_constdet_one_flavor_rat_monomial5d_w.cc,v 2.1 2006-01-16 00:33:52 bjoo Exp $
/*! @file
 * @brief One-flavor collection of even-odd preconditioned 5D ferm monomials
 */

#include "update/molecdyn/monomial/prec_constdet_one_flavor_rat_monomial5d_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"
#include "update/molecdyn/monomial/genapprox.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/invert/minvcg_array.h"

#include "actions/ferm/fermacts/prec_dwf_fermact_array_w.h"
#include "actions/ferm/fermacts/prec_ovdwf_fermact_array_w.h"
#include "actions/ferm/fermacts/prec_nef_fermact_array_w.h"
#include "actions/ferm/fermacts/prec_zolo_nef_fermact_array_w.h"
#include "actions/ferm/fermacts/prec_ovlap_contfrac5d_fermact_array_w.h"
#include "actions/ferm/fermacts/prec_ovext_fermact_array_w.h"

#if 0
#include "actions/ferm/fermacts/prec_stout_fermact_array_w.h"
#endif
namespace Chroma 
{ 
 
  namespace EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5DEnv 
  {
    //! Does the work
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomial(const string& name, 
							    XMLReader& xml, const string& path)
    {
      QDPIO::cout << "Create Fractional Monomial: " << name << endl;
      return new EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5D(
	name, OneFlavorWilsonTypeFermRatMonomial5DParams(xml, path));
    }
    
    //! Does the work
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomial(const string& name, 
							    XMLReader& xml, const string& path,
							    int expNumPower, int expDenPower) 
    {
      QDPIO::cout << "Create Monomial: " << name << endl;
      return new EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5D(
	name, OneFlavorWilsonTypeFermRatMonomial5DParams(xml, path, 
							 expNumPower, expDenPower));
    }
    

    //----------------------------------------------------------------------
    // One flavor
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialDWF1(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecDWFermActArrayEnv::name, xml, path, 1, 1);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialOvDWF1(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecOvDWFermActArrayEnv::name, xml, path, 1, 1);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialNEF1(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecNEFFermActArrayEnv::name, xml, path, 1, 1);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialZoloNEF1(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecZoloNEFFermActArrayEnv::name, xml, path, 1, 1);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialContFrac1(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecOvlapContFrac5DFermActArrayEnv::name, xml, path, 1, 1);
    }

    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialOvExt1(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecOvExtFermActArrayEnv::name, xml, path, 1, 1);
    }

#if 0
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialStout1(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecStoutWilsonTypeFermAct5DEnv::name, xml, path, 1, 1);
    }
#endif

    //----------------------------------------------------------------------
    // Three flavor
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialDWF3(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecDWFermActArrayEnv::name, xml, path, 3, 1);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialOvDWF3(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecOvDWFermActArrayEnv::name, xml, path, 3, 1);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialNEF3(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecNEFFermActArrayEnv::name, xml, path, 3, 1);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialZoloNEF3(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecZoloNEFFermActArrayEnv::name, xml, path, 3, 1);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialContFrac3(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecOvlapContFrac5DFermActArrayEnv::name, xml, path, 3, 1);
    }

    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialOvExt3(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecOvExtFermActArrayEnv::name, xml, path, 3, 1);
    }

#if 0
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialStout3(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecStoutWilsonTypeFermAct5DEnv::name, xml, path, 3, 1);
    }
#endif

    //----------------------------------------------------------------------
    // Generic fractional flavor
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialDWF(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecDWFermActArrayEnv::name, xml, path);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialOvDWF(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecOvDWFermActArrayEnv::name, xml, path);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialNEF(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecNEFFermActArrayEnv::name, xml, path);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialZoloNEF(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecZoloNEFFermActArrayEnv::name, xml, path);
    }
    
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialContFrac(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecOvlapContFrac5DFermActArrayEnv::name, xml, path);
    }

    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialOvExt(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecOvExtFermActArrayEnv::name, xml, path);
    }

#if 0
    //! Callback function for the factory
    Monomial< multi1d<LatticeColorMatrix>,
	      multi1d<LatticeColorMatrix> >* createMonomialStout(XMLReader& xml, const string& path) 
    {
      return createMonomial(EvenOddPrecStoutWilsonTypeFermAct5DEnv::name, xml, path);
    }
#endif

    //------------------------------------------------------
    //! Register one flavor
    bool registerOne(const string& prefix, const string& suffix)
    {
      bool foo = true;

      // Use a pattern to register all the qualifying fermacts
      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecDWFermActArrayEnv::name+suffix, 
							   createMonomialDWF1);

      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecOvDWFermActArrayEnv::name+suffix, 
							   createMonomialOvDWF1);

      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecNEFFermActArrayEnv::name+suffix, 
							   createMonomialNEF1);

      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecZoloNEFFermActArrayEnv::name+suffix, 
							   createMonomialZoloNEF1);

      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecOvlapContFrac5DFermActArrayEnv::name+suffix, 
							   createMonomialContFrac1);

      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecOvExtFermActArrayEnv::name+suffix, 
							   createMonomialOvExt1);

#if 0
      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecStoutWilsonTypeFermAct5DEnv::name+suffix, 
							   createMonomialStout1);
#endif

      return foo;
    }

    //------------------------------------------------------
    //! Register three flavor objects
    bool registerThree(const string& prefix, const string& suffix)
    {
      bool foo = true;

      // Use a pattern to register all the qualifying fermacts
      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecDWFermActArrayEnv::name+suffix, 
							   createMonomialDWF3);

      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecOvDWFermActArrayEnv::name+suffix, 
							   createMonomialOvDWF3);

      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecNEFFermActArrayEnv::name+suffix, 
							   createMonomialNEF3);

      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecZoloNEFFermActArrayEnv::name+suffix, 
							   createMonomialZoloNEF3);

      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecOvlapContFrac5DFermActArrayEnv::name+suffix, 
							   createMonomialContFrac3);

      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecOvExtFermActArrayEnv::name+suffix, 
							   createMonomialOvExt3);

#if 0
      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecStoutWilsonTypeFermAct5DEnv::name+suffix, 
							   createMonomialStout3);
#endif
      return foo;
    }

    //------------------------------------------------------
    //! Register generic fractional flavor
    bool registerGeneric(const string& prefix, const string& suffix)
    {
      bool foo = true;

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

      foo &= EvenOddPrecOvExtFermActArrayEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecOvExtFermActArrayEnv::name+suffix, 
							   createMonomialOvExt);

#if 0
      foo &= EvenOddPrecStoutWilsonTypeFermAct5DEnv::registered;
      foo &= TheMonomialFactory::Instance().registerObject(prefix+EvenOddPrecStoutWilsonTypeFermAct5DEnv::name+suffix, 
							   createMonomialStout);
#endif
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
  } //end namespace EvenOddPrecConstDet OneFlavorWilsonFermRatMonomialEnv


  // Constructor
  EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5D::EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5D(
    const string& name_,
    const OneFlavorWilsonTypeFermRatMonomial5DParams& param) 
  {
    inv_param = param.inv_param;
    nthRoot   = param.nthRoot;
    nthRootPV = param.nthRootPV;

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

    const FermionAction<LatticeFermion>* tmp_act = TheFermionActionFactory::Instance().createObject(fermact_string, fermact_reader, "/FermionAction");
  

    const EvenOddPrecWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >* downcast=dynamic_cast<const EvenOddPrecWilsonTypeFermAct5D< LatticeFermion, multi1d<LatticeColorMatrix> >*>(tmp_act);

    // Check success of the downcast 
    if( downcast == 0x0 ) {
      QDPIO::cerr << "Unable to downcast FermAct to EvenOddPrecWilsonTypeFermAct5D in EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5D()" << endl;
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

    QDPIO::cout << "EvenOddPrecConstDetOneFlavorWilsonTypeFermRatMonomial5D: " << fermact_string << endl;
  }


}; //end namespace Chroma


