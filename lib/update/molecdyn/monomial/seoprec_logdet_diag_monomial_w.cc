/*! \file
 *  \brief Symmetric even-odd preconditioned log(det(A_ee)) and log(det(A_oo))
 */

#include "update/molecdyn/monomial/seoprec_logdet_diag_monomial_w.h"
#include "update/molecdyn/monomial/monomial_factory.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"

namespace Chroma 
{ 

  namespace SymEvenOddPrecLogDetDiagMonomial4DEnv 
  {
    namespace
    {
      //! Callback function for the factory
      Monomial< multi1d<LatticeColorMatrix>,
		multi1d<LatticeColorMatrix> >* createMonomial(XMLReader& xml, const std::string& path) 
      {
	return new SymEvenOddPrecLogDetDiagMonomial4D(SymEvenOddPrecLogDetDiagMonomialParams(xml, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = std::string("N_FLAVOR_LOGDET_DIAG_FERM_MONOMIAL");
 
    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheMonomialFactory::Instance().registerObject(name, createMonomial);
	registered = true;
      }
      return success;
    }
  }
 

 
  SymEvenOddPrecLogDetDiagMonomialParams::SymEvenOddPrecLogDetDiagMonomialParams(XMLReader& in, const std::string& path) 
  {
    XMLReader paramtop(in, path);
    
    fermact = readXMLGroup(paramtop, "FermionAction", "FermAct");
    read(paramtop,"num_flavors", num_flavors);

    QDPIO::cout << "SymEvenOddPrecLogDetDiagMonomialParams: read \n" << fermact.id << std::endl;
  }

  void read(XMLReader& r, const std::string& path,  SymEvenOddPrecLogDetDiagMonomialParams& p) 
  {
    SymEvenOddPrecLogDetDiagMonomialParams tmp(r, path);
    p = tmp;
  }

  void write(XMLWriter& xml, const std::string& path, const SymEvenOddPrecLogDetDiagMonomialParams& p)
  {
    // Not implemented
  }


  SymEvenOddPrecLogDetDiagMonomial4D::SymEvenOddPrecLogDetDiagMonomial4D(const SymEvenOddPrecLogDetDiagMonomialParams& p) : 
    num_flavors(p.num_flavors) 
  {
    START_CODE();

    // Grok the fermact out of the XML
    std::istringstream is(p.fermact.xml);
    XMLReader fermact_reader(is);
    QDPIO::cout << "SymEvenOddPrecLogDetTwoFlavorWilsonTypeFermMonomial: construct " << p.fermact.id << std::endl;

    WilsonTypeFermAct<T,P,Q>* tmp_act = 
      TheWilsonTypeFermActFactory::Instance().createObject(p.fermact.id, fermact_reader, p.fermact.path);

    SymEvenOddPrecLogDetWilsonTypeFermAct<T,P,Q>* downcast = 
      dynamic_cast<SymEvenOddPrecLogDetWilsonTypeFermAct<T,P,Q>*>(tmp_act);

    // Check success of the downcast 
    if( downcast == 0x0 ) 
    {
      QDPIO::cerr << "Unable to downcast FermAct to SymEvenOddPrecLogDetWilsonTypeFermAct in " 
		  << __func__ << std::endl;
      QDP_abort(1);
    }
    
    fermact = downcast; 
    
    END_CODE();
  }


} // End namespace
