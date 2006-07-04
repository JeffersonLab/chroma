// $Id: two_flavor_hasenbusch_monomial_params_w.cc,v 3.2 2006-07-04 02:55:52 edwards Exp $
/*! @file
 * @brief Two-flavor Hasenbusch monomial params
 */

#include "update/molecdyn/monomial/two_flavor_hasenbusch_monomial_params_w.h"
#include "io/param_io.h"

namespace Chroma 
{ 
 
  // Read the parameters
  TwoFlavorHasenbuschWilsonTypeFermMonomialParams::TwoFlavorHasenbuschWilsonTypeFermMonomialParams(XMLReader& xml_in, const string& path)
  {
    // Get the top of the parameter XML tree
    XMLReader paramtop(xml_in, path);
    
    try {
      // Read the inverter Parameters
      inv_param = readXMLGroup(paramtop, "InvertParam", "invType");
      fermact = readXMLGroup(paramtop, "FermionAction", "FermAct");
      fermact_prec = readXMLGroup(paramtop, "FermionActionPrec", "FermAct");

      if( paramtop.count("./ChronologicalPredictor") == 0 ) 
      {
	predictor.xml = "";
      }
      else {
	predictor = readXMLGroup(paramtop, "ChronologicalPredictor", "Name");
      }
      
    }
    catch(const string& s) {
      QDPIO::cerr << "Caught Exception while reading parameters: " << s <<endl;
      QDP_abort(1);
    }

    QDPIO::cout << "TwoFlavorHasenbuschWilsonTypeFermMonomialParams: read \n" << fermact.id << endl;
  }

  //! Read Parameters
  void read(XMLReader& xml, const std::string& path,
	    TwoFlavorHasenbuschWilsonTypeFermMonomialParams& params) 
  {
    TwoFlavorHasenbuschWilsonTypeFermMonomialParams tmp(xml, path);
    params = tmp;
  }

  //! Write Parameters
  void write(XMLWriter& xml, const std::string& path,
	     const TwoFlavorHasenbuschWilsonTypeFermMonomialParams& params) 
  {
    // Not implemented
  }
 
} //end namespace Chroma


