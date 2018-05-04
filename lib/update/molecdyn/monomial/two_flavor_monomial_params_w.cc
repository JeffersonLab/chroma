/*! @file
 * @brief Two-flavor monomial params
 */

#include "update/molecdyn/monomial/two_flavor_monomial_params_w.h"
#include "io/param_io.h"


namespace Chroma 
{ 
 
  // Read the parameters
  TwoFlavorWilsonTypeFermMonomialParams::TwoFlavorWilsonTypeFermMonomialParams(XMLReader& xml_in, const std::string& path)
  {
    // Get the top of the parameter XML tree
    XMLReader paramtop(xml_in, path);
    
    try {
      // Read the inverter Parameters
      inv_param = readXMLGroup(paramtop, "InvertParam", "invType");
      fermact = readXMLGroup(paramtop, "FermionAction", "FermAct");

      if( paramtop.count("./ChronologicalPredictor") == 0 ) 
      {
	predictor.xml="";
      }
      else {
	predictor = readXMLGroup(paramtop, "ChronologicalPredictor", "Name");
      }
      
    }
    catch(const std::string& s) {
      QDPIO::cerr << "Caught Exception while reading parameters: " << s <<std::endl;
      QDP_abort(1);
    }

    QDPIO::cout << "TwoFlavorWilsonTypeFermMonomialParams: read \n" << fermact.id << std::endl;
  }

  //! Read Parameters
  void read(XMLReader& xml, const std::string& path,
	    TwoFlavorWilsonTypeFermMonomialParams& params) 
  {
    TwoFlavorWilsonTypeFermMonomialParams tmp(xml, path);
    params = tmp;
  }

  //! Write Parameters
  void write(XMLWriter& xml, const std::string& path,
	     const TwoFlavorWilsonTypeFermMonomialParams& params) 
  {
    // Not implemented
  }
  
} //end namespace Chroma


