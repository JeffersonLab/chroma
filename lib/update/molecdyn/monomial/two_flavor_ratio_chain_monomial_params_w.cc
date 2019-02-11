/*! @file
 * @brief Two-flavor RatioChain monomial params
 */

#include "update/molecdyn/monomial/two_flavor_ratio_chain_monomial_params_w.h"
#include "io/param_io.h"

namespace Chroma 
{ 
 
  // Read the parameters
  TwoFlavorRatioChainWilsonTypeFermMonomialParams::TwoFlavorRatioChainWilsonTypeFermMonomialParams(XMLReader& xml_in, const std::string& path)
  {
    // Get the top of the parameter XML tree
    XMLReader paramtop(xml_in, path);
    
    try {
      // Read the inverter Parameters
      read(paramtop, "Action", fermact);

      if( paramtop.count("./ChronologicalPredictor") == 0 ) 
      {
	predictor.xml = "";
      }
      else {
	predictor = readXMLGroup(paramtop, "ChronologicalPredictor", "Name");
      }
      read(paramtop, "mu", mu);      
    }
    catch(const std::string& s) {
      QDPIO::cerr << "Caught Exception while reading parameters: " << s <<std::endl;
      QDP_abort(1);
    }
  }

  //! Read Parameters
  void read(XMLReader& xml, const std::string& path,
	    TwoFlavorRatioChainWilsonTypeFermMonomialParams& params) 
  {
    TwoFlavorRatioChainWilsonTypeFermMonomialParams tmp(xml, path);
    params = tmp;
  }

  //! Write Parameters
  void write(XMLWriter& xml, const std::string& path,
	     const TwoFlavorRatioChainWilsonTypeFermMonomialParams& params) 
  {
    write(xml, "Action", params.fermact);
    write(xml, "mu", params.mu);

    xml << params.predictor.xml;
  }
 
} //end namespace Chroma


