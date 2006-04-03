// $Id: two_flavor_monomial_params_w.cc,v 3.0 2006-04-03 04:59:09 edwards Exp $
/*! @file
 * @brief Two-flavor monomial params
 */

#include "update/molecdyn/monomial/two_flavor_monomial_params_w.h"
#include "io/param_io.h"


namespace Chroma 
{ 
 
  // Read the parameters
  TwoFlavorWilsonTypeFermMonomialParams::TwoFlavorWilsonTypeFermMonomialParams(XMLReader& xml_in, const string& path)
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

    QDPIO::cout << "TwoFlavorWilsonTypeFermMonomialParams: read \n" << ferm_act << endl;
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


