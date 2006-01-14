// $Id: two_flavor_hasenbusch_monomial_params_w.cc,v 2.1 2006-01-14 05:22:32 edwards Exp $
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
      read(paramtop, "./InvertParam", inv_param);
      XMLReader xml_tmp(paramtop, "./FermionAction");
      std::ostringstream os;
      xml_tmp.print(os);
      ferm_act = os.str();

      XMLReader xml_tmp_prec(paramtop, "./FermionActionPrec");
      std::ostringstream os_prec;
      xml_tmp_prec.print(os_prec);
      ferm_act_prec = os_prec.str();

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

    QDPIO::cout << "TwoFlavorHasenbuschWilsonTypeFermMonomialParams: read \n" << ferm_act << endl;
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


