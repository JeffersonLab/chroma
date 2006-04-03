// $Id: one_flavor_rat_monomial5d_params_w.cc,v 3.0 2006-04-03 04:59:09 edwards Exp $
/*! @file
 * @brief One-flavor monomial params
 */

#include "update/molecdyn/monomial/one_flavor_rat_monomial5d_params_w.h"
#include "io/param_io.h"

namespace Chroma 
{ 
 
  //! Remez input
  void read(XMLReader& xml, const string& path, 
	    OneFlavorWilsonTypeFermRatMonomial5DParams::Remez_t& input)
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
  OneFlavorWilsonTypeFermRatMonomial5DParams::OneFlavorWilsonTypeFermRatMonomial5DParams(
    XMLReader& xml_in, const string& path)
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

    QDPIO::cout << "OneFlavorWilsonTypeFermRatMonomial5DParams: read " << ferm_act << endl;
  }

  // Read the parameters
  OneFlavorWilsonTypeFermRatMonomial5DParams::OneFlavorWilsonTypeFermRatMonomial5DParams(
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

    QDPIO::cout << "OneFlavorWilsonTypeFermRatMonomial5DParams: read " << ferm_act << endl;
  }

  //! Read Parameters
  void read(XMLReader& xml, const std::string& path,
	    OneFlavorWilsonTypeFermRatMonomial5DParams& params) 
  {
    OneFlavorWilsonTypeFermRatMonomial5DParams tmp(xml, path);
    params = tmp;
  }

  //! Write Parameters
  void write(XMLWriter& xml, const std::string& path,
	     const OneFlavorWilsonTypeFermRatMonomial5DParams& params) 
  {
    // Not implemented
  }

} //end namespace Chroma


