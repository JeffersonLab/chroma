// $Id: one_flavor_rat_monomial_params_w.cc,v 3.4 2008-05-24 04:29:31 edwards Exp $
/*! @file
 * @brief One-flavor monomial params
 */

#include "update/molecdyn/monomial/one_flavor_rat_monomial_params_w.h"

namespace Chroma 
{ 

  // Read the parameters
  OneFlavorWilsonTypeFermRatMonomialParams::OneFlavorWilsonTypeFermRatMonomialParams(XMLReader& xml_in, const string& path)
  {
    // Get the top of the parameter XML tree
    XMLReader paramtop(xml_in, path);
    
    try 
    {
      read(paramtop, "num_pf", num_pf);
      read(paramtop, "Action", numer);
    }
    catch(const string& s) 
    {
      QDPIO::cerr << "Caught Exception while reading parameters: " << s <<endl;
      QDP_abort(1);
    }
  }

  //! Read Parameters
  void read(XMLReader& xml, const std::string& path,
	    OneFlavorWilsonTypeFermRatMonomialParams& params) 
  {
    OneFlavorWilsonTypeFermRatMonomialParams tmp(xml, path);
    params = tmp;
  }

  //! Write Parameters
  void write(XMLWriter& xml, const std::string& path,
	     const OneFlavorWilsonTypeFermRatMonomialParams& params) 
  {
    push(xml, path);

    write(xml, "num_pf", params.num_pf);
    write(xml, "Action", params.numer);

    pop(xml);
  }
  
} //end namespace Chroma


