// -*- C++ -*-
// $Id: two_flavor_monomial_params_w.h,v 3.1 2006-07-03 15:26:10 edwards Exp $
/*! @file
 * @brief Two-flavor monomial params
 */

#ifndef __two_flavor_monomial_params_w_h__
#define __two_flavor_monomial_params_w_h__

#include "chromabase.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{

  // Parameter structure
  /*! @ingroup monomial */
  struct TwoFlavorWilsonTypeFermMonomialParams 
  {
    // Base Constructor
    TwoFlavorWilsonTypeFermMonomialParams();

    // Read monomial from some root path
    TwoFlavorWilsonTypeFermMonomialParams(XMLReader& in, const std::string&  path);
    GroupXML_t inv_param; // Inverter Parameters
    string ferm_act;
    string predictor_xml;   // The Chrono Predictor XML

  };

  void read(XMLReader& xml, const string& path, TwoFlavorWilsonTypeFermMonomialParams& param);

  void write(XMLWriter& xml, const string& path, const TwoFlavorWilsonTypeFermMonomialParams& params);

} //end namespace chroma

#endif
