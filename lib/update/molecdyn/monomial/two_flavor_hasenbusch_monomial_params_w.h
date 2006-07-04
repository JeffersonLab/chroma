// -*- C++ -*-
// $Id: two_flavor_hasenbusch_monomial_params_w.h,v 3.2 2006-07-04 02:55:52 edwards Exp $
/*! @file
 * @brief Two-flavor Hasenbusch monomial params
 */

#ifndef __two_flavor_hasenbusch_monomial_params_w_h__
#define __two_flavor_hasenbusch_monomial_params_w_h__

#include "chromabase.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{

  // Parameter structure
  /*! @ingroup monomial */
  struct TwoFlavorHasenbuschWilsonTypeFermMonomialParams 
  {
    // Base Constructor
    TwoFlavorHasenbuschWilsonTypeFermMonomialParams();

    // Read monomial from some root path
    TwoFlavorHasenbuschWilsonTypeFermMonomialParams(XMLReader& in, const std::string&  path);
    GroupXML_t inv_param; // Inverter Parameters
    GroupXML_t fermact;
    GroupXML_t fermact_prec;
    GroupXML_t predictor;   // The Chrono Predictor XML
  };

  void read(XMLReader& xml, const string& path, TwoFlavorHasenbuschWilsonTypeFermMonomialParams& param);

  void write(XMLWriter& xml, const string& path, const TwoFlavorHasenbuschWilsonTypeFermMonomialParams& params);

} //end namespace chroma

#endif
