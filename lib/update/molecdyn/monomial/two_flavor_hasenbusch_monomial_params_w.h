// -*- C++ -*-
// $Id: two_flavor_hasenbusch_monomial_params_w.h,v 3.0 2006-04-03 04:59:09 edwards Exp $
/*! @file
 * @brief Two-flavor Hasenbusch monomial params
 */

#ifndef __two_flavor_hasenbusch_monomial_params_w_h__
#define __two_flavor_hasenbusch_monomial_params_w_h__

#include "chromabase.h"
#include "io/param_io.h"

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
    InvertParam_t inv_param; // Inverter Parameters
    string ferm_act;
    string ferm_act_prec;
    string predictor_xml;   // The Chrono Predictor XML

  };

  void read(XMLReader& xml, const string& path, TwoFlavorHasenbuschWilsonTypeFermMonomialParams& param);

  void write(XMLWriter& xml, const string& path, const TwoFlavorHasenbuschWilsonTypeFermMonomialParams& params);

} //end namespace chroma

#endif
