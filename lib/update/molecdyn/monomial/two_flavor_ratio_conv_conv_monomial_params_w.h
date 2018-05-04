// -*- C++ -*-
/*! @file
 * @brief Two-flavor RatioConvConv monomial params
 */

#ifndef __two_flavor_ratio_conv_conv_monomial_params_w_h__
#define __two_flavor_ratio_conv_conv_monomial_params_w_h__

#include "update/molecdyn/monomial/comp_approx.h"

namespace Chroma 
{

  // Parameter structure
  /*! @ingroup monomial */
  struct TwoFlavorRatioConvConvWilsonTypeFermMonomialParams 
  {
    // Base Constructor
    TwoFlavorRatioConvConvWilsonTypeFermMonomialParams();

    // Read monomial from some root path
    TwoFlavorRatioConvConvWilsonTypeFermMonomialParams(XMLReader& in, const std::string&  path);
    CompActionInv_t   numer;         /*!< Fermion action and invert params for numerator */
    CompActionInv_t   denom;         /*!< Fermion action for denominator */
    GroupXML_t        predictor;     /*!< The Chrono Predictor XML */
  };

  void read(XMLReader& xml, const std::string& path, TwoFlavorRatioConvConvWilsonTypeFermMonomialParams& param);

  void write(XMLWriter& xml, const std::string& path, const TwoFlavorRatioConvConvWilsonTypeFermMonomialParams& params);

} //end namespace chroma

#endif
