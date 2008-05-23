// -*- C++ -*-
// $Id: two_flavor_ratio_conv_rat_monomial_params_w.h,v 3.1 2008-05-23 21:31:34 edwards Exp $
/*! @file
 * @brief Two-flavor ratio of conventional fermion action monomial params
 */

#ifndef __two_flavor_ratio_conv_rat_monomial_params_w_h__
#define __two_flavor_ratio_conv_rat_monomial_params_w_h__

#include "update/molecdyn/monomial/comp_approx.h"

namespace Chroma 
{

  // Parameter structure
  /*! @ingroup monomial */
  struct TwoFlavorRatioConvRatWilsonTypeFermMonomialParams 
  {
    // Base Constructor
    TwoFlavorRatioConvRatWilsonTypeFermMonomialParams();

    // Read monomial from some root path
    TwoFlavorRatioConvRatWilsonTypeFermMonomialParams(XMLReader& in, const std::string&  path);
    CompActionInv_t   numer;         /*!< Fermion action and invert params for numerator */
    CompApprox_t      denom;         /*!< Fermion action and rat. structure for denominator */
    GroupXML_t        predictor;     /*!< The Chrono Predictor XML */
  };

  /*! @ingroup monomial */
  void read(XMLReader& xml, const string& path, TwoFlavorRatioConvRatWilsonTypeFermMonomialParams& param);

  /*! @ingroup monomial */
  void write(XMLWriter& xml, const string& path, const TwoFlavorRatioConvRatWilsonTypeFermMonomialParams& params);

} //end namespace chroma

#endif
