// -*- C++ -*-
/*! @file
 * @brief Two-flavor RatioConvConv monomial params
 */

#ifndef __two_flavor_ratio_chain_monomial_params_w_h__
#define __two_flavor_ratio_chain_monomial_params_w_h__

#include "update/molecdyn/monomial/comp_approx.h"

namespace Chroma 
{

  // Parameter structure
  /*! @ingroup monomial */
  struct TwoFlavorRatioChainWilsonTypeFermMonomialParams 
  {
    // Base Constructor
    TwoFlavorRatioChainWilsonTypeFermMonomialParams();

    // Read monomial from some root path
    TwoFlavorRatioChainWilsonTypeFermMonomialParams(XMLReader& in, const std::string&  path);

    CompActionInv_t   fermact;         /*!< Fermion action and invert params for numerator */
    multi1d<Real>     mu;
    GroupXML_t        predictor;     /*!< The Chrono Predictor XML */
  };

  void read(XMLReader& xml, const std::string& path, TwoFlavorRatioChainWilsonTypeFermMonomialParams& param);

  void write(XMLWriter& xml, const std::string& path, const TwoFlavorRatioChainWilsonTypeFermMonomialParams& params);

} //end namespace chroma

#endif
