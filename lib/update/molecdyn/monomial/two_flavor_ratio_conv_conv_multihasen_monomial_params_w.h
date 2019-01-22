// -*- C++ -*-
/*! @file
 * @brief Two-flavor RatioConvConv monomial params
 */

#ifndef __TWO_FLAVOR_RATIO_CONV_CONV_MULTIHASEN_MONOMIAL_PARAMS_W_H__
#define __TWO_FLAVOR_RATIO_CONV_CONV_MULTIHASEN_MONOMIAL_PARAMS_W_H__
#include "update/molecdyn/monomial/comp_approx.h"

namespace Chroma 
{

  // Parameter structure
  /*! @ingroup monomial */
  struct TwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomialParams 
  {
    // Base Constructor
    TwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomialParams();

    // Read monomial from some root path
    TwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomialParams(XMLReader& in, const std::string&  path);
    CompActionInv_t   fermactInv;         /*!< Fermion action and invert params */
    GroupXML_t        predictor;		  /*!< The Chrono Predictor XML */
	multi1d<Real> mu;					  /*!< Shifted mass term */
	int numHasenTerms;					  /*!< Number of Hasenbusch terms */
  };

  void read(XMLReader& xml, const std::string& path, TwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomialParams& param);

  void write(XMLWriter& xml, const std::string& path, const TwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomialParams& params);

} //end namespace chroma

#endif
