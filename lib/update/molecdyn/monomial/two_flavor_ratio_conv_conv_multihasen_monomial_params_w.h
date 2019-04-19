// -*- C++ -*-
/*! @file
 * @brief Two-flavor RatioConvConv monomial params
 */

#ifndef __TWO_FLAVOR_RATIO_CONV_CONV_MULTIHASEN_MONOMIAL_PARAMS_W_H__
#define __TWO_FLAVOR_RATIO_CONV_CONV_MULTIHASEN_MONOMIAL_PARAMS_W_H__
#include "update/molecdyn/monomial/comp_approx.h"
#include "actions/ferm/invert/syssolver_mrhs_twisted_params.h"

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
	SysSolverMRHSTwistedParams inverter;  /*!< Multi-right hand solver parameter */
	GroupXML_t fermact;					  /*!< Fermion action */
	GroupXML_t predictor;				  /*!< The Chrono Predictor XML */
	Real base_twist;					  /*!< Shifted mass of base operator */
	int numHasenTerms;					  /*!< Number of Hasenbusch terms */
  };

  void read(XMLReader& xml, const std::string& path, TwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomialParams& param);

  void write(XMLWriter& xml, const std::string& path, const TwoFlavorRatioConvConvMultihasenWilsonTypeFermMonomialParams& params);

} //end namespace chroma

#endif
