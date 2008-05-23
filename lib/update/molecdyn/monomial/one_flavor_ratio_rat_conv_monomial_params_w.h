// -*- C++ -*-
// $Id: one_flavor_ratio_rat_conv_monomial_params_w.h,v 3.1 2008-05-23 21:31:33 edwards Exp $
/*! @file
 * @brief One-flavor ratio of determinants rational monomial params
 */

#ifndef __one_flavor_ratio_rat_conv_monomial_params_w_h__
#define __one_flavor_ratio_rat_conv_monomial_params_w_h__

#include "update/molecdyn/monomial/comp_approx.h"

namespace Chroma 
{

  // Parameter structure
  /*! @ingroup monomial */
  struct OneFlavorWilsonTypeFermRatioRatConvMonomialParams 
  {
    // Base Constructor
    OneFlavorWilsonTypeFermRatioRatConvMonomialParams();

    // Read monomial from some root path
    OneFlavorWilsonTypeFermRatioRatConvMonomialParams(XMLReader& in, const std::string& path);

    // Params for numerator and denominator fermion actions
    CompApprox_t    numer;         /*!< Fermion action and rat. structure for numerator */
    CompAction_t    denom;         /*!< Fermion action for denominator */
    int             num_pf;        /*!< Use "num_pf" copies of pseudo-fermions for chi^dag*f(M^dag*M)*chi  */
  };

  void read(XMLReader& xml, const string& path, OneFlavorWilsonTypeFermRatioRatConvMonomialParams& param);

  void write(XMLWriter& xml, const string& path, const OneFlavorWilsonTypeFermRatioRatConvMonomialParams& params);

} //end namespace chroma

#endif
