// -*- C++ -*-
// $Id: one_flavor_ratio_rat_rat_monomial_params_w.h,v 3.1 2008-05-23 21:31:34 edwards Exp $
/*! @file
 * @brief One-flavor ratio of determinants rational monomial params
 */

#ifndef __one_flavor_ratio_rat_rat_monomial_params_w_h__
#define __one_flavor_ratio_rat_rat_monomial_params_w_h__

#include "update/molecdyn/monomial/comp_approx.h"

namespace Chroma 
{

  // Parameter structure
  /*! @ingroup monomial */
  struct OneFlavorWilsonTypeFermRatioRatRatMonomialParams 
  {
    // Base Constructor
    OneFlavorWilsonTypeFermRatioRatRatMonomialParams();

    // Read monomial from some root path
    OneFlavorWilsonTypeFermRatioRatRatMonomialParams(XMLReader& in, const std::string& path);

    // Params for numerator and denominator fermion actions
    CompApprox_t    numer;         /*!< Fermion action and rat. structure for numerator */
    CompApprox_t    denom;         /*!< Fermion action and rat. structure for denominator */
    int             num_pf;        /*!< Use "num_pf" copies of pseudo-fermions for chi^dag*f(M^dag*M)*chi  */
  };

  /*! @ingroup monomial */
  void read(XMLReader& xml, const string& path, OneFlavorWilsonTypeFermRatioRatRatMonomialParams& param);

  /*! @ingroup monomial */
  void write(XMLWriter& xml, const string& path, const OneFlavorWilsonTypeFermRatioRatRatMonomialParams& params);

} //end namespace chroma

#endif
