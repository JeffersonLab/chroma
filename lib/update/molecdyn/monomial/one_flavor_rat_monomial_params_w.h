// -*- C++ -*-
// $Id: one_flavor_rat_monomial_params_w.h,v 3.3 2008-05-23 21:31:33 edwards Exp $
/*! @file
 * @brief One-flavor monomial params
 */

#ifndef __one_flavor_rat_monomial_params_w_h__
#define __one_flavor_rat_monomial_params_w_h__

#include "chromabase.h"
#include "update/molecdyn/monomial/comp_approx.h"

namespace Chroma 
{

  // Parameter structure
  /*! @ingroup monomial */
  struct OneFlavorWilsonTypeFermRatMonomialParams 
  {
    // Base Constructor
    OneFlavorWilsonTypeFermRatMonomialParams();

    // Read monomial from some root path
    OneFlavorWilsonTypeFermRatMonomialParams(XMLReader& in, const std::string& path);

    // Params for each major group - action/heatbath & force
    CompApprox_t    numer;         /*!< Fermion action and rat. structure for numerator */
    int             num_pf;        /*!< Use "num_pf" copies of pseudo-fermions for chi^dag*f(M^dag*M)*chi  */
  };

  void read(XMLReader& xml, const string& path, OneFlavorWilsonTypeFermRatMonomialParams& param);

  void write(XMLWriter& xml, const string& path, const OneFlavorWilsonTypeFermRatMonomialParams& params);

} //end namespace chroma

#endif
