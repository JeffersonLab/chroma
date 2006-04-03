// -*- C++ -*-
// $Id: one_flavor_rat_monomial_params_w.h,v 3.0 2006-04-03 04:59:09 edwards Exp $
/*! @file
 * @brief One-flavor monomial params
 */

#ifndef __one_flavor_rat_monomial_params_w_h__
#define __one_flavor_rat_monomial_params_w_h__

#include "chromabase.h"
#include "invtype.h"

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
    OneFlavorWilsonTypeFermRatMonomialParams(XMLReader& in, const std::string& path,
					     int expNumPower_, int expDenPower_);

    InvertParam_t   inv_param;     // Inverter Parameters
    std::string     ferm_act;
    int             expNumPower;   // (M^dag*M)^{expNumPower / (2*expDenPower)}
    int             expDenPower;   // (M^dag*M)^{expNumPower / (2*expDenPower)}
    int             nthRoot;       // Use "n" copies of nth-root 1-flavor

    struct Remez_t   // eigenvalue bounds of M^dag*M
    {
      Real lowerMin;
      Real upperMax;
      int  forceDegree;
      int  actionDegree;
      int  digitPrecision;
    } remez;
  };

  void read(XMLReader& xml, const string& path, OneFlavorWilsonTypeFermRatMonomialParams& param);

  void write(XMLWriter& xml, const string& path, const OneFlavorWilsonTypeFermRatMonomialParams& params);

} //end namespace chroma

#endif
