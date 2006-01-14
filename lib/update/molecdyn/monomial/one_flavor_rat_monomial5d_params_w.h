// -*- C++ -*-
// $Id: one_flavor_rat_monomial5d_params_w.h,v 2.1 2006-01-14 05:22:32 edwards Exp $
/*! @file
 * @brief One-flavor monomial params
 */

#ifndef __one_flavor_rat_monomial5d_params_w_h__
#define __one_flavor_rat_monomial5d_params_w_h__

#include "chromabase.h"
#include "io/param_io.h"

namespace Chroma 
{

  // Parameter structure
  /*! @ingroup monomial */
  struct OneFlavorWilsonTypeFermRatMonomial5DParams 
  {
    // Base Constructor
    OneFlavorWilsonTypeFermRatMonomial5DParams();

    // Read monomial from some root path
    OneFlavorWilsonTypeFermRatMonomial5DParams(XMLReader& in, const std::string& path);
    OneFlavorWilsonTypeFermRatMonomial5DParams(XMLReader& in, const std::string& path,
							  int expNumPower_, int expDenPower_);

    InvertParam_t   inv_param;     // Inverter Parameters
    std::string     ferm_act;
    int             expNumPower;   // (M^dag*M)^{expNumPower / (2*expDenPower)}
    int             expDenPower;   // (M^dag*M)^{expNumPower / (2*expDenPower)}
    int             nthRoot;       // Use "n" copies of nth-root 1-flavor
    int             nthRootPV;     // Use "n" copies of nth-root 1 flavor PV

    struct Remez_t   // eigenvalue bounds of M^dag*M
    {
      Real lowerMin;
      Real upperMax;
      Real lowerMinPV;
      Real upperMaxPV;

      int  degree;
      int  degreePV;
      
      int digitPrecision;
    } remez;
  };

  void read(XMLReader& xml, const string& path, OneFlavorWilsonTypeFermRatMonomial5DParams& param);

  void write(XMLWriter& xml, const string& path, const OneFlavorWilsonTypeFermRatMonomial5DParams& params);

} //end namespace chroma

#endif
