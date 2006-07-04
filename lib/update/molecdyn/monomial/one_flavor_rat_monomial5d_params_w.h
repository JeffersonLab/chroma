// -*- C++ -*-
// $Id: one_flavor_rat_monomial5d_params_w.h,v 3.2 2006-07-04 02:55:51 edwards Exp $
/*! @file
 * @brief One-flavor monomial params
 */

#ifndef __one_flavor_rat_monomial5d_params_w_h__
#define __one_flavor_rat_monomial5d_params_w_h__

#include "chromabase.h"
#include "io/xml_group_reader.h"

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

    GroupXML_t      inv_param;     /*!< Inverter Parameters */
    GroupXML_t      fermact;       /*!< Fermion action params */
    int             expNumPower;   /*!< (M^dag*M)^{expNumPower / (2*expDenPower)} */
    int             expDenPower;   /*!< (M^dag*M)^{expNumPower / (2*expDenPower)} */
    int             nthRoot;       /*!< Use "n" copies of nth-root 1-flavor */
    int             nthRootPV;     /*!< Use "n" copies of nth-root 1 flavor PV */

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
