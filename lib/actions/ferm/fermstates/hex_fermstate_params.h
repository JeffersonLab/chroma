// -*- C++ -*-
/* \file
 * \brief Stout params
 */

#ifndef _hex_fermstate_params_h_
#define _hex_fermstate_params_h_

#include "chromabase.h"

namespace Chroma 
{
  //! Params for hex-links
  /*! \ingroup fermstates */
  struct HexFermStateParams
  {
    //! Default constructor
    HexFermStateParams();
    HexFermStateParams(XMLReader& in, const std::string& path);

    int            n_smear;
  };

  void read(XMLReader& xml, const std::string& path, HexFermStateParams& p);
  void write(XMLWriter& xml, const std::string& path, const HexFermStateParams& p);
}

#endif
