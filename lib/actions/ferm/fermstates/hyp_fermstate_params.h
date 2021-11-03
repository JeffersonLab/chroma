// -*- C++ -*-
/* \file
 * \brief Stout params
 */

#ifndef _hyp_fermstate_params_h_
#define _hyp_fermstate_params_h_

#include "chromabase.h"

namespace Chroma 
{
  //! Params for hyp-links
  /*! \ingroup fermstates */
  struct HypFermStateParams
  {
    //! Default constructor
    HypFermStateParams();
    HypFermStateParams(XMLReader& in, const std::string& path);

    int            n_smear;
  };

  void read(XMLReader& xml, const std::string& path, HypFermStateParams& p);
  void write(XMLWriter& xml, const std::string& path, const HypFermStateParams& p);
}

#endif
