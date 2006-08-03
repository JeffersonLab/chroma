// -*- C++ -*-
// $Id: stout_fermstate_params.h,v 1.2 2006-08-03 21:14:36 edwards Exp $
/* \file
 * \brief Stout params
 */

#ifndef _stout_fermstate_params_h_
#define _stout_fermstate_params_h_

#include "chromabase.h"

namespace Chroma 
{
  //! Params for stout-links
  /*! \ingroup fermacts */
  struct StoutFermStateParams
  {
    //! Default constructor
    StoutFermStateParams() {} 
    StoutFermStateParams(XMLReader& in, const std::string& path);

    multi2d<Real>  rho;
    multi1d<bool>  smear_in_this_dirP; // inelegant?
    int            n_smear;
    int            orthog_dir;
  };

  void read(XMLReader& xml, const std::string& path, StoutFermStateParams& p);
  void write(XMLWriter& xml, const std::string& path, const StoutFermStateParams& p);
}

#endif
