// -*- C++ -*-
// $Id: stout_fermact_params.h,v 2.1 2005-10-02 03:08:49 bjoo Exp $
/* \file
 * \brief Shared things between all stout fermion actions 
 */

#ifndef _stout_fermact_params_h_
#define _stout_fermact_params_h_

#include "chromabase.h"

namespace Chroma 
{
  struct StoutFermActParams
  {
    //! Default constructor
    StoutFermActParams() {} 
    StoutFermActParams(XMLReader& in, const std::string& path);
    Real rho;
    int n_smear;
    std::string internal_fermact;
  };

  void read(XMLReader& xml, const std::string& path, StoutFermActParams& p);
  void write(XMLWriter& xml, const std::string& path, const StoutFermActParams& p);
}

#endif
