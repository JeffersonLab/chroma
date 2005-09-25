// -*- C++ -*-
// $Id: simple_fermbc.h,v 2.0 2005-09-25 21:04:27 edwards Exp $
/*! \file
 *  \brief Simple fermionic BC
 */

#ifndef __simple_fermbc_h__
#define __simple_fermbc_h__

#include "fermbc.h"


namespace Chroma
{
  //! Params for simple fermbc
  struct SimpleFermBCParams
  {
    SimpleFermBCParams() {}
    SimpleFermBCParams(XMLReader& in, const std::string& path);

    multi1d<int> boundary;
  };

  //! Name
  namespace SimpleFermBCEnv
  {
    extern const std::string name;
  };


  // Reader/writers
  void read(XMLReader& xml, const string& path, SimpleFermBCParams& param);
  void write(XMLWriter& xml, const string& path, const SimpleFermBCParams& param);
}

#endif
