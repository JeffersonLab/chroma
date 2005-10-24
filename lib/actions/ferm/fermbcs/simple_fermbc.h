// -*- C++ -*-
// $Id: simple_fermbc.h,v 2.1 2005-10-24 05:53:55 edwards Exp $
/*! \file
 *  \brief Simple fermionic BC
 */

#ifndef __simple_fermbc_h__
#define __simple_fermbc_h__

#include "fermbc.h"


namespace Chroma
{
  //! Params for simple fermbc
  /*! \ingroup fermbc */
  struct SimpleFermBCParams
  {
    SimpleFermBCParams() {}
    SimpleFermBCParams(XMLReader& in, const std::string& path);

    multi1d<int> boundary;
  };

  // Reader/writers
  /*! \ingroup fermbc */
  void read(XMLReader& xml, const std::string& path, SimpleFermBCParams& param);
  /*! \ingroup fermbc */
  void write(XMLWriter& xml, const std::string& path, const SimpleFermBCParams& param);
}

#endif
