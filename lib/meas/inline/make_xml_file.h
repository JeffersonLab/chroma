// -*- C++ -*-
// $Id: make_xml_file.h,v 1.1 2005-04-19 17:11:07 edwards Exp $
/*! \file
 * \brief Make xml file writer
 */

#ifndef __make_xml_file_h__
#define __make_xml_file_h__

#include "chroma.h"

namespace Chroma
{
  //! Return a xml file writer for inline measurements
  /*! \ingroup inline */
  XMLFileWriter* makeXMLFileWriter(std::string xml_file, unsigned long update_no);
}

#endif
