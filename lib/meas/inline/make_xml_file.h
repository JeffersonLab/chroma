// -*- C++ -*-
// $Id: make_xml_file.h,v 1.2 2005-04-19 20:05:22 edwards Exp $
/*! \file
 * \brief Make xml file writer
 */

#ifndef __make_xml_file_h__
#define __make_xml_file_h__

#include "chromabase.h"

namespace Chroma
{
  //! Return a xml file name for inline measurements
  /*! \ingroup inline */
  std::string makeXMLFileName(std::string xml_file, unsigned long update_no);
}

#endif
