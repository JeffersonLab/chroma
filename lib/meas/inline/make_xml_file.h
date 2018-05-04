// -*- C++ -*-
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
