// $Id: make_xml_file.cc,v 1.2 2005-04-19 20:05:22 edwards Exp $
/*! \file
 *  \brief Make xml file writer
 */

#include "meas/inline/make_xml_file.h"

namespace Chroma
{
  // Return a xml file name for inline measurements
  string makeXMLFileName(std::string xml_file, unsigned long update_no)
  {
    if (xml_file == "")
    {
      QDPIO::cerr << __func__ << ": empty xml file" << endl;
      QDP_abort(1);
    }

    // Could/should allow file pattern
    string xml = xml_file;

    // Return xml
    return xml;
  }
}
