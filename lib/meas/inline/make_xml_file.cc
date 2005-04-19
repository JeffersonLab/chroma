// $Id: make_xml_file.cc,v 1.1 2005-04-19 17:11:07 edwards Exp $
/*! \file
 *  \brief Make xml file writer
 */

#include "meas/inline/make_xml_file.h"

namespace Chroma
{
  //! Return a xml file writer for inline measurements
  XMLFileWriter* makeXMLFileWriter(std::string xml_file, unsigned long update_no)
  {
    if (xml_file == "")
    {
      QDPIO::cerr << __func__ << ": empty xml file" << endl;
      QDP_abort(1);
    }

    // Could allow file pattern
    string xml = xml_file;

    // Make the writer and bolt
    return new XMLFileWriter(xml);
  }
}
