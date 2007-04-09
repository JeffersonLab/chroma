// $Id: make_xml_file.cc,v 3.1 2007-04-09 18:41:57 bjoo Exp $
/*! \file
 *  \brief Make xml file writer
 */

#include "meas/inline/make_xml_file.h"
#include <ostream>
using namespace std;

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

    std::string xml;
    // Could/should allow file pattern
    if( update_no == 0 ) { 
      xml = xml_file;
    }
    else { 
      std::ostringstream os;
      os << xml_file << "." << update_no;
      xml = os.str();
    }    
    // Return xml
    return xml;
  }
}
