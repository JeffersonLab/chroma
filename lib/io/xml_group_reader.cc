// $Id: xml_group_reader.cc,v 1.1 2006-05-24 21:07:27 edwards Exp $
/*! \file
 *  \brief Read an XML group as a string
 */

#include "io/xml_group_reader.h"

namespace Chroma 
{

  // Read group and return as a string
  GroupXML_t readXMLGroup(XMLReader& xml_in, 
			  const std::string& path, const std::string& type_name)
  {
    GroupXML_t group;

//    QDPIO::cout << __func__ << ": here is the current context XX";
//    xml_in.printCurrentContext(cout);
//    QDPIO::cout << "XX" << endl;

    try
    {
      XMLReader xml_tmp(xml_in, path);
      std::ostringstream os;
      xml_tmp.print(os);
      read(xml_tmp, type_name, group.id);
      group.xml = os.str();
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << __func__ << ": caught exception reading XML: " << e << endl;
      QDP_abort(1);
    }

    return group;
  }

}  // end namespace Chroma
