// $Id: xml_group_reader.cc,v 1.3 2007-06-11 03:26:43 edwards Exp $
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
      group.path = "/" + path;
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << __func__ << ": caught exception reading XML: " << e << endl;
      QDP_abort(1);
    }

    return group;
  }


  // Read group and return as a string
  multi1d<GroupXML_t> readXMLArrayGroup(XMLReader& xml_in, 
					const std::string& path, 
					const std::string& type_name)
  {
    multi1d<GroupXML_t> group;

//    QDPIO::cout << __func__ << ": here is the current context XX";
//    xml_in.printCurrentContext(cout);
//    QDPIO::cout << "XX" << endl;

    try
    {
      XMLReader xml_tmp(xml_in, path);
      group.resize(xml_tmp.count("elem"));

      for(int i=0; i < group.size(); i++) 
      {
	// Create the query for the element 
	std::ostringstream element_xpath;
	element_xpath << "elem[" << (i+1) << "]";

	XMLReader xml_elem(xml_tmp, element_xpath.str());
	std::ostringstream os;
	xml_elem.print(os);
	read(xml_elem, type_name, group[i].id);
	group[i].xml = os.str();
	group[i].path = "/elem";
      }
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << __func__ << ": caught exception reading XML: " << e << endl;
      QDP_abort(1);
    }

    return group;
  }

}  // end namespace Chroma
