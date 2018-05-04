// -*- C++ -*-
/*! \file
 *  \brief Read an XML group as a std::string
 */

#ifndef __xml_group_reader_h__
#define __xml_group_reader_h__

#include "chromabase.h"

namespace Chroma 
{

  //! Hold group xml and type id
  /*! \ingroup io */
  struct GroupXML_t
  {
    std::string  xml;     /*!< xml holding group */
    std::string  id;      /*!< typeid within group */
    std::string  path;    /*!< pathname of root */
  };


  //! Read group and return as a std::string
  /*! \ingroup io */
  GroupXML_t readXMLGroup(XMLReader& xml, 
			  const std::string& path, 
			  const std::string& type_name);

  //! Read group and return as a std::string
  /*! \ingroup io */
  multi1d<GroupXML_t> readXMLArrayGroup(XMLReader& xml, 
					const std::string& path, 
					const std::string& type_name);

} //end namespace chroma

#endif
