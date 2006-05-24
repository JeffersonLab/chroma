// -*- C++ -*-
// $Id: xml_group_reader.h,v 1.1 2006-05-24 21:07:28 edwards Exp $
/*! \file
 *  \brief Read an XML group as a string
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
  };


  //! Read group and return as a string
  /*! \ingroup io */
  GroupXML_t readXMLGroup(XMLReader& xml, 
			  const std::string& path, const std::string& type_name);

} //end namespace chroma

#endif
