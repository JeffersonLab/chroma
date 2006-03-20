// -*- C++ -*-
// $Id: default_gauge_field.h,v 2.1 2006-03-20 04:22:03 edwards Exp $
/*! \file
 * \brief Functions to set and get default gauge field
 */

#ifndef __default_gauge_field_h__
#define __default_gauge_field_h__

#include "chromabase.h"

namespace Chroma 
{ 
  //! Namespace to support default gauge field manipulations
  /*! \ingroup inlineio */
  namespace InlineDefaultGaugeField
  {
    //! Set the default gauge field
    /*! \ingroup inlineio */
    void set(const multi1d<LatticeColorMatrix>& u,
	     XMLBufferWriter& gauge_xml);
  
    //! Get the default gauge field
    /*! \ingroup inlineio */
    void get(multi1d<LatticeColorMatrix>& u,
	     XMLBufferWriter& file_xml,
	     XMLBufferWriter& record_xml);
  
    //! Get the default gauge field named object id
    /*! \ingroup inlineio */
    std::string getId();
 
    //! Helper function to read the Id from an XML input
    /*! 
     * \ingroup inlineio 
     *
     * \return either the input id or if missing then the default gauge id
     */
    std::string readGaugeId(XMLReader& xml_in, const std::string path);
  
    //! Reset the default gauge field state
    /*! \ingroup inlineio */
    void reset();
  }  
  
}

#endif
