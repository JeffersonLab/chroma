// $Id: default_gauge_field.cc,v 3.4 2007-09-23 15:42:04 edwards Exp $
/*! \file
 * \brief Functions to set and get default gauge field
 */

#include "chromabase.h"

#include "meas/inline/io/named_objmap.h"
#include "meas/inline/io/default_gauge_field.h"

#include <stdio.h>
#include <typeinfo>

namespace Chroma 
{ 
  namespace InlineDefaultGaugeField
  {
    //! Private namespace
    /*! \ingroup inlineio */
    namespace
    { 
      bool initP = false;
      const std::string private_id = "default_gauge_field";
//      std::string private_id = "Supercalifragilisticexpialidocious";
    }


    //! Get the default gauge field named object id
    /*! \ingroup inlineio */
    std::string getId()
    {
      return private_id;
    }
  

    //! Reset the default gauge field state
    /*! \ingroup inlineio */
    void reset()
    {
      START_CODE();

      // If initialized, then destroy and reset
      if (initP)
      {
	// Only delete if the object is still present
	if ( TheNamedObjMap::Instance().check(private_id) )
	{
	  try
	  {
	    // Now erase the object
	    TheNamedObjMap::Instance().erase(private_id);
	  }
	  catch( std::bad_cast ) 
	  {
	    QDPIO::cerr << __func__ << ": cast error" 
			<< endl;
	    QDP_abort(1);
	  }
	  catch (const string& e) 
	  {
	    QDPIO::cerr << __func__ << ": error message: " << e 
			<< endl;
	    QDP_abort(1);
	  }
	}
      }

      initP = false;

      END_CODE();
    }
  


    //! Set the default gauge field
    /*! \ingroup inlineio */
    void set(const multi1d<LatticeColorMatrix>& u,
	     XMLBufferWriter& record_xml)
    {
      START_CODE();

      // This might be too restrictive - might prefer automatic resetting
      if (initP)
      {
	QDPIO::cerr << __func__ << ": gauge field already initialized" << endl;
	QDP_abort(1);
      }

      // Generate a private file name suitable for a buffer name
//      private_id = tmpnam(NULL);

      // Set the default gauge field
      try
      {
	XMLBufferWriter file_xml;
	push(file_xml, "gauge");
	write(file_xml, "id", int(0));
	pop(file_xml);

	// Save the gauge field in this private buffer
	TheNamedObjMap::Instance().create< multi1d<LatticeColorMatrix> >(private_id);
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(private_id) = u;
	TheNamedObjMap::Instance().get(private_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(private_id).setRecordXML(record_xml);
      }
      catch (std::bad_cast) 
      {
	QDPIO::cerr << __func__ << ": cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << __func__ << ": error message: " << e 
		    << endl;
	QDP_abort(1);
      }
    
      initP = true;

      END_CODE();
    }


  
    //! Get the default gauge field
    /*! \ingroup inlineio */
    void get(multi1d<LatticeColorMatrix>& u,
	     XMLBufferWriter& file_xml,
	     XMLBufferWriter& record_xml)
    {
      START_CODE();

      if (! initP)
      {
	QDPIO::cerr << __func__ << ": gauge field not initialized" << endl;
	QDP_abort(1);
      }

      // Extract the default gauge field from a named object
      try
      {
	u = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(private_id);
	TheNamedObjMap::Instance().get(private_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(private_id).setRecordXML(record_xml);
      }
      catch (std::bad_cast) 
      {
	QDPIO::cerr << __func__ << ": cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << __func__ << ": error message: " << e 
		    << endl;
	QDP_abort(1);
      }
    
      END_CODE();
    }
  
  }

}
