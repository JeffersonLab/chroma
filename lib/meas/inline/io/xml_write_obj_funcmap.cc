// $Id: xml_write_obj_funcmap.cc,v 3.1 2006-09-20 20:28:03 edwards Exp $
/*! \file
 *  \brief Write object function map
 */

#include "named_obj.h"
#include "meas/inline/io/xml_write_obj_funcmap.h"
#include "meas/inline/io/named_objmap.h"
#include "util/ferm/eigeninfo.h"

namespace Chroma
{
 
  //! IO function map environment
  /*! \ingroup inlineio */
  namespace XMLWriteObjCallMapEnv
  { 
    // Anonymous namespace
    namespace
    {
      //! Write a propagator
      void XMLWriteLatProp(const string& buffer_id,
			   const string& file)
      {
	LatticePropagator obj;
	XMLBufferWriter file_xml, record_xml;

	obj = TheNamedObjMap::Instance().getData<LatticePropagator>(buffer_id);
	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);
    
	XMLFileWriter to(file);
	push(to,"LatticePropagator");
	write(to,"FileXML",file_xml);
	write(to,"RecordXML",record_xml);
	write(to,"Object",obj);
	pop(to);
	to.close();
      }


      //! Write a fermion
      void XMLWriteLatFerm(const string& buffer_id,
			   const string& file)
      {
	LatticeFermion obj;
	XMLBufferWriter file_xml, record_xml;

	obj = TheNamedObjMap::Instance().getData<LatticeFermion>(buffer_id);
	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);
    
	XMLFileWriter to(file);
	push(to,"LatticeFermion");
	write(to,"FileXML",file_xml);
	write(to,"RecordXML",record_xml);
	write(to,"Object",obj);
	pop(to);
	to.close();
      }


      //! Write a propagator
      void XMLWriteLatStagProp(const string& buffer_id,
			       const string& file)
      {
	LatticeStaggeredPropagator obj;
	XMLBufferWriter file_xml, record_xml;

	obj = TheNamedObjMap::Instance().getData<LatticeStaggeredPropagator>(buffer_id);
	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);
    
	XMLFileWriter to(file);
	push(to,"LatticeStaggeredPropagator");
	write(to,"FileXML",file_xml);
	write(to,"RecordXML",record_xml);
	write(to,"Object",obj);
	pop(to);
	to.close();
      }


      //! Write a gauge field in floating precision
      void XMLWriteArrayLatColMat(const string& buffer_id,
				  const string& file)
      {
	multi1d<LatticeColorMatrix> obj;
	XMLBufferWriter file_xml, record_xml;

	obj = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(buffer_id);
	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);
    
	XMLFileWriter to(file);
	push(to,"Multi1dLatticeColorMatrix");
	write(to,"FileXML",file_xml);
	write(to,"RecordXML",record_xml);
	write(to,"Object",obj);
	pop(to);
	to.close();
      }

      //! Local registration flag
      bool registered = false;

    }  // end anonymous namespace


    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheXMLWriteObjFuncMap::Instance().registerFunction(string("LatticePropagator"), 
								      XMLWriteLatProp);
	success &= TheXMLWriteObjFuncMap::Instance().registerFunction(string("LatticeFermion"), 
								      XMLWriteLatFerm);
	success &= TheXMLWriteObjFuncMap::Instance().registerFunction(string("LatticeStaggeredPropagator"), 
								      XMLWriteLatStagProp);
	success &= TheXMLWriteObjFuncMap::Instance().registerFunction(string("Multi1dLatticeColorMatrix"), 
								      XMLWriteArrayLatColMat);

	registered = true;
      }
      return success;
    }
  }

}
