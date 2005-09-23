// $Id: readobj_funcmap.cc,v 1.1 2005-09-23 03:43:09 edwards Exp $
/*! \file
 *  \brief Read object function map
 */

#include "named_obj.h"
#include "meas/inline/io/readobj_funcmap.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma
{
 
  //! IO function map
  /*! \ingroup inlineio */
  namespace ReadObjCallMap
  {
    //! Read a propagator
    void readQIOLatProp(const string& buffer_id,
			const string& file, 
			QDP_serialparallel_t serpar)
    {
      LatticePropagator obj;
      XMLReader file_xml, record_xml;

      QDPFileReader to(file_xml,file,serpar);
      read(to,record_xml,obj);
      close(to);

      TheNamedObjMap::Instance().create<LatticePropagator>(buffer_id);
      TheNamedObjMap::Instance().getData<LatticePropagator>(buffer_id) = obj;
      TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
      TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
    }


    //! Read a single prec propagator
    void readQIOLatPropF(const string& buffer_id,
			 const string& file, 
			 QDP_serialparallel_t serpar)
    {
      LatticePropagatorF obj;
      XMLReader file_xml, record_xml;

      QDPFileReader to(file_xml,file,serpar);
      read(to,record_xml,obj);
      close(to);

      TheNamedObjMap::Instance().create<LatticePropagator>(buffer_id);
      TheNamedObjMap::Instance().getData<LatticePropagator>(buffer_id) = obj;
      TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
      TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
    }


    //! Read a double prec propagator
    void readQIOLatPropD(const string& buffer_id,
			 const string& file, 
			 QDP_serialparallel_t serpar)
    {
      LatticePropagatorD obj;
      XMLReader file_xml, record_xml;

      QDPFileReader to(file_xml,file,serpar);
      read(to,record_xml,obj);
      close(to);

      TheNamedObjMap::Instance().create<LatticePropagator>(buffer_id);
      TheNamedObjMap::Instance().getData<LatticePropagator>(buffer_id) = obj;
      TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
      TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
    }



    //! Read a fermion
    void readQIOLatFerm(const string& buffer_id,
			const string& file, 
			QDP_serialparallel_t serpar)
    {
      LatticeFermion obj;
      XMLReader file_xml, record_xml;

      QDPFileReader to(file_xml,file,serpar);
      read(to,record_xml,obj);
      close(to);

      TheNamedObjMap::Instance().create<LatticeFermion>(buffer_id);
      TheNamedObjMap::Instance().getData<LatticeFermion>(buffer_id) = obj;
      TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
      TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
    }


    //! Read a single prec fermion
    void readQIOLatFermF(const string& buffer_id,
			 const string& file, 
			 QDP_serialparallel_t serpar)
    {
      LatticeFermionF obj;
      XMLReader file_xml, record_xml;

      QDPFileReader to(file_xml,file,serpar);
      read(to,record_xml,obj);
      close(to);

      TheNamedObjMap::Instance().create<LatticeFermion>(buffer_id);
      TheNamedObjMap::Instance().getData<LatticeFermion>(buffer_id) = obj;
      TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
      TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
    }


#if 0
    // RGE: FOR SOME REASON, QDP CANNOT CAST A DOUBLE TO FLOATING HERE. NEED TO FIX.

    //! Read a double prec fermion
    void readQIOLatFermD(const string& buffer_id,
			 const string& file, 
			 QDP_serialparallel_t serpar)
    {
      LatticeFermionD obj;
      XMLReader file_xml, record_xml;

      QDPFileReader to(file_xml,file,serpar);
      read(to,record_xml,obj);
      close(to);

      TheNamedObjMap::Instance().create<LatticeFermion>(buffer_id);
      TheNamedObjMap::Instance().getData<LatticeFermion>(buffer_id) = obj;
      TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
      TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
    }
#endif

  }  // end namespace ReadObjCallMap


  //! IO function map environment
  /*! \ingroup inlineio */
  namespace ReadObjCallMapEnv
  { 
    bool registerAll(void) 
    {
      bool success = true;
      success &= TheReadObjFuncMap::Instance().registerFunction(string("LatticePropagator"), 
								ReadObjCallMap::readQIOLatProp);
      success &= TheReadObjFuncMap::Instance().registerFunction(string("LatticePropagatorF"), 
								ReadObjCallMap::readQIOLatPropF);
      success &= TheReadObjFuncMap::Instance().registerFunction(string("LatticePropagatorD"), 
								ReadObjCallMap::readQIOLatPropD);

      success &= TheReadObjFuncMap::Instance().registerFunction(string("LatticeFermion"), 
								ReadObjCallMap::readQIOLatFerm);
      success &= TheReadObjFuncMap::Instance().registerFunction(string("LatticeFermionF"), 
								ReadObjCallMap::readQIOLatFermF);
//      success &= TheReadObjFuncMap::Instance().registerFunction(string("LatticeFermionD"), 
//								ReadObjCallMap::readQIOLatFermD);

      return success;
    }

    bool registered = registerAll();
  };

}
