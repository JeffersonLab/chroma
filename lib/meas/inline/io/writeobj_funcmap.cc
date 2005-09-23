// $Id: writeobj_funcmap.cc,v 1.1 2005-09-23 03:43:10 edwards Exp $
/*! \file
 *  \brief Write object function map
 */

#include "named_obj.h"
#include "meas/inline/io/writeobj_funcmap.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma
{
 
  //! IO function map
  /*! \ingroup inlineio */
  namespace WriteObjCallMap
  {
    //! Write a propagator
    void writeQIOLatProp(const string& buffer_id,
			 const string& file, 
			 QDP_volfmt_t volfmt, QDP_serialparallel_t serpar)
    {
      LatticePropagator obj;
      XMLBufferWriter file_xml, record_xml;

      obj = TheNamedObjMap::Instance().getData<LatticePropagator>(buffer_id);
      TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
      TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);
    
      QDPFileWriter to(file_xml,file,volfmt,serpar,QDPIO_OPEN);
      write(to,record_xml,obj);
      close(to);
    }


    //! Write a single prec propagator
    void writeQIOLatPropF(const string& buffer_id,
			  const string& file, 
			  QDP_volfmt_t volfmt, QDP_serialparallel_t serpar)
    {
      LatticePropagatorF obj;
      XMLBufferWriter file_xml, record_xml;

      obj = TheNamedObjMap::Instance().getData<LatticePropagator>(buffer_id);
      TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
      TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);
    
      QDPFileWriter to(file_xml,file,volfmt,serpar,QDPIO_OPEN);
      write(to,record_xml,obj);
      close(to);
    }


    //! Write a double prec propagator
    void writeQIOLatPropD(const string& buffer_id,
			  const string& file, 
			  QDP_volfmt_t volfmt, QDP_serialparallel_t serpar)
    {
      LatticePropagatorD obj;
      XMLBufferWriter file_xml, record_xml;

      obj = TheNamedObjMap::Instance().getData<LatticePropagator>(buffer_id);
      TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
      TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);
    
      QDPFileWriter to(file_xml,file,volfmt,serpar,QDPIO_OPEN);
      write(to,record_xml,obj);
      close(to);
    }


    //! Write a fermion
    void writeQIOLatFerm(const string& buffer_id,
			 const string& file, 
			 QDP_volfmt_t volfmt, QDP_serialparallel_t serpar)
    {
      LatticeFermion obj;
      XMLBufferWriter file_xml, record_xml;

      obj = TheNamedObjMap::Instance().getData<LatticeFermion>(buffer_id);
      TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
      TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);
    
      QDPFileWriter to(file_xml,file,volfmt,serpar,QDPIO_OPEN);
      write(to,record_xml,obj);
      close(to);
    }


    //! Write a single prec fermion
    void writeQIOLatFermF(const string& buffer_id,
			  const string& file, 
			  QDP_volfmt_t volfmt, QDP_serialparallel_t serpar)
    {
      LatticeFermionF obj;
      XMLBufferWriter file_xml, record_xml;

      obj = TheNamedObjMap::Instance().getData<LatticeFermion>(buffer_id);
      TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
      TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);
    
      QDPFileWriter to(file_xml,file,volfmt,serpar,QDPIO_OPEN);
      write(to,record_xml,obj);
      close(to);
    }


#if 0
    // RGE: FOR SOME REASON, QDP CANNOT CAST A DOUBLE TO FLOATING HERE. NEED TO FIX.

    //! Write a double prec fermion
    void writeQIOLatFermD(const string& buffer_id,
			  const string& file, 
			  QDP_volfmt_t volfmt, QDP_serialparallel_t serpar)
    {
      LatticeFermionD obj;
      XMLBufferWriter file_xml, record_xml;

      obj = TheNamedObjMap::Instance().getData<LatticeFermion>(buffer_id);
      TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
      TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);
    
      QDPFileWriter to(file_xml,file,volfmt,serpar,QDPIO_OPEN);
      write(to,record_xml,obj);
      close(to);
    }
#endif

  }  // end namespace WriteObjCallMap


  //! IO function map environment
  /*! \ingroup inlineio */
  namespace WriteObjCallMapEnv
  { 
    bool registerAll(void) 
    {
      bool success = true;
      success &= TheWriteObjFuncMap::Instance().registerFunction(string("LatticePropagator"), 
								 WriteObjCallMap::writeQIOLatProp);
      success &= TheWriteObjFuncMap::Instance().registerFunction(string("LatticePropagatorF"), 
								 WriteObjCallMap::writeQIOLatPropF);
      success &= TheWriteObjFuncMap::Instance().registerFunction(string("LatticePropagatorD"), 
								 WriteObjCallMap::writeQIOLatPropD);

      success &= TheWriteObjFuncMap::Instance().registerFunction(string("LatticeFermion"), 
								 WriteObjCallMap::writeQIOLatFerm);
      success &= TheWriteObjFuncMap::Instance().registerFunction(string("LatticeFermionF"), 
								 WriteObjCallMap::writeQIOLatFermF);
//      success &= TheWriteObjFuncMap::Instance().registerFunction(string("LatticeFermionD"), 
//								 WriteObjCallMap::writeQIOLatFermD);

      return success;
    }

    bool registered = registerAll();
  };

}
