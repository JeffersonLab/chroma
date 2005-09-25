// $Id: qio_read_obj_funcmap.cc,v 1.2 2005-09-25 20:41:09 edwards Exp $
/*! \file
 *  \brief Read object function map
 */

#include "named_obj.h"
#include "meas/inline/io/qio_read_obj_funcmap.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma
{
 
  //! IO function map
  /*! \ingroup inlineio */
  namespace QIOReadObjCallMap
  {
    //! Read a propagator
    void QIOReadLatProp(const string& buffer_id,
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
    void QIOReadLatPropF(const string& buffer_id,
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
    void QIOReadLatPropD(const string& buffer_id,
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
    void QIOReadLatFerm(const string& buffer_id,
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
    void QIOReadLatFermF(const string& buffer_id,
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
    void QIOReadLatFermD(const string& buffer_id,
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

    //! Read a propagator
    void QIOReadArrayLatColMat(const string& buffer_id,
			       const string& file, 
			       QDP_serialparallel_t serpar)
    {
      multi1d<LatticeColorMatrix> obj;
      XMLReader file_xml, record_xml;

      QDPFileReader to(file_xml,file,serpar);
      read(to,record_xml,obj);
      close(to);

      TheNamedObjMap::Instance().create< multi1d<LatticeColorMatrix> >(buffer_id);
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(buffer_id) = obj;
      TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
      TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
    }

  }  // end namespace QIOReadObjCallMap


  //! IO function map environment
  /*! \ingroup inlineio */
  namespace QIOReadObjCallMapEnv
  { 
    bool registerAll(void) 
    {
      bool success = true;
      success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("LatticePropagator"), 
								   QIOReadObjCallMap::QIOReadLatProp);
      success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("LatticePropagatorF"), 
								   QIOReadObjCallMap::QIOReadLatPropF);
      success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("LatticePropagatorD"), 
								   QIOReadObjCallMap::QIOReadLatPropD);

      success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("LatticeFermion"), 
								   QIOReadObjCallMap::QIOReadLatFerm);
      success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("LatticeFermionF"), 
								   QIOReadObjCallMap::QIOReadLatFermF);
//      success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("LatticeFermionD"), 
//								   QIOReadObjCallMap::QIOReadLatFermD);

      success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("Multi1dLatticeColorMatrix"), 
								   QIOReadObjCallMap::QIOReadArrayLatColMat);

      return success;
    }

    bool registered = registerAll();
  };

}
