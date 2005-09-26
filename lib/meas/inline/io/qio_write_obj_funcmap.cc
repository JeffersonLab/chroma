// $Id: qio_write_obj_funcmap.cc,v 2.1 2005-09-26 05:14:13 edwards Exp $
/*! \file
 *  \brief Write object function map
 */

#include "named_obj.h"
#include "meas/inline/io/qio_write_obj_funcmap.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma
{
 
  //! IO function map
  /*! \ingroup inlineio */
  namespace QIOWriteObjCallMap
  {
    //! Write a propagator
    void QIOWriteLatProp(const string& buffer_id,
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
    void QIOWriteLatPropF(const string& buffer_id,
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
    void QIOWriteLatPropD(const string& buffer_id,
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
    void QIOWriteLatFerm(const string& buffer_id,
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


#if 0
    // RGE: FOR SOME REASON, QDP CANNOT CAST A DOUBLE TO FLOATING HERE. NEED TO FIX.

    //! Write a single prec fermion
    void QIOWriteLatFermF(const string& buffer_id,
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


    //! Write a double prec fermion
    void QIOWriteLatFermD(const string& buffer_id,
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


    //! Write a gauge field in floating precision
    void QIOWriteArrayLatColMat(const string& buffer_id,
				const string& file, 
				QDP_volfmt_t volfmt, QDP_serialparallel_t serpar)
    {
      multi1d<LatticeColorMatrix> obj;
      XMLBufferWriter file_xml, record_xml;

      obj = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(buffer_id);
      TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
      TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);
    
      QDPFileWriter to(file_xml,file,volfmt,serpar,QDPIO_OPEN);
      write(to,record_xml,obj);
      close(to);
    }

    //! Write a gauge field in single precision
    void QIOWriteArrayLatColMatF(const string& buffer_id,
				 const string& file, 
				 QDP_volfmt_t volfmt, QDP_serialparallel_t serpar)
    {
      XMLBufferWriter file_xml, record_xml;

      multi1d<LatticeColorMatrix>& obj 
	= TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(buffer_id);
      multi1d<LatticeColorMatrixF> u_f(obj.size());
      for(int mu=0; mu < obj.size(); ++mu)
	u_f[mu] = obj[mu];

      TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
      TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);
    
      QDPFileWriter to(file_xml,file,volfmt,serpar,QDPIO_OPEN);
      write(to,record_xml,u_f);
      close(to);
    }

    //! Write a gauge field in double precision
    void QIOWriteArrayLatColMatD(const string& buffer_id,
				 const string& file, 
				 QDP_volfmt_t volfmt, QDP_serialparallel_t serpar)
    {
      XMLBufferWriter file_xml, record_xml;

      multi1d<LatticeColorMatrix>& obj 
	= TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(buffer_id);
      multi1d<LatticeColorMatrixD> u_f(obj.size());
      for(int mu=0; mu < obj.size(); ++mu)
	u_f[mu] = obj[mu];

      TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
      TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);
    
      QDPFileWriter to(file_xml,file,volfmt,serpar,QDPIO_OPEN);
      write(to,record_xml,u_f);
      close(to);
    }


  }  // end namespace WriteObjCallMap


  //! IO function map environment
  /*! \ingroup inlineio */
  namespace QIOWriteObjCallMapEnv
  { 
    bool registerAll(void) 
    {
      bool success = true;
      success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("LatticePropagator"), 
								    QIOWriteObjCallMap::QIOWriteLatProp);
      success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("LatticePropagatorF"), 
								    QIOWriteObjCallMap::QIOWriteLatPropF);
      success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("LatticePropagatorD"), 
								    QIOWriteObjCallMap::QIOWriteLatPropD);

      success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("LatticeFermion"), 
								    QIOWriteObjCallMap::QIOWriteLatFerm);
//      success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("LatticeFermionF"), 
//								    QIOWriteObjCallMap::QIOWriteLatFermF);
//      success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("LatticeFermionD"), 
//								 QIOWriteObjCallMap::QIOWriteLatFermD);

      success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("Multi1dLatticeColorMatrix"), 
								    QIOWriteObjCallMap::QIOWriteArrayLatColMat);
      success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("Multi1dLatticeColorMatrixF"), 
								    QIOWriteObjCallMap::QIOWriteArrayLatColMatF);
      success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("Multi1dLatticeColorMatrixD"), 
								    QIOWriteObjCallMap::QIOWriteArrayLatColMatD);
      return success;
    }

    bool registered = registerAll();
  };

}
