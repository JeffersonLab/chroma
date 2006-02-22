// $Id: qio_write_obj_funcmap.cc,v 2.3 2006-02-22 23:48:05 bjoo Exp $
/*! \file
 *  \brief Write object function map
 */

#include "named_obj.h"
#include "meas/inline/io/qio_write_obj_funcmap.h"
#include "meas/inline/io/named_objmap.h"
#include "util/ferm/eigeninfo.h"

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


    //! Write out an EigenInfo Type
    void QIOWriteEigenInfo(const string& buffer_id,
			   const string& file,
			   QDP_volfmt_t volfmt, QDP_serialparallel_t serpar)
    {

      // This is needed for QIO writing
      XMLBufferWriter file_xml, record_xml;

      // A shorthand for the object
      EigenInfo& obj=TheNamedObjMap::Instance().getData< EigenInfo >(buffer_id);

      // get the file XML and record XML out of the named buffer
      TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
      TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);

      // Write out the largest EV in BINARY so as not to loose precision
      // BUG!?? :For some reason I need to write this as a 1 element array and 
      // read it as a 1 element  array. 

      // BUG!?? I need to know my array length in advance to read back the
      //   data so I need to dump the number of Ev pairs

      // BUG!?? What is more because I cannot write primitive types using QIO 
      // (such as int or float) as it is integratedi into QDP++, I have
      // to write the number of evalues/vectors into a record XML for this
      // largest ev. 
      multi1d<Real64> largestD(1);
      largestD[0]=obj.getLargest();

      // SHorthand
      multi1d<Real> evals = obj.getEvalues();

      // Open file
      QDPFileWriter to(file_xml,file,volfmt,serpar,QDPIO_OPEN);

      // Make up an XML doc with the number of ev pairs in it
      XMLBufferWriter largest_record_xml;
      push(largest_record_xml, "NumElem");
      write(largest_record_xml, "n_vec", evals.size());
      pop(largest_record_xml);

      // Write largest EV plus XML doc
      write(to, largest_record_xml, largestD);
      
      // Convert ev-s to double if they are not that already 
      multi1d<Real64> evalsD(evals.size());

      for (int i=0; i<evals.size(); i++)
	evalsD[i]=Real64(evals[i]);

      // Write them with record XML
      write(to, record_xml, evalsD);

      // Now write the evectors 1 by 1 to avoid double storing
      multi1d<LatticeFermion>& evecs=obj.getEvectors();
      for (int i=0; i<evecs.size(); i++)
      {
	XMLBufferWriter record_xml_dummy;
	push(record_xml_dummy, "dummy_record_xml");
	pop(record_xml_dummy);
	
	// evecD=LatticeFermionD(evecs[i]);
	write(to, record_xml_dummy, evecs[i]);
      }

      // Done! That was unnecessarily painful
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
      success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("EigenInfo"), 
								    QIOWriteObjCallMap::QIOWriteEigenInfo);
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
