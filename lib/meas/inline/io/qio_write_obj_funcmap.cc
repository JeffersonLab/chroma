// $Id: qio_write_obj_funcmap.cc,v 3.14 2009-01-30 03:42:40 kostas Exp $
/*! \file
 *  \brief Write object function map
 */

#include "named_obj.h"
#include "meas/inline/io/qio_write_obj_funcmap.h"
#include "meas/inline/io/named_objmap.h"

#include "meas/hadron/diquark_w.h"
#include "util/info/unique_id.h"
#include "util/ferm/eigeninfo.h"
#include "util/ferm/subset_vectors.h"
#include "util/ferm/key_prop_colorvec.h"
#include "handle.h"
#include "qdp_map_obj_memory.h"


#include "actions/ferm/invert/containers.h"

namespace Chroma
{
 
  //! IO function map environment
  /*! \ingroup inlineio */
  namespace QIOWriteObjCallMapEnv
  { 
    // Anonymous namespace
    namespace
    {
      //------------------------------------------------------------------------
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


      //------------------------------------------------------------------------
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


      //------------------------------------------------------------------------
      //! Write a propagator
      void QIOWriteLatStagProp(const string& buffer_id,
			       const string& file, 
			       QDP_volfmt_t volfmt, QDP_serialparallel_t serpar)
      {
	LatticeStaggeredPropagator obj;
	XMLBufferWriter file_xml, record_xml;

	obj = TheNamedObjMap::Instance().getData<LatticeStaggeredPropagator>(buffer_id);
	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);
    
	QDPFileWriter to(file_xml,file,volfmt,serpar,QDPIO_OPEN);
	write(to,record_xml,obj);
	close(to);
      }


      //! Write a single prec propagator
      void QIOWriteLatStagPropF(const string& buffer_id,
				const string& file, 
				QDP_volfmt_t volfmt, QDP_serialparallel_t serpar)
      {
	LatticeStaggeredPropagatorF obj;
	XMLBufferWriter file_xml, record_xml;

	obj = TheNamedObjMap::Instance().getData<LatticeStaggeredPropagator>(buffer_id);
	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);
    
	QDPFileWriter to(file_xml,file,volfmt,serpar,QDPIO_OPEN);
	write(to,record_xml,obj);
	close(to);
      }


      //! Write a double prec propagator
      void QIOWriteLatStagPropD(const string& buffer_id,
				const string& file, 
				QDP_volfmt_t volfmt, QDP_serialparallel_t serpar)
      {
	LatticeStaggeredPropagatorD obj;
	XMLBufferWriter file_xml, record_xml;

	obj = TheNamedObjMap::Instance().getData<LatticeStaggeredPropagator>(buffer_id);
	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);
    
	QDPFileWriter to(file_xml,file,volfmt,serpar,QDPIO_OPEN);
	write(to,record_xml,obj);
	close(to);
      }


      //------------------------------------------------------------------------
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


      //------------------------------------------------------------------------
      //! Write a QQDiquark object in floating precision
      void QIOWriteQQDiquarkContract(const string& buffer_id,
				     const string& file, 
				     QDP_volfmt_t volfmt, QDP_serialparallel_t serpar)
      {
	XMLBufferWriter file_xml, record_xml;

	const QQDiquarkContract_t& obj = 
	  TheNamedObjMap::Instance().getData<QQDiquarkContract_t>(buffer_id);
	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);
    
	// QIO cannot deal with multiNd arrays, so we have to flatten to a 
	// 1-D array. This is okay since a QQDiquarkContract object has a fixed
	// size.

	const int sz_qq = Ns*Ns*Ns*Ns*Nc*Nc;
	multi1d<LatticeComplex> obj_1d(sz_qq);

	int cnt = 0;
	multi1d<int> ranks(6);
	const multi1d<int>& sz = obj.comp.size();
	for(ranks[0]=0; ranks[0] < sz[0]; ++ranks[0])
	  for(ranks[1]=0; ranks[1] < sz[1]; ++ranks[1])
	    for(ranks[2]=0; ranks[2] < sz[2]; ++ranks[2])
	      for(ranks[3]=0; ranks[3] < sz[3]; ++ranks[3])
		for(ranks[4]=0; ranks[4] < sz[4]; ++ranks[4])
		  for(ranks[5]=0; ranks[5] < sz[5]; ++ranks[5])
		  {
		    // Sanity check - the size better match
		    if (cnt >= sz_qq)
		    {
		      QDPIO::cerr << __func__ << ": size mismatch for multi1Nd object" << endl;
		      QDP_abort(1);
		    }

		    obj_1d[cnt++] = obj.comp[ranks];
		  }


	QDPFileWriter to(file_xml,file,volfmt,serpar,QDPIO_OPEN);
	write(to,record_xml,obj_1d);
	close(to);
      }


      //-----------------------------------------------------------------------
      //! Write out an EigenInfo Type
      template<typename T>
      void QIOWriteEigenInfo(const string& buffer_id,
			     const string& file,
			     QDP_volfmt_t volfmt, QDP_serialparallel_t serpar)
      {
	// This is needed for QIO writing
	XMLBufferWriter file_xml, record_xml;

	// A shorthand for the object
	EigenInfo<T>& obj=TheNamedObjMap::Instance().getData< EigenInfo<T> >(buffer_id);

	// get the file XML and record XML out of the named buffer
	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);

	// RGE: I BELIEVE ALL THE COMPLAINTS BELOW ARE NOW FIXED IN QDP++,
	// BUT I DO NOT WANT TO REMOVE THE CODE YET - CANNOT BE BOTHERED DEBUGGING

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
	largestD[0] = obj.getLargest();

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
	  evalsD[i] = Real64(evals[i]);

	// Write them with record XML
	write(to, record_xml, evalsD);

	// Now write the evectors 1 by 1 to avoid double storing
	// This is EigenInfo not SubsetVectors
	multi1d<T>& evecs=obj.getEvectors();
	for (int i=0; i<evecs.size(); i++)
	{
	  XMLBufferWriter record_xml_dummy;
	  push(record_xml_dummy, "dummy_record_xml");
	  pop(record_xml_dummy);
	
	  // evec=<T>(evecs[i]);
	  write(to, record_xml_dummy, evecs[i]);
	}

	// Done! That was unnecessarily painful
	close(to);
      }

      //------------------------------------------------------------------------
      //! Write out an RitzPairs Type
      void QIOWriteRitzPairsLatticeFermion(const string& buffer_id,
					   const string& file,
					   QDP_volfmt_t volfmt, QDP_serialparallel_t serpar)
      {
	// A shorthand for the object
	const LinAlg::RitzPairs<LatticeFermion>& obj = 
	  TheNamedObjMap::Instance().getData<LinAlg::RitzPairs< LatticeFermion> >(buffer_id);

	// File XML
	XMLBufferWriter file_xml;
	push(file_xml, "RitzPairs");
	write(file_xml, "id", uniqueId());
	write(file_xml, "Nmax", obj.evec.size());
	write(file_xml, "Neig", obj.Neig);
	pop(file_xml);

	// Open file
	QDPFileWriter to(file_xml,file,volfmt,serpar,QDPIO_OPEN);

	// Write a record for each eigenvalue (in xml) and eigenvector
	for(int i=0; i < obj.Neig; ++i)
	{
	  XMLBufferWriter record_xml;
	  push(record_xml, "Eigenvector");
	  write(record_xml, "eigenNum", i);
	  write(record_xml, "eigenValue", obj.eval.vec[i]);
	  pop(record_xml);

	  write(to, record_xml, obj.evec.vec[i]);
	}

	// Close
	close(to);
      }

      //----------------------------------------------------------------------
      void QIOWriteSubsetVectors(const string& buffer_id,
				 const string& file,
				 QDP_volfmt_t volfmt, QDP_serialparallel_t serpar)
      {
	// A shorthand for the object
	QDP::MapObject<int,EVPair<LatticeColorVector> >& obj =
	  *(TheNamedObjMap::Instance().getData< Handle< QDP::MapObject<int,EVPair<LatticeColorVector> > > >(buffer_id));

	// Yuk. Could read this back in.
	int decay_dir = Nd-1;

	// Write number of EVs to XML
	XMLBufferWriter file_xml;

	push(file_xml, "AllVectors");
	write(file_xml, "n_vec", obj.size());
	write(file_xml, "decay_dir", decay_dir);
	pop(file_xml);

	// Open file
	QDPFileWriter to(file_xml,file,volfmt,serpar,QDPIO_OPEN);

	// Loop and read evecs
	for(int n=0; n < obj.size(); n++)
	{
	  EVPair<LatticeColorVector> write_pair;
	  obj.get(n, write_pair);

	  XMLBufferWriter record_xml;
	  push(record_xml, "VectorInfo");
	  write(record_xml, "weights", write_pair.eigenValue.weights);
	  pop(record_xml);
	  write(to, record_xml, write_pair.eigenVector);
	}

	// Done
	close(to);
      }

      //------------------------------------------------------------------------
      //! Write out a MapObject Type
      template<typename K, typename V>
      void QIOWriteMapObjMemory(const string& buffer_id,
				const string& file,
				QDP_volfmt_t volfmt, QDP_serialparallel_t serpar)
      {
	// This is needed for QIO writing
	XMLBufferWriter file_xml, record_xml;

	// A shorthand for the object
	MapObjectMemory<K,V>& obj = 
	  dynamic_cast<MapObjectMemory<K,V>&>(*(TheNamedObjMap::Instance().getData< Handle<QDP::MapObject<K,V> > >(buffer_id)));

	// Get the file XML and record XML out of the named buffer
	TheNamedObjMap::Instance().get(buffer_id).getFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).getRecordXML(record_xml);

	// Open file
	QDPFileWriter to(file_xml,file,volfmt,serpar,QDPIO_OPEN);

	// Use the iterators to run through the object, saving each 
	// in a separate record
	std::vector<K> keys;
	obj.keys(keys);

	for(typename std::vector<K>::const_iterator mm = keys.begin();
	    mm != keys.end();
	    ++mm)
	{
	  XMLBufferWriter local_record_xml;
	  write(local_record_xml, "MapEntry", *mm);
	
	  write(to, local_record_xml, obj[*mm]);
	}

	// Close and bolt
	close(to);
      }

      //! Local registration flag
      bool registered = false;

    }  // end namespace WriteObjCallMap


    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("LatticePropagator"), 
								      QIOWriteLatProp);
	success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("LatticePropagatorF"), 
								      QIOWriteLatPropF);
	success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("LatticePropagatorD"), 
								      QIOWriteLatPropD);

	success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("LatticeFermion"), 
								      QIOWriteLatFerm);
//      success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("LatticeFermionF"), 
//								    QIOWriteLatFermF);
//      success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("LatticeFermionD"), 
//								 QIOWriteLatFermD);

	success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("LatticeStaggeredPropagator"), 
								      QIOWriteLatStagProp);
	success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("LatticeStaggeredPropagatorF"), 
								      QIOWriteLatStagPropF);
	success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("LatticeStaggeredPropagatorD"), 
								      QIOWriteLatStagPropD);

	success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("Multi1dLatticeColorMatrix"), 
								      QIOWriteArrayLatColMat);
	success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("Multi1dLatticeColorMatrixF"), 
								      QIOWriteArrayLatColMatF);
	success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("Multi1dLatticeColorMatrixD"), 
								      QIOWriteArrayLatColMatD);

	success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("QQDiquarkContract"), 
								      QIOWriteQQDiquarkContract);

	success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("EigenInfoLatticeFermion"), 
								      QIOWriteEigenInfo<LatticeFermion>);

	success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("RitzPairsLatticeFermion"), 
								      QIOWriteRitzPairsLatticeFermion);
	
	success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("SubsetVectorsLatticeColorVector"), 
								      QIOWriteSubsetVectors);

	success &= TheQIOWriteObjFuncMap::Instance().registerFunction(string("MapObjMemoryKeyPropColorVecLatticeFermion"), 
								      QIOWriteMapObjMemory<KeyPropColorVec_t,LatticeFermion>);

	registered = true;
      }
      return success;
    }
  }

}
