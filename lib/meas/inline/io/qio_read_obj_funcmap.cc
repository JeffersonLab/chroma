// $Id: qio_read_obj_funcmap.cc,v 3.16 2009-01-30 03:42:39 kostas Exp $
/*! \file
 *  \brief Read object function map
 */

#include "named_obj.h"
#include "meas/inline/io/qio_read_obj_funcmap.h"
#include "meas/inline/io/named_objmap.h"

#include "meas/hadron/diquark_w.h"
#include "util/ferm/eigeninfo.h"
#include "util/ferm/subset_vectors.h"
#include "util/ferm/key_prop_colorvec.h"
#include "util/ferm/key_grid_prop.h"
#include "util/ferm/key_block_prop.h"
#include "handle.h"
#include "util/ferm/map_obj.h"
#include "util/ferm/map_obj/map_obj_memory.h"
#include "actions/ferm/invert/containers.h"

namespace Chroma
{
 
  //! IO function map environment
  /*! \ingroup inlineio */
  namespace QIOReadObjCallMapEnv
  { 
    // Anonymous namespace
    namespace
    {
      //------------------------------------------------------------------------
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



      //------------------------------------------------------------------------
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


#if 0
      // RGE: FOR SOME REASON, QDP CANNOT CAST A DOUBLE TO FLOATING HERE. NEED TO FIX.

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

      //------------------------------------------------------------------------
      //! Read a gauge field in floating precision
      void QIOReadArrayLatColMat(const string& buffer_id,
				 const string& file, 
				 QDP_serialparallel_t serpar)
      {
	multi1d<LatticeColorMatrix> obj;
	XMLReader file_xml, record_xml;

	obj.resize(Nd);    // BAD BAD BAD - FIX THIS

	QDPFileReader to(file_xml,file,serpar);
	read(to,record_xml,obj);
	close(to);

	TheNamedObjMap::Instance().create< multi1d<LatticeColorMatrix> >(buffer_id);
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(buffer_id) = obj;
	TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
      }

      //! Read a single prec Gauge fields
      void QIOReadArrayLatColMatF(const string& buffer_id,
				  const string& file, 
				  QDP_serialparallel_t serpar)
      {
	multi1d<LatticeColorMatrixF> obj;
	XMLReader file_xml, record_xml;

	obj.resize(Nd);    // BAD BAD BAD - FIX THIS

	QDPFileReader to(file_xml,file,serpar);
	read(to,record_xml,obj);
	close(to);

	TheNamedObjMap::Instance().create< multi1d<LatticeColorMatrix> >(buffer_id);
	(TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(buffer_id)).resize(Nd);
	for(int i=0; i < Nd; i++) {
	    (TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(buffer_id))[i] = obj[i];
	}
	TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
      }

      //! Read a Double Prec Gauge Field
      void QIOReadArrayLatColMatD(const string& buffer_id,
				  const string& file, 
				  QDP_serialparallel_t serpar)
      {
	multi1d<LatticeColorMatrixD> obj;
	XMLReader file_xml, record_xml;

	obj.resize(Nd);    // BAD BAD BAD - FIX THIS

	QDPFileReader to(file_xml,file,serpar);
	read(to,record_xml,obj);
	close(to);

	TheNamedObjMap::Instance().create< multi1d<LatticeColorMatrix> >(buffer_id);
	(TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(buffer_id)).resize(Nd);
	for(int i=0; i < Nd;i++) {
	    (TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(buffer_id))[i] = obj[i];
	}

	TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
      }


      //------------------------------------------------------------------------
      //! Read a QQDiquark object in floating precision
      void QIOReadQQDiquarkContract(const string& buffer_id,
				    const string& file, 
				    QDP_serialparallel_t serpar)
      {
	// Read the 1-d flattened array version
	XMLReader file_xml, record_xml;
	const int sz_qq = Ns*Ns*Ns*Ns*Nc*Nc;
	multi1d<LatticeComplex> obj_1d(sz_qq);

	QDPFileReader to(file_xml,file,serpar);
	read(to,record_xml,obj_1d);
	close(to);

	// Create the named object
	TheNamedObjMap::Instance().create<QQDiquarkContract_t>(buffer_id);
	TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);

	// Convert the 1-d object into an N-d object
	QQDiquarkContract_t& obj = TheNamedObjMap::Instance().getData<QQDiquarkContract_t>(buffer_id);
	multi1d<int> sz(6);
	sz = Ns;
	sz[4] = Nc;  // cf
	sz[5] = Nc;  // ci
	obj.comp.resize(sz);

	int cnt = 0;
	multi1d<int> ranks(6);
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

		    obj.comp[ranks] = obj_1d[cnt++];
		  }

      }


      //------------------------------------------------------------------------
      void QIOReadEigenInfoLatticeFermion(const string& buffer_id,
			    const string& file,
			    QDP_serialparallel_t serpar)
      {
	multi1d<Real64> evalsD;

	// RGE: I BELIEVE ALL THE COMPLAINTS BELOW ARE NOW FIXED IN QDP++,
	// BUT I DO NOT WANT TO REMOVE THE CODE YET - CANNOT BE BOTHERED DEBUGGING

	// BUG? Need to read these as an array even though there is only one
	// and also I need to know in advance the array size.
	multi1d<Real64> largestD(1);

	// Would like this to be double, but casting double prec to single prec doesn't work.
	LatticeFermion evecD;

	// Stuff 
	XMLReader file_xml, record_xml, largest_record_xml;

	// Open file
	QDPFileReader to(file_xml,file,serpar);

	// Create object
	TheNamedObjMap::Instance().create< EigenInfo<LatticeFermion> >(buffer_id);

	// Read largest ev plus XML file containing number of eval/evec pairs
	read(to, largest_record_xml, largestD);

	// Extract number of EVs from XML
	int size;
	try { 
	  read(largest_record_xml, "/NumElem/n_vec", size);
	}
	catch(const std::string& e) { 
	  QDPIO::cerr<< "Caught Exception while reading XML: " << e << endl;
	  QDP_abort(1);
	}

	// Set largest EV
	TheNamedObjMap::Instance().getData< EigenInfo<LatticeFermion> >(buffer_id).getLargest()=largestD[0];

	// REsize eval arrays so that IO works correctly
	evalsD.resize(size);

	// Read evals
	read(to, record_xml, evalsD);

	QDPIO::cout << "Read " << evalsD.size() << "Evals " << endl;
	for(int i=0; i < evalsD.size(); i++) { 
	  QDPIO::cout << "Eval["<<i<<"] = " << Real(evalsD[i]) << endl;
	}

	// Resize eval array 
	TheNamedObjMap::Instance().getData< EigenInfo<LatticeFermion> >(buffer_id).getEvalues().resize(evalsD.size());
      
	// Downcast evals to Real() precision
	for (int i=0; i<evalsD.size(); i++)
	  TheNamedObjMap::Instance().getData< EigenInfo<LatticeFermion> >(buffer_id).getEvalues()[i] = Real(evalsD[i]);


	// Resize evec arrays
	// THIS IS AN EIGEN INFO - NOT SUBSET VECTORS
	TheNamedObjMap::Instance().getData< EigenInfo<LatticeFermion> >(buffer_id).getEvectors().resize(evalsD.size());

	// Loop and read evecs
	for (int i=0; i<evalsD.size(); i++)
	{
	  XMLReader record_xml_dummy;
	  read(to, record_xml_dummy, evecD);
	  LatticeFermion evecR;
	  evecR=evecD;

	  TheNamedObjMap::Instance().getData< EigenInfo<LatticeFermion> >(buffer_id).getEvectors()[i]=evecR;
	}
 
	// Set File and Record XML throw away dummy XMLs
	TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);

	// Done - That too was unnecessarily painful
	close(to);
      }


      //------------------------------------------------------------------------
      template<typename T>
      void QIOReadSubsetVectors(const string& buffer_id,
				const string& file,
				QDP_serialparallel_t serpar)
      {
	// Stuff 
	XMLReader file_xml, record_xml, largest_record_xml;

	// Open file
	QDPFileReader to(file_xml,file,serpar);

	// Create object
	TheNamedObjMap::Instance().create< SubsetVectors<T> >(buffer_id);
	SubsetVectors<T>& obj = TheNamedObjMap::Instance().getData< SubsetVectors<T> >(buffer_id);

	// Extract number of EVs from XML
	int N, decay_dir;
	try { 
	  XMLReader vec_xml(file_xml, "/AllVectors");

	  read(vec_xml, "n_vec", N);
	  read(vec_xml, "decay_dir", decay_dir);
	}
	catch(const std::string& e) { 
	  QDPIO::cerr<< "Caught Exception while reading XML: " << e << endl;
	  QDP_abort(1);
	}

	// Resize arrays
	const int Lt = QDP::Layout::lattSize()[decay_dir];
	obj.resizeEvectors(N);
	obj.resizeEvalues(N);

	// Loop and read evecs
	for(int n=0; n < N; n++)
	{
	  XMLReader record_xml_dummy;

	  read(to, record_xml_dummy, obj.getEvector(n));
	  read(record_xml_dummy, "/VectorInfo/weights", obj.getEvalue(n).weights);
	}
 
	// Set File and Record XML throw away dummy XMLs
	TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);

	// Done
	close(to);
      }

      //-----------------------------------------------------------------------
      //! Read a RitzPairs Type
      void QIOReadRitzPairsLatticeFermion(const string& buffer_id,
					  const string& file,
					  QDP_serialparallel_t serpar)
      {
	// File XML
	XMLReader file_xml;

	// Open file
	QDPFileReader to(file_xml,file,serpar);

	// Create the named object
	TheNamedObjMap::Instance().create< LinAlg::RitzPairs<LatticeFermion> >(buffer_id);
	TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);

	// A shorthand for the object
	LinAlg::RitzPairs<LatticeFermion>& obj = 
	  TheNamedObjMap::Instance().getData<LinAlg::RitzPairs< LatticeFermion> >(buffer_id);

	XMLReader record_xml;
	TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
	
	int Nmax;
	read(file_xml, "/RitzPairs/Nmax", Nmax);
	read(file_xml, "/RitzPairs/Neig", obj.Neig);

	obj.evec.resize(Nmax);
	obj.eval.resize(Nmax);

	if (obj.Neig > Nmax)
	{
	  QDPIO::cerr << __func__ << ": error, found Neig > Nmax" << endl;
	  QDP_abort(1);
	}

	// Read a record for each eigenvalue (in xml) and eigenvector
	for(int i=0; i < obj.Neig; ++i)
	{
	  XMLReader record_xml;
	  read(to, record_xml, obj.evec.vec[i]);

	  read(record_xml, "/Eigenvector/eigenValue", obj.eval.vec[i]);
	}

	// Close
	close(to);
      }

      //! Read a MapObject Type
      template<typename K, typename V>
      void QIOReadMapObjMemory(const string& buffer_id,
			       const string& file,
			       QDP_serialparallel_t serpar)
      {
	// This is needed for QIO reading
	XMLReader file_xml;

	// Open file
	QDPFileReader to(file_xml,file,serpar);

	// Create object if it does not exist. Otherwise, we can continue reading into it.
	if (! TheNamedObjMap::Instance().check(buffer_id))
        {
	  // Make a memory map object -- later make with factory?
	  Handle< MapObject<K,V> > obj_map_handle( new MapObjectMemory<K,V> );
 
	  // Create slot
	  TheNamedObjMap::Instance().create< Handle<MapObject<K,V> > >(buffer_id);
	  // Insert Handle
	  TheNamedObjMap::Instance().getData< Handle<MapObject<K,V> > >(buffer_id) = obj_map_handle;

	  // Set up file xml
  	  TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
        }

	// A referecnce for the object
	MapObject<K,V>& obj = *(TheNamedObjMap::Instance().getData< Handle<MapObject<K,V> > >(buffer_id));

	// The number of records should be recorded in the header
	int num_records;
	if(file_xml.count("/PropColorVectors")!=0)
	  read(file_xml, "/PropColorVectors/num_records", num_records);
	else if (file_xml.count("/GridProp")!=0)
	  read(file_xml, "/GridProp/num_records", num_records);
	else if (file_xml.count("/BlockProp")!=0)
	  read(file_xml, "/BlockProp/num_records", num_records);
	else{
	  QDPIO::cerr << __func__ << ": Can't find num_records "<< endl;
	  QDP_abort(1);
	}

	obj.openWrite(); // Prepare Object for insertion

	// Use the iterators to run through the object, reading each record
	for(int i=0; i < num_records; ++i)
	{
	  XMLReader local_record_xml;

	  // The key and value
	  K key;
	  V val;

	  read(to, local_record_xml, val);

	  // Extract the key
	  read(local_record_xml, "/MapEntry", key);

	  // Insert the key and value into the map
	  obj.insert(key, val);
	}

	obj.openRead(); // Done writing into object - switch to READ mode

	// Sanity check
	if (obj.size() != num_records)
	{
	  QDPIO::cerr << __func__ << ": Error: number of object records unexpected:"
		      << "  num_records= " << num_records
		      << "  obj.size= " << obj.size()
		      << endl;
	  QDP_abort(1);
	}

	// Setup a record xml
	XMLBufferWriter record_xml;
	if(file_xml.count("/PropColorVectors")!=0)
	  push(record_xml, "PropColorVector");
        else if (file_xml.count("/GridProp")!=0)
	  push(record_xml, "GridProp");
        else if (file_xml.count("/BlockProp")!=0)
	  push(record_xml, "BlockProp");
        else{
	  push(record_xml, "NamedObjectMap");
        }


	write(record_xml, "num_records", obj.size());  // This should be the size of num_records
	pop(record_xml);

	// Set File and Record XML throw away dummy XMLs
	TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);

	// Close and bolt
	close(to);
      }

      //! Local registration flag
      bool registered = false;

    }  // end namespace


    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("LatticePropagator"), 
								     QIOReadLatProp);
	success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("LatticePropagatorF"), 
								     QIOReadLatPropF);
	success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("LatticePropagatorD"), 
								     QIOReadLatPropD);

	success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("LatticeFermion"), 
								   QIOReadLatFerm);
//      success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("LatticeFermionF"), 
//								   QIOReadLatFermF);
//      success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("LatticeFermionD"), 
//								   QIOReadLatFermD);

	success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("Multi1dLatticeColorMatrix"), 
								     QIOReadArrayLatColMat);

	success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("Multi1dLatticeColorMatrixF"), 
								     QIOReadArrayLatColMatF);

	success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("Multi1dLatticeColorMatrixD"), 
								     QIOReadArrayLatColMatD);

	success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("QQDiquarkContract"), 
								     QIOReadQQDiquarkContract);
	
	success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("EigenInfoLatticeFermion"),
								     QIOReadEigenInfoLatticeFermion);

	success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("SubsetVectorsLatticeColorVector"),
								     QIOReadSubsetVectors<LatticeColorVector>);

	success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("RitzPairsLatticeFermion"), 
								     QIOReadRitzPairsLatticeFermion);
	
	success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("MapObjMemoryKeyPropColorVecLatticeFermion"), 
								     QIOReadMapObjMemory<KeyPropColorVec_t,LatticeFermion>);

	success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("MapObjMemoryKeyGridPropLatticeFermion"), QIOReadMapObjMemory<KeyGridProp_t,LatticeFermion>);

	success &= TheQIOReadObjFuncMap::Instance().registerFunction(string("MapObjMemoryKeyBlockPropLatticeFermion"), QIOReadMapObjMemory<KeyBlockProp_t,LatticeFermion>);

	registered = true;
      }
      return success;
    }
  }

}
