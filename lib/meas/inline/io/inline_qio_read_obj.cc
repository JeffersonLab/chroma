/*! \file
 * \brief Inline task to read an object from a named buffer
 *
 * Named object writing
 */

#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/inline_qio_read_obj.h"
#include "meas/inline/io/named_objmap.h"

#include "util/ferm/map_obj/map_obj_factory_w.h"
#include "util/ferm/map_obj/map_obj_aggregate_w.h"

#include "meas/hadron/diquark_w.h"
#include "util/ferm/eigeninfo.h"
#include "util/ferm/subset_vectors.h"
#include "util/ferm/key_prop_colorvec.h"
#include "handle.h"
#include "actions/ferm/invert/containers.h"

namespace Chroma 
{ 
  namespace InlineQIOReadNamedObjEnv 
  { 
    //! IO function map environment
    /*! \ingroup inlineio */
    namespace QIOReadObjectEnv
    { 
      /*! @ingroup inlineio */
      class QIOReadObject
      {
      public:
	~QIOReadObject() {}

	//! Read the object
	virtual void operator()(QDP_serialparallel_t serpar) = 0;
      };


      namespace
      {
	//! QIO read factory (foundry)
	/*! \ingroup inlineio */
	typedef SingletonHolder< 
	  ObjectFactory<QIOReadObject, 
			std::string,
			TYPELIST_1(const Params&),
			QIOReadObject* (*)(const Params&), 
			StringFactoryError> >
	TheQIOReadObjectFactory;
      }

      // Anonymous namespace
      namespace
      {
	//------------------------------------------------------------------------
	class QIOReadLatProp : public QIOReadObject
	{
	private:
	  Params params;

	public:
	  QIOReadLatProp(const Params& p) : params(p) {}

	  //! Read a propagator
	  void operator()(QDP_serialparallel_t serpar) {
	    LatticePropagator obj;
	    XMLReader file_xml, record_xml;

	    QDPFileReader to(file_xml,params.file.file_name,serpar);
	    read(to,record_xml,obj);
	    close(to);

	    TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.object_id);
	    TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.object_id) = obj;
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);
	  }
	};

	// Call back
	QIOReadObject* qioReadLatProp(const Params& p)
	{
	  return new QIOReadLatProp(p);
	}


	//------------------------------------------------------------------------
	//! Read a single prec propagator
	class QIOReadLatPropF : public QIOReadObject
	{
	private:
	  Params params;

	public:
	  QIOReadLatPropF(const Params& p) : params(p) {}

	  //! Read a propagator
	  void operator()(QDP_serialparallel_t serpar) {
	    LatticePropagatorF obj;
	    XMLReader file_xml, record_xml;

	    QDPFileReader to(file_xml,params.file.file_name,serpar);
	    read(to,record_xml,obj);
	    close(to);

	    TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.object_id);
	    TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.object_id) = obj;
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);
	  }
	};

	// Call back
	QIOReadObject* qioReadLatPropF(const Params& p)
	{
	  return new QIOReadLatPropF(p);
	}


	//------------------------------------------------------------------------
	//! Read a double prec propagator
	class QIOReadLatPropD : public QIOReadObject
	{
	private:
	  Params params;

	public:
	  QIOReadLatPropD(const Params& p) : params(p) {}

	  //! Read a propagator
	  void operator()(QDP_serialparallel_t serpar) {
	    LatticePropagatorD obj;
	    XMLReader file_xml, record_xml;

	    QDPFileReader to(file_xml,params.file.file_name,serpar);
	    read(to,record_xml,obj);
	    close(to);

	    TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.object_id);
	    TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.object_id) = obj;
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);
	  }
	};

	// Call back
	QIOReadObject* qioReadLatPropD(const Params& p)
	{
	  return new QIOReadLatPropD(p);
	}


	//------------------------------------------------------------------------
	//------------------------------------------------------------------------
	class QIOReadSubsetVectorsLCV : public QIOReadObject
	{
	private:
	  Params params;
	  Handle< MapObject<int,EVPair<LatticeColorVector> > > obj;

	public:
	  QIOReadSubsetVectorsLCV(const Params& params_) : params(params_) {
	    try
	    {
	      std::istringstream  xml_nam(params.named_obj_xml.xml);
	      XMLReader  namtop(xml_nam);
	      GroupXML_t colorvec_obj = readXMLGroup(namtop, "MapObject", "MapObjType");

	      // Create the entry
	      TheNamedObjMap::Instance().create< Handle< MapObject<int,EVPair<LatticeColorVector> > > >(params.named_obj.object_id);
	      TheNamedObjMap::Instance().getData< Handle< MapObject<int,EVPair<LatticeColorVector> > > >(params.named_obj.object_id) =
		TheMapObjIntKeyColorEigenVecFactory::Instance().createObject(colorvec_obj.id,
									     namtop,
									     colorvec_obj.path);
	    }
	    catch (std::bad_cast)
	    {
	      QDPIO::cerr << "QIOReadSubsetVectorsLCV: caught dynamic cast error" << endl;
	      QDP_abort(1);
	    }
	    catch (const string& e) 
	    {
	      QDPIO::cerr << "QIOReadSubsetVectorsLCV: error creating prop: " << e << endl;
	      QDP_abort(1);
	    }

	    // Cast should be valid now
	    obj = TheNamedObjMap::Instance().getData< Handle< MapObject<int,EVPair<LatticeColorVector> > > >(params.named_obj.object_id);
	  }

	  //! Read a propagator
	  void operator()(QDP_serialparallel_t serpar) {
	    // Stuff 
	    XMLReader file_xml, record_xml, largest_record_xml;

	    // Open file
	    QDPFileReader to(file_xml,params.file.file_name,serpar);

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

	    // Loop and read evecs
	    QDPIO::cout << "About to read " << N << " evectors from QIO..."<<endl;
	    obj->openWrite();
	    for(int n=0; n < N; n++)
	    {
	      XMLReader record_xml_dummy;
	      EVPair<LatticeColorVector> read_pair;
	  
	      read(to, record_xml_dummy, read_pair.eigenVector);
	  
	      read_pair.eigenValue.weights.resize(Lt);
	      read(record_xml_dummy, "/VectorInfo/weights", read_pair.eigenValue.weights);
	      QDPIO::cout << "Inserting Pair " << n << " into Map" << endl;
	      obj->insert(n, read_pair);

	    }
	    obj->openRead();

	    // Set File and Record XML throw away dummy XMLs
	    // BJOO: Yes that's all very well. but RecordXML has to hold something
	    XMLBufferWriter dummy;
	    push(dummy,"DummyRecordXML");
	    write(dummy, "mapSize", N);
	    pop(dummy);
	    record_xml.open(dummy);

	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);
	    
	    // Done
	    close(to);
	  }
	};

	// Call back
	QIOReadObject* qioReadSubsetVectorsLCV(const Params& p)
	{
	  return new QIOReadSubsetVectorsLCV(p);
	}



	//------------------------------------------------------------------------
	//------------------------------------------------------------------------
	//------------------------------------------------------------------------
	//------------------------------------------------------------------------
	//------------------------------------------------------------------------
	//------------------------------------------------------------------------
#if 0
	//------------------------------------------------------------------------
	//! Read a fermion
	void QIOReadLatFerm(const string& buffer_id,
			    const string& file, 
			    QDP_serialparallel_t serpar)
	{
	  LatticeFermion obj;
	  XMLReader file_xml, record_xml;

	  QDPFileReader to(file_xml,params.file.file_name,serpar);
	  read(to,record_xml,obj);
	  close(to);

	  TheNamedObjMap::Instance().create<LatticeFermion>(params.named_obj.object_id);
	  TheNamedObjMap::Instance().getData<LatticeFermion>(params.named_obj.object_id) = obj;
	  TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	  TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);
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

	  QDPFileReader to(file_xml,params.file.file_name,serpar);
	  read(to,record_xml,obj);
	  close(to);

	  TheNamedObjMap::Instance().create<LatticeFermion>(params.named_obj.object_id);
	  TheNamedObjMap::Instance().getData<LatticeFermion>(params.named_obj.object_id) = obj;
	  TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	  TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);
	}


	//! Read a double prec fermion
	void QIOReadLatFermD(const string& buffer_id,
			     const string& file, 
			     QDP_serialparallel_t serpar)
	{
	  LatticeFermionD obj;
	  XMLReader file_xml, record_xml;

	  QDPFileReader to(file_xml,params.file.file_name,serpar);
	  read(to,record_xml,obj);
	  close(to);

	  TheNamedObjMap::Instance().create<LatticeFermion>(params.named_obj.object_id);
	  TheNamedObjMap::Instance().getData<LatticeFermion>(params.named_obj.object_id) = obj;
	  TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	  TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);
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

	  QDPFileReader to(file_xml,params.file.file_name,serpar);
	  read(to,record_xml,obj);
	  close(to);

	  TheNamedObjMap::Instance().create< multi1d<LatticeColorMatrix> >(params.named_obj.object_id);
	  TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.object_id) = obj;
	  TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	  TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);
	}

	//! Read a single prec Gauge fields
	void QIOReadArrayLatColMatF(const string& buffer_id,
				    const string& file, 
				    QDP_serialparallel_t serpar)
	{
	  multi1d<LatticeColorMatrixF> obj;
	  XMLReader file_xml, record_xml;

	  obj.resize(Nd);    // BAD BAD BAD - FIX THIS

	  QDPFileReader to(file_xml,params.file.file_name,serpar);
	  read(to,record_xml,obj);
	  close(to);

	  TheNamedObjMap::Instance().create< multi1d<LatticeColorMatrix> >(params.named_obj.object_id);
	  (TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.object_id)).resize(Nd);
	  for(int i=0; i < Nd; i++) {
	    (TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.object_id))[i] = obj[i];
	  }
	  TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	  TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);
	}

	//! Read a Double Prec Gauge Field
	void QIOReadArrayLatColMatD(const string& buffer_id,
				    const string& file, 
				    QDP_serialparallel_t serpar)
	{
	  multi1d<LatticeColorMatrixD> obj;
	  XMLReader file_xml, record_xml;

	  obj.resize(Nd);    // BAD BAD BAD - FIX THIS

	  QDPFileReader to(file_xml,params.file.file_name,serpar);
	  read(to,record_xml,obj);
	  close(to);

	  TheNamedObjMap::Instance().create< multi1d<LatticeColorMatrix> >(params.named_obj.object_id);
	  (TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.object_id)).resize(Nd);
	  for(int i=0; i < Nd;i++) {
	    (TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.object_id))[i] = obj[i];
	  }

	  TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	  TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);
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

	  QDPFileReader to(file_xml,params.file.file_name,serpar);
	  read(to,record_xml,obj_1d);
	  close(to);

	  // Create the named object
	  TheNamedObjMap::Instance().create<QQDiquarkContract_t>(params.named_obj.object_id);
	  TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	  TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);

	  // Convert the 1-d object into an N-d object
	  QQDiquarkContract_t& obj = TheNamedObjMap::Instance().getData<QQDiquarkContract_t>(params.named_obj.object_id);
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
	  QDPFileReader to(file_xml,params.file.file_name,serpar);

	  // Create object
	  TheNamedObjMap::Instance().create< EigenInfo<LatticeFermion> >(params.named_obj.object_id);

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
	  TheNamedObjMap::Instance().getData< EigenInfo<LatticeFermion> >(params.named_obj.object_id).getLargest()=largestD[0];

	  // REsize eval arrays so that IO works correctly
	  evalsD.resize(size);

	  // Read evals
	  read(to, record_xml, evalsD);

	  QDPIO::cout << "Read " << evalsD.size() << "Evals " << endl;
	  for(int i=0; i < evalsD.size(); i++) { 
	    QDPIO::cout << "Eval["<<i<<"] = " << Real(evalsD[i]) << endl;
	  }

	  // Resize eval array 
	  TheNamedObjMap::Instance().getData< EigenInfo<LatticeFermion> >(params.named_obj.object_id).getEvalues().resize(evalsD.size());
      
	  // Downcast evals to Real() precision
	  for (int i=0; i<evalsD.size(); i++)
	    TheNamedObjMap::Instance().getData< EigenInfo<LatticeFermion> >(params.named_obj.object_id).getEvalues()[i] = Real(evalsD[i]);


	  // Resize evec arrays
	  // THIS IS AN EIGEN INFO - NOT SUBSET VECTORS
	  TheNamedObjMap::Instance().getData< EigenInfo<LatticeFermion> >(params.named_obj.object_id).getEvectors().resize(evalsD.size());

	  // Loop and read evecs
	  for (int i=0; i<evalsD.size(); i++)
	  {
	    XMLReader record_xml_dummy;
	    read(to, record_xml_dummy, evecD);
	    LatticeFermion evecR;
	    evecR=evecD;

	    TheNamedObjMap::Instance().getData< EigenInfo<LatticeFermion> >(params.named_obj.object_id).getEvectors()[i]=evecR;
	  }
 
	  // Set File and Record XML throw away dummy XMLs
	  TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	  TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);

	  // Done - That too was unnecessarily painful
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
	  QDPFileReader to(file_xml,params.file.file_name,serpar);

	  // Create the named object
	  TheNamedObjMap::Instance().create< LinAlg::RitzPairs<LatticeFermion> >(params.named_obj.object_id);
	  TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);

	  // A shorthand for the object
	  LinAlg::RitzPairs<LatticeFermion>& obj = 
	    TheNamedObjMap::Instance().getData<LinAlg::RitzPairs< LatticeFermion> >(params.named_obj.object_id);

	  XMLReader record_xml;
	  TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);
	
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
	  QDPFileReader to(file_xml,params.file.file_name,serpar);

	  // Create object if it does not exist. Otherwise, we can continue reading into it.
	  if (! TheNamedObjMap::Instance().check(params.named_obj.object_id))
	  {
	    // Make a memory map object -- later make with factory?
	    Handle< MapObject<K,V> > obj_map_handle( new MapObjectMemory<K,V> );
 
	    // Create slot
	    TheNamedObjMap::Instance().create< Handle<MapObject<K,V> > >(params.named_obj.object_id);
	    // Insert Handle
	    TheNamedObjMap::Instance().getData< Handle<MapObject<K,V> > >(params.named_obj.object_id) = obj_map_handle;

	    // Set up file xml
	    TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	  }

	  // A referecnce for the object
	  MapObject<K,V>& obj = *(TheNamedObjMap::Instance().getData< Handle<MapObject<K,V> > >(params.named_obj.object_id));

	  // The number of records should be recorded in the header
	  int num_records;
	  if(file_xml.count("/PropColorVectors")!=0)
	    read(file_xml, "/PropColorVectors/num_records", num_records);
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
	  else{
	    push(record_xml, "NamedObjectMap");
	  }


	  write(record_xml, "num_records", obj.size());  // This should be the size of num_records
	  pop(record_xml);

	  // Set File and Record XML throw away dummy XMLs
	  TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	  TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);

	  // Close and bolt
	  close(to);
	}
	//------------------------------------------------------------------------
	//------------------------------------------------------------------------
	//------------------------------------------------------------------------
	//------------------------------------------------------------------------
	//------------------------------------------------------------------------
	//------------------------------------------------------------------------
	//------------------------------------------------------------------------
#endif

	//! Local registration flag
	bool registered = false;

      }  // end namespace


      //! Register all the factories
      bool registerAll() 
      {
	bool success = true; 
	if (! registered)
	{
	  success &= TheQIOReadObjectFactory::Instance().registerObject(string("LatticePropagator"),   
									qioReadLatProp);
	  success &= TheQIOReadObjectFactory::Instance().registerObject(string("LatticePropagatorF"), 
									qioReadLatPropF);
	  success &= TheQIOReadObjectFactory::Instance().registerObject(string("LatticePropagatorD"), 
									qioReadLatPropD);

	  success &= TheQIOReadObjectFactory::Instance().registerObject(string("SubsetVectorsLatticeColorVector"),
									qioReadSubsetVectorsLCV);

#if 0
	  success &= TheQIOReadObjectFactory::Instance().registerObject(string("MapObjMemoryKeyPropColorVecLatticeFermion"), 
									qioReadMapObjMemory<KeyPropColorVec_t,LatticeFermion>);

	  success &= TheQIOReadObjectFactory::Instance().registerObject(string("LatticeFermion"), 
									qioReadLatFerm);
//      success &= TheQIOReadObjectFactory::Instance().registerObject(string("LatticeFermionF"), 
//								   qioReadLatFermF);
//      success &= TheQIOReadObjectFactory::Instance().registerObject(string("LatticeFermionD"), 
//								   qioReadLatFermD);

	  success &= TheQIOReadObjectFactory::Instance().registerObject(string("Multi1dLatticeColorMatrix"), 
									qioReadArrayLatColMat);

	  success &= TheQIOReadObjectFactory::Instance().registerObject(string("Multi1dLatticeColorMatrixF"), 
									qioReadArrayLatColMatF);

	  success &= TheQIOReadObjectFactory::Instance().registerObject(string("Multi1dLatticeColorMatrixD"), 
									qioReadArrayLatColMatD);

	  success &= TheQIOReadObjectFactory::Instance().registerObject(string("QQDiquarkContract"), 
									qioReadQQDiquarkContract);
	
	  success &= TheQIOReadObjectFactory::Instance().registerObject(string("EigenInfoLatticeFermion"),
									qioReadEigenInfoLatticeFermion);

	  success &= TheQIOReadObjectFactory::Instance().registerObject(string("RitzPairsLatticeFermion"), 
									qioReadRitzPairsLatticeFermion);
	
#endif

	  registered = true;
	}
	return success;
      }
    } // namespace QIOReadObjectEnv


    //------------------------------------------------------------------------
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMeas(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;

      const std::string name = "QIO_READ_NAMED_OBJECT";
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	success &= QIOReadObjectEnv::registerAll();
	success &= MapObjectWilson4DEnv::registerAll();
	registered = true;
      }
      return success;
    }


    //------------------------------------------------------------------------
    struct NamedObject_t
    {
      std::string   object_id;
      std::string   object_type;
    };


    //! Object buffer
    void read(XMLReader& xml, const string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "object_id", input.object_id);
      read(inputtop, "object_type", input.object_type);
    }

    //! File output
    void read(XMLReader& xml, const string& path, Params::File_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "file_name", input.file_name);
    }


    // Param stuff
    Params::Params() { frequency = 0; }

    Params::Params(XMLReader& xml_in, const std::string& path) 
    {
      try 
      {
	XMLReader paramtop(xml_in, path);

	if (paramtop.count("Frequency") == 1)
	  read(paramtop, "Frequency", frequency);
	else
	  frequency = 1;

	// Named objects - original copy
	named_obj_xml = readXMLGroup(paramtop, "NamedObject", "object_type");

	// Named objects - make a concrete copy. Useful for digging out basic stuff.
	read(paramtop, "NamedObject", named_obj);

	// File name
	read(paramtop, "File", file);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << __func__ << ": caught Exception reading XML: " << e << endl;
	QDP_abort(1);
      }
    }


    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      START_CODE();

      push(xml_out, "qio_read_named_obj");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": object reader" << endl;
      StopWatch swatch;

      // Read the object
      // ONLY SciDAC output format is supported in this task
      // Other tasks could support other disk formats
      try
      {	
	QDPIO::cout << "Attempt to read object name = " << params.named_obj.object_id << endl;
	write(xml_out, "object_id", params.named_obj.object_id);

	// Create the object reader
	Handle<QIOReadObjectEnv::QIOReadObject> qioReadObject(
	  QIOReadObjectEnv::TheQIOReadObjectFactory::Instance().createObject(params.named_obj.object_type,
									     params));

	// Read the object
	swatch.start();

	(*qioReadObject)(QDPIO_SERIAL);

	swatch.stop();

	QDPIO::cout << "Object successfully read: time= " 
		    << swatch.getTimeInSeconds() 
		    << " secs" << endl;
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << name << ": caught dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << name << ": map call failed: " << e 
		    << endl;
	QDP_abort(1);
      }
    
      QDPIO::cout << name << ": ran successfully" << endl;

      pop(xml_out);  // qio_read_named_obj

      END_CODE();
    } 

  }

}
