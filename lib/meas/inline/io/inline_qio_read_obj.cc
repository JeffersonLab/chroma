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

	    // Determine if number of eigenvectors to read is specified in input xml
	    int num_vecs = 0;
	    std::istringstream xml_nam(params.named_obj_xml.xml);
	    XMLReader param_xml_top(xml_nam);
	    XMLReader param_xml(param_xml_top, "MapObject");
	    if(param_xml.count("num_vecs") > 0)
	      read(param_xml, "num_vecs", num_vecs);
	    else
	      num_vecs = N;
	    if(num_vecs > N)
	      {
		QDPIO::cerr<< "Error: number of vectors to read, num_vecs= " << num_vecs << ", is greater than number in file, N= " << N << endl;
		QDP_abort(1);
	      }

	    // Loop and read evecs
	    QDPIO::cout << "About to read " << num_vecs << " out of " << N << " evectors from QIO..."<<endl;
	    for(int n=0; n < num_vecs; n++)
	    {
	      XMLReader record_xml_dummy;
	      EVPair<LatticeColorVector> read_pair;
	  
	      read(to, record_xml_dummy, read_pair.eigenVector);
	  
	      read_pair.eigenValue.weights.resize(Lt);
	      read(record_xml_dummy, "/VectorInfo/weights", read_pair.eigenValue.weights);
	      QDPIO::cout << "Inserting Pair " << n << " into Map" << endl;
	      obj->insert(n, read_pair);

	    }
	    obj->flush();

	    // Set File and Record XML throw away dummy XMLs
	    // BJOO: Yes that's all very well. but RecordXML has to hold something
	    XMLBufferWriter dummy;
	    push(dummy,"DummyRecordXML");
	    write(dummy, "mapSize", num_vecs);
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
	//------------------------------------------------------------------------
	//! Read a fermion
	class QIOReadLatFerm : public QIOReadObject
	{
	private:
	  Params params;

	public:
	  QIOReadLatFerm(const Params& p) : params(p) {}

	  //! Read a propagator
	  void operator()(QDP_serialparallel_t serpar) {
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
	};

	// Call back
	QIOReadObject* qioReadLatFerm(const Params& p)
	{
	  return new QIOReadLatFerm(p);
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
	class QIOReadArrayLatColMat : public QIOReadObject
	{
	private:
	  Params params;

	public:
	  QIOReadArrayLatColMat(const Params& p) : params(p) {}

	  //! Read a propagator
	  void operator()(QDP_serialparallel_t serpar) {
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
	};

	// Call back
	QIOReadObject* qioReadArrayLatColMat(const Params& p)
	{
	  return new QIOReadArrayLatColMat(p);
	}


	//------------------------------------------------------------------------
	//! Read a single prec Gauge fields
	class QIOReadArrayLatColMatF : public QIOReadObject
	{
	private:
	  Params params;

	public:
	  QIOReadArrayLatColMatF(const Params& p) : params(p) {}

	  //! Read a propagator
	  void operator()(QDP_serialparallel_t serpar) {
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
	};

	// Call back
	QIOReadObject* qioReadArrayLatColMatF(const Params& p)
	{
	  return new QIOReadArrayLatColMatF(p);
	}


	//------------------------------------------------------------------------
	//! Read a Double Prec Gauge Field
	class QIOReadArrayLatColMatD : public QIOReadObject
	{
	private:
	  Params params;

	public:
	  QIOReadArrayLatColMatD(const Params& p) : params(p) {}

	  //! Read a propagator
	  void operator()(QDP_serialparallel_t serpar) {
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
	};


	// Call back
	QIOReadObject* qioReadArrayLatColMatD(const Params& p)
	{
	  return new QIOReadArrayLatColMatD(p);
	}


	//------------------------------------------------------------------------
	//! Read a QQDiquark object in floating precision
	class QIOReadQQDiquarkContract : public QIOReadObject
	{
	private:
	  Params params;

	public:
	  QIOReadQQDiquarkContract(const Params& p) : params(p) {}

	  //! Read a propagator
	  void operator()(QDP_serialparallel_t serpar) {
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
	};

	// Call back
	QIOReadObject* qioReadQQDiquarkContract(const Params& p)
	{
	  return new QIOReadQQDiquarkContract(p);
	}


	//------------------------------------------------------------------------
	class QIOReadEigenInfoLatticeFermion : public QIOReadObject
	{
	private:
	  Params params;

	public:
	  QIOReadEigenInfoLatticeFermion(const Params& p) : params(p) {}

	  //! Read a propagator
	  void operator()(QDP_serialparallel_t serpar) {
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
	};

	// Call back
	QIOReadObject* qioReadEigenInfoLatticeFermion(const Params& p)
	{
	  return new QIOReadEigenInfoLatticeFermion(p);
	}


	//-----------------------------------------------------------------------
	//! Read a RitzPairs Type
	class QIOReadRitzPairsLatticeFermion : public QIOReadObject
	{
	private:
	  Params params;

	public:
	  QIOReadRitzPairsLatticeFermion(const Params& p) : params(p) {}

	  //! Read a propagator
	  void operator()(QDP_serialparallel_t serpar) {
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
	};

	// Call back
	QIOReadObject* qioReadRitzPairsLatticeFermion(const Params& p)
	{
	  return new QIOReadRitzPairsLatticeFermion(p);
	}



	//------------------------------------------------------------------------
	//------------------------------------------------------------------------

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
      if( inputtop.count("parallel_io") == 1 ) { 
	read(inputtop, "parallel_io", input.parallel_io);
      }
      else { 
	input.parallel_io = false; 
      }
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
	QDP_serialparallel_t serpar;
	if ( params.file.parallel_io ) { 
	  QDPIO::cout << "Attempting Parallel IO read" << endl;
	  serpar = QDPIO_PARALLEL;
	}
	else { 
	  serpar = QDPIO_SERIAL;
	}

	// Read the object
	swatch.start();

	(*qioReadObject)(serpar);

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
