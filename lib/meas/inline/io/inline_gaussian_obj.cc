// $Id: inline_gaussian_obj.cc,v 3.7 2008-09-13 19:56:40 edwards Exp $
/*! \file
 * \brief Inline task to gaussian init a named object
 *
 * Named object initialization
 */

#include "singleton.h"
#include "funcmap.h"
#include "chromabase.h"
#include "handle.h"

#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/inline/io/inline_gaussian_obj.h"

#include "util/ferm/subset_vectors.h"
#include "util/ferm/key_prop_colorvec.h"
#include "util/ferm/map_obj.h"
#include "util/ferm/map_obj/map_obj_memory.h"

#include "util/ft/sftmom.h"
#include "meas/eig/gramschm.h"

namespace Chroma 
{ 
  //! IO function map environment
  /*! \ingroup inlineio */
  namespace GaussianInitObjCallMapEnv
  { 
    // Anonymous namespace
    namespace
    {
      struct DumbDisambiguator {};

      //! GaussianInit init function map
      /*! \ingroup inlineio */
      typedef SingletonHolder< 
	FunctionMap<DumbDisambiguator,
		    void,
		    std::string,
		    TYPELIST_1(const string&),
		    void (*)(const string& buffer_id),
		    StringFunctionMapError> >
      TheGaussianInitObjFuncMap;


      //! Init a propagator
      template<typename T>
      void gaussianInitObj(const string& buffer_id)
      {
	TheNamedObjMap::Instance().create<T>(buffer_id);
	XMLBufferWriter file_xml, record_xml;

	push(file_xml,"FileXML");
	pop(file_xml);

	push(record_xml,"RecordXML");
	pop(record_xml);

	gaussian(TheNamedObjMap::Instance().getData<T>(buffer_id));
	TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
      }

      //! Init a propagator
      void gaussianInitMulti1dLatColMat(const string& buffer_id)
      {
	multi1d<LatticeColorMatrix> u(Nd);
	for(int mu=0; mu < u.size(); ++mu)
	  gaussian(u[mu]);

	XMLBufferWriter file_xml, record_xml;

	push(file_xml,"FileXML");
	pop(file_xml);

	push(record_xml,"RecordXML");
	pop(record_xml);

	TheNamedObjMap::Instance().create< multi1d<LatticeColorMatrix> >(buffer_id);
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(buffer_id) = u;
	TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
      }

      //! Init a faky gaussian subset of vectors over T-s
      template<typename T>
      void gaussianInitSubsetVectors(const string& buffer_id)
      {
	// A shorthand for the object
	TheNamedObjMap::Instance().create< SubsetVectors<T> >(buffer_id);
	SubsetVectors<T>& obj = TheNamedObjMap::Instance().getData< SubsetVectors<T> >(buffer_id);

	// To get a time direction subset
	int decay_dir = Nd-1;
	SftMom phases(0, false, decay_dir);
	const int Lt = phases.numSubsets();

	// Put 4 in here just for fun
	int N = 4;
	obj.getEvectors().resize(N);
	obj.getEvalues().resize(N);
	obj.getDecayDir() = decay_dir;

 	for(int n=0; n < N; ++n)
	{
	  obj.getEvalues()[n].weights.resize(Lt);
	  gaussian(obj.getEvectors()[n]);
	}

	//
	// Orthogonalize the vectors
	//
	// Do this the inefficient way - loop over all subsets
	for(int t=0; t < Lt; ++t)
	{
	  const Subset& s = phases.getSet()[t];

	  for(int n=0; n < N; ++n)
	  {
	    random(obj.getEvalues()[n].weights[t]);

	    // Convenience
	    T& v = obj.getEvectors()[n];

	    if (n > 0)
	    {
	      // Orthogonalize this vector against the previous "n" of them
	      GramSchm(v, obj.getEvectors(), n, s);
	    }

	    // Normalize
	    v[s] /= sqrt(norm2(v, s));
	  }
	}

	// I haven't figure out what to put in here
	XMLBufferWriter file_xml, record_xml;

	push(file_xml,"FileXML");
	pop(file_xml);

	push(record_xml,"RecordXML");
	pop(record_xml);

	TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
      }


      //! Init a 3 unit vectors in an eigeninfo struct of T-s
      /*! Impose site level orthogonality */
      void unitInitSubsetVectorsLatColVec(const string& buffer_id)
      {
	// A shorthand for the object
	TheNamedObjMap::Instance().create< SubsetVectors<LatticeColorVector> >(buffer_id);
	SubsetVectors<LatticeColorVector>& obj=TheNamedObjMap::Instance().getData< SubsetVectors<LatticeColorVector> >(buffer_id);
	// Time extent
	int decay_dir = Nd-1;
	const int Lt = QDP::Layout::lattSize()[decay_dir];

	// Use Nc vectors. There are only Nc site-level orthog. vectors in SU(N)
	int N = Nc;
	obj.getEvectors().resize(N);
	obj.getEvalues().resize(N);
	obj.getDecayDir() = decay_dir;

	for(int n=0; n < N; ++n)
	{
	  ColorVector vec = zero;
	  pokeColor(vec, cmplx(Real(1),Real(0)), n);

	  obj.getEvectors()[n] = vec;   // Same for all sites
	  obj.getEvalues()[n].weights.resize(Lt);
	  obj.getEvalues()[n].weights = zero;  // There are all zero for the constant field
	}

	// I haven't figure out what to put in here
	XMLBufferWriter file_xml, record_xml;

	push(file_xml,"FileXML");
	pop(file_xml);

	push(record_xml,"RecordXML");
	pop(record_xml);

	TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
      }


      //! Init a faky gaussian MapObject Type of struct of keys and vals
      void gaussianInitMapObjKeyPropColorVecLatFerm(const string& buffer_id)
      {
	// Create the object
	{ 
	  // Make a new object - for now memory. Later create with factory 
	  Handle<MapObject<KeyPropColorVec_t,LatticeFermion> > new_obj_handle( new MapObjectMemory<KeyPropColorVec_t,LatticeFermion>() );

	  // Create slot
	  TheNamedObjMap::Instance().create< Handle<MapObject<KeyPropColorVec_t,LatticeFermion> > >(buffer_id);

	  // Insert
	  TheNamedObjMap::Instance().getData< Handle<MapObject<KeyPropColorVec_t,LatticeFermion> > >(buffer_id) = new_obj_handle;
	}
	
	// Now make a working reference
	MapObject<KeyPropColorVec_t,LatticeFermion>& obj = 
	  *(TheNamedObjMap::Instance().getData< Handle<MapObject<KeyPropColorVec_t,LatticeFermion> > >(buffer_id));

	// Fix these. Put two sources in here.
	int N = 8;
	multi1d<int> t_sources(2);
	t_sources[0] = 0;
	t_sources[1] = QDP::Layout::lattSize()[Nd-1] / 2;

	obj.openWrite(); // Prepare MapObject for writing into.

	// Loop over each source and create a fake set of propagators
	for(int tt=0; tt < t_sources.size(); ++tt)
	{
	  int t0 = t_sources[tt];

	  for(int colorvec_src=0; colorvec_src < N; ++colorvec_src)
	  {
	    for(int spin_src=0; spin_src < Ns; ++spin_src)
	    {
	      // The key
	      KeyPropColorVec_t key;
	      key.t_source     = t0;
	      key.colorvec_src = colorvec_src;
	      key.spin_src     = spin_src;

	      // The value
	      LatticeFermion val;
	      gaussian(val);

	      // Insert the key and value into the map
	      obj.insert(key, val);
	    } // spin_src
	  } // colorvec_src
	} // tt

	obj.openRead(); // Done inserting

	// Write the meta-data for this operator
	{
	  XMLBufferWriter file_xml;

	  push(file_xml, "PropColorVectors");
	  write(file_xml, "num_records", obj.size()); 
	  push(file_xml, "Params");
	  pop(file_xml);
	  push(file_xml, "Config_info");
	  pop(file_xml);
	  pop(file_xml); // PropColorVectors

	  XMLBufferWriter record_xml;
	  push(record_xml, "PropColorVector");
	  write(record_xml, "num_records", obj.size()); 
	  pop(record_xml);

	  // Write the propagator xml info
	  TheNamedObjMap::Instance().get(buffer_id).setFileXML(file_xml);
	  TheNamedObjMap::Instance().get(buffer_id).setRecordXML(record_xml);
	}
      }


      //! Local registration flag
      bool registered = false;

    }  // end namespace GaussianInitCallMap

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheGaussianInitObjFuncMap::Instance().registerFunction(string("LatticePropagator"), 
									  gaussianInitObj<LatticePropagator>);
	success &= TheGaussianInitObjFuncMap::Instance().registerFunction(string("LatticeFermion"), 
									  gaussianInitObj<LatticeFermion>);
	success &= TheGaussianInitObjFuncMap::Instance().registerFunction(string("LatticeStaggeredPropagator"), 
									  gaussianInitObj<LatticeStaggeredPropagator>);
	success &= TheGaussianInitObjFuncMap::Instance().registerFunction(string("LatticeStaggeredFermion"), 
									  gaussianInitObj<LatticeStaggeredFermion>);
	success &= TheGaussianInitObjFuncMap::Instance().registerFunction(string("SubsetVectorsLatticeColorVector"), 
									  gaussianInitSubsetVectors<LatticeColorVector>);
	success &= TheGaussianInitObjFuncMap::Instance().registerFunction(string("UnitSubsetVectorsLatticeColorVector"), 
									  unitInitSubsetVectorsLatColVec);
	success &= TheGaussianInitObjFuncMap::Instance().registerFunction(string("MapObjectKeyPropColorVecLatticeFermion"), 
									  gaussianInitMapObjKeyPropColorVecLatFerm);

	registered = true;
      }
      return success;
    }
  }  // end CallMap namespace


  namespace InlineGaussianInitNamedObjEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineGaussianInitNamedObj(InlineGaussianInitNamedObjParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "GAUSSIAN_INIT_NAMED_OBJECT";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// Gaussian init functions
	success &= GaussianInitObjCallMapEnv::registerAll();
	// Inline measurement
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
  }


  //! Object buffer
  void write(XMLWriter& xml, const string& path, const InlineGaussianInitNamedObjParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "object_id", input.object_id);
    write(xml, "object_type", input.object_type);

    pop(xml);
  }


  //! Object buffer
  void read(XMLReader& xml, const string& path, InlineGaussianInitNamedObjParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "object_id", input.object_id);
    read(inputtop, "object_type", input.object_type);
  }


  // Param stuff
  InlineGaussianInitNamedObjParams::InlineGaussianInitNamedObjParams() { frequency = 0; }

  InlineGaussianInitNamedObjParams::InlineGaussianInitNamedObjParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;

      // Ids
      read(paramtop, "NamedObject", named_obj);
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << __func__ << ": caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  void
  InlineGaussianInitNamedObjParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    // Ids
    Chroma::write(xml_out, "NamedObject", named_obj);

    pop(xml_out);
  }


  void 
  InlineGaussianInitNamedObj::operator()(unsigned long update_no,
					 XMLWriter& xml_out) 
  {
    START_CODE();

    push(xml_out, "gaussian_init_named_obj");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlineGaussianInitNamedObjEnv::name << ": gaussian init an object of type " 
		<< params.named_obj.object_type << endl;

    // Gaussian the object
    QDPIO::cout << "Attempt to list all object names" << endl;
    try
    {
      // Gaussian init the object
      GaussianInitObjCallMapEnv::TheGaussianInitObjFuncMap::Instance().callFunction(params.named_obj.object_type,
										    params.named_obj.object_id);
    }
    catch (std::bad_cast) 
    {
      QDPIO::cerr << InlineGaussianInitNamedObjEnv::name << ": cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineGaussianInitNamedObjEnv::name << ": error message: " << e 
		  << endl;
      QDP_abort(1);
    }
    
    QDPIO::cout << InlineGaussianInitNamedObjEnv::name << ": ran successfully" << endl;

    pop(xml_out);  // gaussian_init_named_obj

    END_CODE();
  } 

}
