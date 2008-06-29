// $Id: inline_gaussian_obj.cc,v 3.4 2008-06-29 20:17:08 edwards Exp $
/*! \file
 * \brief Inline task to gaussian init a named object
 *
 * Named object initialization
 */

#include "singleton.h"
#include "funcmap.h"
#include "chromabase.h"

#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/inline/io/inline_gaussian_obj.h"

#include "util/ferm/eigeninfo.h"

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

      //! Init a faky gaussian eigeninfo struct of T-s
      template<typename T>
      void gaussianInitEigenInfo(const string& buffer_id)
      {
	// A shorthand for the object
	TheNamedObjMap::Instance().create< EigenInfo<T> >(buffer_id);
	EigenInfo<T>& obj=TheNamedObjMap::Instance().getData< EigenInfo<T> >(buffer_id);

	// Put 4 in here just for fun
	int N = 4;
	obj.getEvalues().resize(N);
	obj.getEvectors().resize(N);

	for(int m=0; m < N; ++m)
	{
	  random(obj.getEvalues()[m]);
	  gaussian(obj.getEvectors()[m]);
	}

	//
	// Orthogonalize the vectors
	//
	// To get a time direction subset
	SftMom phases(0, false, Nd-1);

	// Do this the inefficient way - loop over all subsets
	for(int t=0; t < phases.numSubsets(); ++t)
	{
	  const Subset& s = phases.getSet()[t];

	  for(int m=0; m < N; ++m)
	  {
	    // Convenience
	    T& v = obj.getEvectors()[m];

	    if (m > 0)
	    {
	      // Orthogonalize this vector against the previous "m" of them
	      GramSchm(v, obj.getEvectors(), m, s);
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
	success &= TheGaussianInitObjFuncMap::Instance().registerFunction(string("EigenInfoLatticeColorVector"), 
									  gaussianInitEigenInfo<LatticeColorVector>);
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
