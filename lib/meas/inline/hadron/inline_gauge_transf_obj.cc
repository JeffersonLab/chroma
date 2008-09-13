// $Id: inline_gauge_transf_obj.cc,v 3.2 2008-09-13 21:39:53 edwards Exp $
/*! \file
 * \brief Inline task gauge transform some fermion object
 *
 * Gauge transform a named object
 */

#include "singleton.h"
#include "funcmap.h"
#include "chromabase.h"

#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/inline/hadron/inline_gauge_transf_obj.h"

#include "util/ferm/subset_vectors.h"

namespace Chroma 
{ 
  //! Object buffer
  void write(XMLWriter& xml, const string& path, const InlineGaugeTransfNamedObjEnv::Params::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_rot_id", input.gauge_rot_id);
    write(xml, "input_id", input.input_id);
    write(xml, "output_id", input.output_id);
    write(xml, "object_type", input.object_type);

    pop(xml);
  }


  //! Object buffer
  void read(XMLReader& xml, const string& path, InlineGaugeTransfNamedObjEnv::Params::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_rot_id", input.gauge_rot_id);
    read(inputtop, "input_id", input.input_id);
    read(inputtop, "output_id", input.output_id);
    read(inputtop, "object_type", input.object_type);
  }


  //! IO function map environment
  /*! \ingroup inlinehadron */
  namespace GaugeTransfObjCallMapEnv
  { 
    // Anonymous namespace
    namespace
    {
      struct DumbDisambiguator {};

      //! GaugeTransf function map
      /*! \ingroup inlinehadron */
      typedef SingletonHolder< 
	FunctionMap<DumbDisambiguator,
		    void,
		    std::string,
		    TYPELIST_3(const string&, const LatticeColorMatrix&, const string&),
		    void (*)(const string& output_id, const LatticeColorMatrix& g, const string& input_id),
		    StringFunctionMapError> >
      TheGaugeTransfObjFuncMap;


      //! Transform a generic object. This works only for non-array objects.
      template<typename T>
      void gaugeTransfObj(const string& output_id, const LatticeColorMatrix& g, const string& input_id)
      {
	// Grab the source
	const T& input_obj = 
	  TheNamedObjMap::Instance().getData<T>(input_id);

	XMLReader input_file_xml, input_record_xml;
	TheNamedObjMap::Instance().get(input_id).getFileXML(input_file_xml);
	TheNamedObjMap::Instance().get(input_id).getRecordXML(input_record_xml);

	// Create space for the target
	TheNamedObjMap::Instance().create<T>(output_id);
	TheNamedObjMap::Instance().get(output_id).setFileXML(input_file_xml);
	TheNamedObjMap::Instance().get(output_id).setRecordXML(input_record_xml);

	// Do the actual rotation
	TheNamedObjMap::Instance().getData<T>(output_id) = g * input_obj;
      }


      //! Transform a subset_vectors object. This only works for non-array objects
      template<typename T>
      void gaugeTransfSubsetVectors(const string& output_id, const LatticeColorMatrix& g, const string& input_id)
      {
	// A shorthand for the object
	SubsetVectors<T>& input_obj = TheNamedObjMap::Instance().getData< SubsetVectors<T> >(input_id);

	XMLReader input_file_xml, input_record_xml;
	TheNamedObjMap::Instance().get(input_id).getFileXML(input_file_xml);
	TheNamedObjMap::Instance().get(input_id).getRecordXML(input_record_xml);

	// Create space for the target
	TheNamedObjMap::Instance().create< SubsetVectors<T> >(output_id);
	SubsetVectors<T>& output_obj = TheNamedObjMap::Instance().getData< SubsetVectors<T> >(output_id);
	TheNamedObjMap::Instance().get(output_id).setFileXML(input_file_xml);
	TheNamedObjMap::Instance().get(output_id).setRecordXML(input_record_xml);

	// I'm lazy. Use the input object to initialize the output object.
	output_obj = input_obj;

	// Do the actual rotation. The evs stay the same
	for(int n=0; n < input_obj.getNumVectors(); n++)
	{
	  output_obj.getEvectors()[n] = g * input_obj.getEvectors()[n];
	}
      }


      //! Local registration flag
      bool registered = false;

    }  // end namespace GaugeTransfCallMap

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheGaugeTransfObjFuncMap::Instance().registerFunction(string("LatticePropagator"), 
									 gaugeTransfObj<LatticePropagator>);
	success &= TheGaugeTransfObjFuncMap::Instance().registerFunction(string("LatticeFermion"), 
									 gaugeTransfObj<LatticeFermion>);
	success &= TheGaugeTransfObjFuncMap::Instance().registerFunction(string("LatticeStaggeredPropagator"), 
									 gaugeTransfObj<LatticeStaggeredPropagator>);
	success &= TheGaugeTransfObjFuncMap::Instance().registerFunction(string("LatticeStaggeredFermion"), 
									 gaugeTransfObj<LatticeStaggeredFermion>);

	success &= TheGaugeTransfObjFuncMap::Instance().registerFunction(string("SubsetVectorsLatticeColorVector"), 
									 gaugeTransfSubsetVectors<LatticeColorVector>);
	registered = true;
      }
      return success;
    }
  }  // end CallMap namespace


  namespace InlineGaugeTransfNamedObjEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMeas(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "GAUGE_TRANSFORM_NAMED_OBJECT";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// Gaussian init functions
	success &= GaugeTransfObjCallMapEnv::registerAll();

	// Inline measurement
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
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
    Params::writeXML(XMLWriter& xml_out, const std::string& path) 
    {
      push(xml_out, path);
    
      // Ids
      write(xml_out, "NamedObject", named_obj);

      pop(xml_out);
    }


    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      START_CODE();

      push(xml_out, "gauge_transf_object");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": gauge transform an object of type "
		  << params.named_obj.object_type << endl;

      // Grab the input object
      try
      {
	const LatticeColorMatrix& g = 
	  TheNamedObjMap::Instance().getData<LatticeColorMatrix>(params.named_obj.gauge_rot_id);

	// Gaussian init the object
	GaugeTransfObjCallMapEnv::TheGaugeTransfObjFuncMap::Instance().callFunction(params.named_obj.object_type,
										    params.named_obj.output_id,
										    g,
										    params.named_obj.input_id);
      }
      catch (std::bad_cast) 
      {
	QDPIO::cerr << name << ": cast error" 
		    << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << name << ": error message: " << e 
		    << endl;
	QDP_abort(1);
      }
    
      QDPIO::cout << name << ": ran successfully" << endl;

      pop(xml_out);  // gaussian_init_named_obj

      END_CODE();
    } 

  } // namespace InlineGaugeTransfNamedObjEnv 

} // namespace Chroma
