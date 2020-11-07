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

#include "qdp_map_obj_memory.h"

namespace Chroma 
{ 
  using namespace QDP;

  //! Object buffer
  void write(XMLWriter& xml, const std::string& path, const InlineGaugeTransfNamedObjEnv::Params::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_rot_id", input.gauge_rot_id);
    write(xml, "input_id", input.input_id);
    write(xml, "output_id", input.output_id);
    write(xml, "object_type", input.object_type);

    pop(xml);
  }


  //! Object buffer
  void read(XMLReader& xml, const std::string& path, InlineGaugeTransfNamedObjEnv::Params::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_rot_id", input.gauge_rot_id);
    read(inputtop, "input_id", input.input_id);
    read(inputtop, "output_id", input.output_id);
    read(inputtop, "object_type", input.object_type);
  }


  //! IO function std::map environment
  /*! \ingroup inlinehadron */
  namespace GaugeTransfObjCallMapEnv
  { 
    // Anonymous namespace
    namespace
    {
      struct DumbDisambiguator {};

      //! GaugeTransf function std::map
      /*! \ingroup inlinehadron */
      typedef SingletonHolder< 
	FunctionMap<DumbDisambiguator,
		    void,
		    std::string,
		    TYPELIST_3(const std::string&, const LatticeColorMatrix&, const std::string&),
		    void (*)(const std::string& output_id, const LatticeColorMatrix& g, const std::string& input_id),
		    StringFunctionMapError> >
      TheGaugeTransfObjFuncMap;


      //! Transform a generic object. This works only for non-array objects.
      template<typename T>
      void gaugeTransfObj(const std::string& output_id, const LatticeColorMatrix& g, const std::string& input_id)
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

      //! Transform a generic object. This works only for non-array objects.
       template<typename T>
       void gaugeTransfGaugeField(const std::string& output_id, const LatticeColorMatrix& g, const std::string& input_id)
       {

    	   T g_prec = g;
    	   // Grab the source
    	   const multi1d<T>& input_obj =
    			   TheNamedObjMap::Instance().getData< multi1d<T> >(input_id);

    	   XMLReader input_file_xml, input_record_xml;
    	   TheNamedObjMap::Instance().get(input_id).getFileXML(input_file_xml);
    	   TheNamedObjMap::Instance().get(input_id).getRecordXML(input_record_xml);

    	   // Create space for the target
    	   TheNamedObjMap::Instance().create< multi1d<T> >(output_id);
    	   TheNamedObjMap::Instance().get(output_id).setFileXML(input_file_xml);
    	   TheNamedObjMap::Instance().get(output_id).setRecordXML(input_record_xml);

    	   // Do the actual rotation
    	   multi1d<T>& dest_gauge = TheNamedObjMap::Instance().getData< multi1d<T> >(output_id);
    	   dest_gauge.resize(Nd);
    	   for(int mu=0; mu < Nd; ++mu) {
    		   //  r_mu = g * u * g^\dagger(x + mu)
    		   dest_gauge[mu] = g_prec*input_obj[mu]*adj(shift(g_prec,FORWARD,mu));
    	   }
       }
      //! Transform a subset_vectors object. This only works for non-array objects
      void gaugeTransfSubsetVectors(const std::string& output_id, const LatticeColorMatrix& g, const std::string& input_id)
      {
	// A shorthand for the object
	Handle< MapObject<int,EVPair<LatticeColorVector> > > input_obj =
	  TheNamedObjMap::Instance().getData< Handle< MapObject<int,EVPair<LatticeColorVector> > > >(input_id);

	XMLReader input_file_xml, input_record_xml;
	TheNamedObjMap::Instance().get(input_id).getFileXML(input_file_xml);
	TheNamedObjMap::Instance().get(input_id).getRecordXML(input_record_xml);

	// Create space for the target
	Handle< QDP::MapObject<int,EVPair<LatticeColorVector> > > output_obj(new MapObjectMemory<int,EVPair<LatticeColorVector> >());

	TheNamedObjMap::Instance().create< Handle< QDP::MapObject<int,EVPair<LatticeColorVector> > >, Handle< QDP::MapObject<int,EVPair<LatticeColorVector> > > >(output_id, output_obj);
	TheNamedObjMap::Instance().get(output_id).setFileXML(input_file_xml);
	TheNamedObjMap::Instance().get(output_id).setRecordXML(input_record_xml);

	 
	// Do the actual rotation. The evs stay the same
	for(int n=0; n < input_obj->size(); n++) {
	  EVPair<LatticeColorVector> pair;
	  input_obj->get(n,pair);

	  EVPair<LatticeColorVector> pair2 = pair; 
	  pair2.eigenVector  = g*pair.eigenVector;
	  output_obj->insert(n, pair2);
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
	success &= TheGaugeTransfObjFuncMap::Instance().registerFunction(std::string("LatticePropagator"), 
									 gaugeTransfObj<LatticePropagator>);
	success &= TheGaugeTransfObjFuncMap::Instance().registerFunction(std::string("LatticeFermion"), 
									 gaugeTransfObj<LatticeFermion>);
	success &= TheGaugeTransfObjFuncMap::Instance().registerFunction(std::string("LatticeStaggeredPropagator"), 
									 gaugeTransfObj<LatticeStaggeredPropagator>);
	success &= TheGaugeTransfObjFuncMap::Instance().registerFunction(std::string("LatticeStaggeredFermion"), 
									 gaugeTransfObj<LatticeStaggeredFermion>);
	success &= TheGaugeTransfObjFuncMap::Instance().registerFunction(std::string("SubsetVectorsLatticeColorVector"), 
									 gaugeTransfSubsetVectors);
	success &= TheGaugeTransfObjFuncMap::Instance().registerFunction(std::string("Multi1dLatticeColorMatrix"),
			gaugeTransfGaugeField<LatticeColorMatrix>);

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
	QDPIO::cerr << __func__ << ": caught Exception reading XML: " << e << std::endl;
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
		  << params.named_obj.object_type << std::endl;

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
		    << std::endl;
	QDP_abort(1);
      }
      catch (const std::string& e) 
      {
	QDPIO::cerr << name << ": error message: " << e 
		    << std::endl;
	QDP_abort(1);
      }
    
      QDPIO::cout << name << ": ran successfully" << std::endl;

      pop(xml_out);  // gaussian_init_named_obj

      END_CODE();
    } 

  } // namespace InlineGaugeTransfNamedObjEnv 

} // namespace Chroma
