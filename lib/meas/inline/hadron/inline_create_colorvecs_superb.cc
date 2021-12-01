/*! \file
 * \brief Construct colorvectors of the laplacian with PRIMME
 */

#include "fermact.h"
#include "meas/inline/hadron/inline_create_colorvecs_superb.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "util/ferm/superb_contractions.h"
#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/gaus_smear.h"
#include "meas/glue/mesplq.h"
#include "util/ferm/map_obj/map_obj_aggregate_w.h"
#include "util/ferm/map_obj/map_obj_factory_w.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

#ifdef BUILD_SB

namespace Chroma 
{ 
  namespace InlineCreateColorVecsSuperbEnv 
  {
    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      if (inputtop.count("colorvec_files") == 1)
	read(inputtop, "colorvec_files", input.colorvec_files);
      read(inputtop, "colorvec_out", input.colorvec_out);
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "colorvec_files", input.colorvec_files);
      write(xml, "colorvec_out", input.colorvec_out);

      pop(xml);
    }

    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::Param_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "num_vecs", input.num_vecs);
      read(inputtop, "decay_dir", input.decay_dir);
      input.write_fingerprint = false;
      if (inputtop.count("write_fingerprint") == 1)
	read(inputtop, "write_fingerprint", input.write_fingerprint);
      input.link_smear = readXMLGroup(inputtop, "LinkSmearing", "LinkSmearingType");

      if (inputtop.count("phase") == 1)
      {
	read(inputtop, "phase", input.phase);
      }
      else
      {
	input.phase.resize(Nd - 1);
	for (int i = 0; i < Nd - 1; ++i)
	  input.phase[i] = 0;
      }
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t& out)
    {
      push(xml, path);

      write(xml, "num_vecs", out.num_vecs);
      write(xml, "decay_dir", out.decay_dir);
      write(xml, "write_fingerprint", out.write_fingerprint);
      write(xml, "phase", out.phase);
      xml << out.link_smear.xml;

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params& input)
    {
      Params tmp(xml, path);
      input = tmp;
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params& input)
    {
      push(xml, path);
    
      write(xml, "Param", input.param);
      write(xml, "NamedObject", input.named_obj);

      pop(xml);
    }
  } // namespace InlineCreateColorVecsSuperbEnv 


  namespace InlineCreateColorVecsSuperbEnv 
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
      
    const std::string name = "CREATE_COLORVECS_SUPERB";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	success &= MapObjectWilson4DEnv::registerAll();
	registered = true;
      }
      return success;
    }


    //-------------------------------------------------------------------------
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

	// Parameters for source construction
	read(paramtop, "Param", param);

	// Read in the output propagator/source configuration info
	read(paramtop, "NamedObject", named_obj);

	// Possible alternate XML file pattern
	if (paramtop.count("xml_file") != 0) 
	{
	  read(paramtop, "xml_file", xml_file);
	}
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
	QDP_abort(1);
      }
    }



    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	std::string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "CreateColorSuperbVecs");
	write(xml_out, "update_no", update_no);
	write(xml_out, "xml_file", xml_file);
	pop(xml_out);

	XMLFileWriter xml(xml_file);
	func(update_no, xml);
      }
      else
      {
	func(update_no, xml_out);
      }
    }


    // Real work done here
    void 
    InlineMeas::func(unsigned long update_no,
		     XMLWriter& xml_out) 
    {
      START_CODE();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      // Test and grab a reference to the gauge field
      multi1d<LatticeColorMatrix> u;
      XMLBufferWriter gauge_xml;
      try
      {
	u = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
	TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
	QDP_abort(1);
      }
      catch (const std::string& e) 
      {
	QDPIO::cerr << name << ": std::map call failed: " << e << std::endl;
	QDP_abort(1);
      }

      push(xml_out, "CreateColorVecsSuperb");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": Create color vectors" << std::endl;

      proginfo(xml_out);    // Print out basic program info

      // Write out the input
      write(xml_out, "Input", params);

      // Write out the config header
      write(xml_out, "Config_info", gauge_xml);

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      // Calculate some gauge invariant observables just for info.
      MesPlq(xml_out, "Observables", u);

      // Parse the phase
      if (params.param.phase.size() != Nd - 1)
      {
	QDPIO::cerr << "phase tag should have " << Nd - 1 << " components" << std::endl;
	QDP_abort(1);
      }
      SB::Coor<Nd - 1> phase;
      for (int i = 0; i < Nd - 1; ++i)
      {
	phase[i] = params.param.phase[i];
	if (std::fabs(phase[i] - params.param.phase[i]) > 0)
	  std::runtime_error("phase should be integer");
      }

      // Compute and store the colorvecs
      const int Lt = Layout::lattSize()[3];
      SB::createColorvecStorage(params.named_obj.colorvec_out, params.param.link_smear, u, 0, Lt,
				params.param.num_vecs, params.param.write_fingerprint,
				params.param.write_fingerprint, phase,
				params.named_obj.colorvec_files);

      pop(xml_out); // CreateColorVecsSuperb

      snoop.stop();
      QDPIO::cout << name << ": total time = " << snoop.getTimeInSeconds() << " secs" << std::endl;

      QDPIO::cout << name << ": ran successfully" << std::endl;

      END_CODE();
    } 

  }

} // namespace Chroma

#endif // BUILD_SB
