// $Id: inline_qpropqio_w.cc,v 3.2 2006-11-17 02:17:31 edwards Exp $
/*! \file
 * \brief Inline measurement of qpropqio
 *
 * Form-factor measurements
 */

#include "meas/inline/hadron/inline_qpropqio_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "io/enum_io/enum_qdpvolfmt_io.h"
#include "io/qprop_io.h"

namespace Chroma 
{ 
  namespace InlineQpropQIOEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineQpropQIO(InlineQpropQIOParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "QPROPQIO";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
  }


  //! Propagator parameters
  void read(XMLReader& xml, const string& path, InlineQpropQIOParams::Prop_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "prop_in_file", input.prop_in_file);
    read(inputtop, "prop_out_file", input.prop_out_file);
    read(inputtop, "prop_out_volfmt", input.prop_out_volfmt);
  }

  //! Propagator parameters
  void write(XMLWriter& xml, const string& path, const InlineQpropQIOParams::Prop_t& input)
  {
    push(xml, path);

    write(xml, "prop_in_file", input.prop_in_file);
    write(xml, "prop_out_file", input.prop_out_file);
    write(xml, "prop_out_volfmt", input.prop_out_volfmt);

    pop(xml);
  }


  // Reader for input parameters
  void read(XMLReader& xml, const string& path, InlineQpropQIOParams::Param_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    switch (version) 
    {
    case 1:
      /**************************************************************************/
      break;

    default :
      /**************************************************************************/

      QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
      QDP_abort(1);
    }
  }


  // Reader for input parameters
  void write(XMLWriter& xml, const string& path, const InlineQpropQIOParams::Param_t& param)
  {
    push(xml, path);

    int version = 1;

    write(xml, "version", version);

    pop(xml);
  }


  // Param stuff
  InlineQpropQIOParams::InlineQpropQIOParams() { frequency = 0; }

  InlineQpropQIOParams::InlineQpropQIOParams(XMLReader& xml_in, const std::string& path) 
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
      read(paramtop, "Prop", prop);
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  void
  InlineQpropQIOParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    // Parameters for source construction
    Chroma::write(xml_out, "Param", param);

    // Write out the output propagator/source configuration info
    Chroma::write(xml_out, "Prop", prop);

    pop(xml_out);
  }

  //--------------------------------------------------------------


  // Function call
  void 
  InlineQpropQIO::operator()(unsigned long update_no,
			     XMLWriter& xml_out) 
  {
    START_CODE();

    push(xml_out, "qpropqio");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << "QPROPQIO: propagator transformation utility" << endl;

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.write(xml_out, "Input");

  
    /*
     * Now read them thangs...
     * NOTE: only SCIDAC format is allowed !!!
     */
    XMLReader prop_in_record_xml, prop_in_file_xml;
    LatticePropagator  prop;
    readQprop(prop_in_file_xml, prop_in_record_xml, prop, 
	      params.prop.prop_in_file, QDPIO_SERIAL);

    push(xml_out,"SciDAC_propagator");
    write(xml_out, "File_xml", prop_in_file_xml);
    write(xml_out, "Record_xml", prop_in_record_xml);
    pop(xml_out);

    // Sanity check - write out the propagator (pion) correlator in the Nd-1 direction
    {
      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, Nd-1);

      multi1d<Double> prop_corr = sumMulti(localNorm2(prop), 
					   phases.getSet());

      push(xml_out, "Prop_correlator");
      write(xml_out, "prop_corr", prop_corr);
      pop(xml_out);
    }

    /*
     * Now write them thangs...
     */ 
    {
      XMLBufferWriter prop_out_file_xml;
      prop_out_file_xml << prop_in_file_xml;

      XMLBufferWriter prop_out_record_xml;
      prop_out_record_xml << prop_in_record_xml;
    
      // Write the source
      writeQprop(prop_out_file_xml, prop_out_record_xml, prop,
		 params.prop.prop_out_file, params.prop.prop_out_volfmt, 
		 QDPIO_SERIAL);
    }

    pop(xml_out);   // qpropqio
        
    QDPIO::cout << "QpropQIO ran successfully" << endl;

    END_CODE();
  } 

};
