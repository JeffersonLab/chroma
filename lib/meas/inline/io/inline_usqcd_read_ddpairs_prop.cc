// $Id: inline_usqcd_read_ddpairs_prop.cc,v 3.3 2008-04-28 20:23:46 bjoo Exp $
/*! \file
 * \brief Inline task to read an object from a named buffer
 *
 * Named object writing
 */

#include "chromabase.h"

#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/inline_usqcd_read_ddpairs_prop.h"
#include "meas/inline/io/named_objmap.h"
#include "util/info/unique_id.h"
#include "util/ferm/transf.h"
#include "util/ft/sftmom.h"

using namespace QDP;

namespace Chroma 
{ 
  namespace InlineUSQCDReadDDPairsPropEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineUSQCDReadDDPairsProp(InlineUSQCDReadDDPairsPropParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "USQCD_READ_DD_PAIRS_PROP";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	// Inline measurement
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);

	registered = true;
      }
      return success;
    }
  }


  //! Object buffer
  void write(XMLWriter& xml, const string& path, const InlineUSQCDReadDDPairsPropParams& p)
  {
    push(xml, path);
    write(xml, "Frequency", p.frequency);
    push(xml, "Param");
    write(xml, "InputFile", p.input_file_name);
    write(xml, "SourceXML", p.source_xml);
    write(xml, "Precision", p.precision);
    pop(xml); // Param
    
    push(xml, "NamedObject");
    write(xml, "source_id", p.source_id); // ID Of source to create
    write(xml, "prop_id", p.prop_id); // ID Of prop to create
    pop(xml);

    if ( p.xml_file != "") { 
      write(xml, "xml_file", p.xml_file);
    }
    pop(xml);
  }

  // Param stuff
  InlineUSQCDReadDDPairsPropParams::InlineUSQCDReadDDPairsPropParams() { frequency = 0; xml_file=""; }

  InlineUSQCDReadDDPairsPropParams::InlineUSQCDReadDDPairsPropParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;


      read(paramtop, "Param/InputFile", input_file_name);
      read(paramtop, "Param/Precision", precision); 

      if ( (precision != "single") && (precision != "double")) { 
	QDPIO::cout << "Precision must be either single or double" << endl;
	QDP_abort(1);
      }

      if ( paramtop.count("Param/SourceXML") == 1 ) { 

	XMLReader prop_header(paramtop, "Param/SourceXML");
	ostringstream dummy;
	prop_header.printCurrentContext(dummy);
	source_xml = dummy.str();
      }
      else { 
	QDPIO::cout << "NO source XML in Parameter file. Will look in Prop Info" << endl;
	source_xml == "";
      }

      
      read(paramtop, "./NamedObject/source_id", source_id);
      read(paramtop, "./NamedObject/prop_id", prop_id);

      if( paramtop.count("xml_file") == 1 ) { 
	read(paramtop, "./xml_file", xml_file);
      }
      else { 
	xml_file="";
      }
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << __func__ << ": caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }




  void 
  InlineUSQCDReadDDPairsProp::operator()(unsigned long update_no,
				    XMLWriter& xml_out) 
  {
    START_CODE();
    if( params.xml_file == "" ) { 
      func(update_no, xml_out);
    }
    else { 
      XMLFileWriter separate_output(params.xml_file);
      func(update_no, separate_output);
    }

    END_CODE();
  } 


  void InlineUSQCDReadDDPairsProp::func(unsigned long update_no, XMLWriter& xml_out)
  {
    START_CODE();

    push(xml_out, "USQCDReadDDPairsProp");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlineUSQCDReadDDPairsPropEnv::name << ":" << endl;
    StopWatch swatch;
    QDPIO::cout << "Attempt to read object name = " << params.input_file_name;

    XMLBufferWriter testout;
    write(testout, "InputParams", params);
   
    QDPIO::cout << "Echo Of Input Params" << endl;
    QDPIO::cout << testout.str();

    /* Do the IO here */
    /* Create the named objects */
    try { 
      TheNamedObjMap::Instance().create<LatticePropagator>(params.source_id);
      TheNamedObjMap::Instance().create<LatticePropagator>(params.prop_id);
    }
    catch(...) { 
      QDPIO::cerr << "Failed to create the source and prop" <<endl;
      QDP_abort(1);
    }

    /* Get references to the data */
    LatticePropagator& source_ref=TheNamedObjMap::Instance().getData<LatticePropagator>(params.source_id);
    LatticePropagator& prop_ref=TheNamedObjMap::Instance().getData<LatticePropagator>(params.prop_id);

    /* Open the file */
    XMLReader file_xml;
    QDPFileReader input_file(file_xml, 
			     params.input_file_name,
			     QDPIO_SERIAL);


    push(xml_out, "FileXML"); 
    xml_out << file_xml;
    pop(xml_out);

    LatticeFermion tmp;
    LatticeFermionF tmpF;
    LatticeFermionD tmpD;

    // Cycle through the spin color components 
    for(int i=0; i < 24; i++) { 

      // For the record XML
      XMLReader rec_xml; 

      // To get the color and spin values
      int col, spin;

      // Is the record under consideration a source?
      bool issource=false;
      
      // Read the record depending on the precision and cast 
      // to base precision
      if( params.precision == "single" ) {
	read(input_file, rec_xml, tmpF);
	tmp=tmpF; // Truncate to my inner floating point if necessary
      }
      else {
	// Double
	read(input_file, rec_xml, tmpD);
	tmp=tmpD; // Truncate to my inner floating point if necessary
      }
       

      // Try and grok source color and spin info
      try { 
	read(rec_xml, "/usqcdSourceInfo/info/spin", spin);
	read(rec_xml, "/usqcdSourceInfo/info/color", col);

	// Found source record spin and color info, so assume it is a source
	issource = true;
	QDPIO::cout << "Record is source record with spin="<< spin << " color= "<< col << endl;

      }
      catch(const std::string& e) { 
	// We didn't find any source info... May still be a prop component
	// Otherwise we didn't find a source
	issource = false;
      }

      // If it is a source, shove it off into the source object
      if( issource ) { 
	FermToProp(tmp, source_ref, col, spin);
      }

      else { 
	
	// The thingie was not a source. Is it a Prop component?
	try { 
	  // Look for propagator spin and color info
	  read(rec_xml, "/usqcdInfo/spin", spin);
	  read(rec_xml, "/usqcdInfo/color", col);
	  QDPIO::cout << "Record is sink record with spin=" << spin << " color = "<< col << endl;

	  // If you find it, it is a prop and we can shove it off 
	  FermToProp(tmp, prop_ref, col, spin);
	  
	}
	catch(const std::string& e) { 

	  // We didn't find it so neither a source nor a prop
	  // and we just quietly barf
	  QDPIO::cout << "Caught exception reading XML" << e << endl;
	  QDP_abort(1);
	}
      }
    }


    // We will also need source and prop record XMLs. To conform
    // with what we do now, I give them a unique ID
    XMLBufferWriter source_rec_xml;
    XMLBufferWriter prop_rec_xml;
    push(source_rec_xml, "source");
    write(source_rec_xml, "id", uniqueId()); 
    pop(source_rec_xml);
    
    push(prop_rec_xml, "propagator");
    write(prop_rec_xml, "id", uniqueId());
    pop(prop_rec_xml);

    // Now we should have all the objects, just deal with the XML
    // First look for a header in the <info> tag
    bool found_header = false;

    try { 
      // Try and open a reader on our header
      XMLReader prop_header_xml(file_xml, "//info/Propagator");

      // If we don't exception then we found the header
      found_header = true;
      QDPIO::cout << "Found header in Prop File <info> record" << endl;

      // Try to grok the Source header in it the same way
      XMLReader source_header_xml(prop_header_xml, "./PropSource");
      QDPIO::cout << "Found source header too" << endl;

      // Set the source and prop XMLs
      TheNamedObjMap::Instance().get(params.prop_id).setFileXML(prop_header_xml);      
      TheNamedObjMap::Instance().get(params.prop_id).setRecordXML(prop_rec_xml);      

      TheNamedObjMap::Instance().get(params.source_id).setFileXML(source_header_xml);      
      TheNamedObjMap::Instance().get(params.source_id).setRecordXML(source_rec_xml);      
      // And we are done
    }
    catch(const std::string& e) {

      // We didn't find the header
      QDPIO::cout <<  "Didn't find header in info tag. Looking at user supplied" << endl ; 
      found_header = false;
    }

    // If we didn't find the header in <info> Look in the source file XML string
    if( ! found_header ) { 
      try { 
	// Turn string into istream
	istringstream is(params.source_xml);
	// Open a reader on it
	XMLReader source_top(is);

	// Try and set the context to be our Header
	XMLReader prop_header_xml(source_top, "//Propagator");

	// If we don't exception we found it
	QDPIO::cout << "Found header in User Supplied SourceXML" << endl;

	// Look for the source XML 
	XMLReader source_header_xml(prop_header_xml, "./PropSource");

	// If we don't exception we found it
	QDPIO::cout << "Found source header too" << endl;
	
	// Set the headers
	TheNamedObjMap::Instance().get(params.prop_id).setFileXML(prop_header_xml);      
	TheNamedObjMap::Instance().get(params.prop_id).setRecordXML(prop_rec_xml);      
	
	TheNamedObjMap::Instance().get(params.source_id).setFileXML(source_header_xml);      
	TheNamedObjMap::Instance().get(params.source_id).setRecordXML(source_rec_xml);      
      }
      catch(const std::string& e) { 
	// We've exceptioned so the XML is neither in <info> tag or in the input  params. Barf with message
	QDPIO::cout << "Caught exception while parsing user supplied possibly fake source XML: " << e << endl;
	QDPIO::cout << "If you got here it means that the necessary metadata for the prop you are trying to read is neither in the supplied binary neither in your parameters. Your file XML is: " << endl;
	ostringstream ostr;
	file_xml.printCurrentContext(ostr);
	QDPIO::cout << ostr << endl;

	QDPIO::cout << "Your source XML parameter in the param file is: " << endl;
	QDPIO::cout << params.source_xml << endl;

	QDP_abort(1);
      }
    }
      


    // Print out debugging info about the source correlator...
    {
      SftMom phases(0, true, Nd-1);
      multi1d<Double> source_corr = sumMulti(localNorm2(source_ref),
					   phases.getSet());

      multi1d<Double> prop_corr = sumMulti(localNorm2(prop_ref),
					   phases.getSet());

      write(xml_out, "source_corr", source_corr);
      write(xml_out, "prop_corr", prop_corr);
    }



    QDPIO::cout << InlineUSQCDReadDDPairsPropEnv::name << ": ran successfully" << endl;

    pop(xml_out);  // qio_read_named_obj

    END_CODE();
}
};
