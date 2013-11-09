// $Id: inline_usqcd_read_ddpairs_prop.cc,v 3.6 2009-03-19 17:10:02 mcneile Exp $
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
#include "io/qprop_io.h"

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
    write(xml, "PropXML", p.prop_xml);
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


      if ( paramtop.count("Param/PropXML") == 1 ) { 

	XMLReader prop_header(paramtop, "Param/PropXML");
	ostringstream dummy;
	prop_header.printCurrentContext(dummy);
	prop_xml = dummy.str();
      }
      else { 
	QDPIO::cout << "NO source XML in Parameter file. Will look in Prop Info" << endl;
	prop_xml = "";
      }

      
      read(paramtop, "./NamedObject/source_id", source_id);
      read(paramtop, "./NamedObject/prop_id", prop_id);

      if( paramtop.count("xml_file") == 1 ) { 
	read(paramtop, "./xml_file", xml_file);
      }
      else { 
	xml_file = "";
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
    swatch.reset();
    swatch.start();

    QDPIO::cout << "Attempt to read object name = " << params.input_file_name << endl;

   

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

    /* Open the file - read thef ile XML*/
    XMLReader file_xml;
    QDPFileReader input_file(file_xml, 
			     params.input_file_name,
			     QDPIO_SERIAL);


    push(xml_out, "FileXML"); 
    xml_out << file_xml;
    pop(xml_out);


    /* Try to bind the prop and source info up front, so you can know if the metadata is OK
       before reading in the whole file */

    /* Strategy: Look first in the file xml, and if you cant find the info there, look in the
       User supplied <PropXML> tags */
    Propagator_t prop_h_info;
    MakeSourceProp_t make_source_header;
	
    try  { 
      // go to the file <info> tag in the file XML
      XMLReader prop_info_xml(file_xml, "/usqcdPropFile/info");
      
      // check for the existence of a <Propagator> subtag
      if( prop_info_xml.count("./Propagator") != 0 ) { 

	// It's there
	QDPIO::cout << "Found header in Prop File <info> record" << endl;
	QDPIO::cout << "Attempting to bind" << endl;

	// Try to bind it
	read(prop_info_xml, "Propagator", prop_h_info);

      }
      else { 

	// The <info> tag doesn't contain what we want
	// Let's look at what the user specified in the <PropXML>
	// Tags...

	// Make a reader out of the prop_xml;
	istringstream is(params.prop_xml);
	XMLReader source_top(is);
      	
	QDPIO::cout << "Looking for <Propagator> tag in the user supplied XML" << endl;
	// Check if the tag is there
	if( source_top.count("/Propagator") == 0 ) { 
	  QDPIO::cout << "<Propagator> tag not found!" << endl;
	  QDP_abort(1);
	}

	// Yes, try to bind it.
	QDPIO::cout << "Attempting to bind" << endl;
	read(source_top, "/Propagator", prop_h_info);
	
      }

      // Make a source header. The info is there in the propagato
      QDPIO::cout << "<Propagator> header XML successfully bound" << endl;
      make_source_header.source_header = prop_h_info.source_header;
      make_source_header.gauge_header = prop_h_info.gauge_header;
    }
    catch(const std::string& e) { 
      // We've exceptioned so the XML is neither in <info> tag or in the input  params. Barf with message
      QDPIO::cout << "Caught exception: " << e << endl;

      QDP_abort(1);
    }
    catch(...) { 
      QDPIO::cout << "Caught generic exception. Most likely I failed to open a reader on th string you provided for  <PropXML>" << endl;

      QDP_abort(1);
    }
    
    {
      LatticeFermion tmpSource;
      LatticeFermion tmpProp;
    
      // Cycle through the spin color components 
      // Assume strict ordering of source spin pairs
      for(int i=0; i < 12; i++) { 

	// For the record XML
	XMLReader rec_src_xml, rec_prop_xml; 
	
	// To get the color and spin values
	int col, spin;
	
	
	// Read - this will do the precision casting into 
	// my precision for me

	// Read source
	read(input_file, rec_src_xml, tmpSource);
       	
	// read prop
	read(input_file, rec_prop_xml, tmpProp);
	
	// Look for spin/color info in the Propagator record.
	// This either has tag <usqcdPropInfo> or <usqcdInfo>
	// depending on your QIO version apparently.
	try { 
	  
	  // Find Spin Color Components
	  if ( rec_prop_xml.count("/usqcdPropInfo") != 0 ) { 
	    //
	    // As per Manual, the record tag is <usqcdPropInfo>
	    // 
	    read(rec_prop_xml, "/usqcdPropInfo/spin", spin);
	    read(rec_prop_xml, "/usqcdPropInfo/color", col);
	    QDPIO::cout << "Record is sink record with spin=" << spin << " color = "<< col << endl;
	    
	    // If you find it, it is a prop and we can shove it off 
	  }
	  else if( rec_prop_xml.count("/usqcdInfo" ) !=  0 ) { 
	    //
	    // As per QIO v2.3.4 and Sergey, the record tag is <usqcdInfo>
	    //
	    read(rec_prop_xml, "/usqcdInfo/spin", spin);
	    read(rec_prop_xml, "/usqcdInfo/color", col);
	    QDPIO::cout << "Record is sink record with spin=" << spin << " color = "<< col << endl;
	  }
	  else {
	    // As yet unforseen cases
	    QDPIO::cout << "Found neither usqcdInfo nor usqcdPropInfo tag in the propagator record XML" << endl;
	    QDPIO::cout << "Aborting" << endl;
	    QDP_abort(1);
	  }
	}
	catch(const std::string& e) { 
	  // We didn't find it so neither a source nor a prop
	  // and we just quietly barf
	  QDPIO::cout << "Caught exception reading XML" << e << endl;
	  QDP_abort(1);
	  
	}

	FermToProp(tmpSource, source_ref, col, spin);
	FermToProp(tmpProp,   prop_ref, col, spin);

      }
    }

    // We will also need source and prop record XMLs. To conform
    // with what we do now, I give them a unique ID
    XMLBufferWriter source_file_xml;
    XMLBufferWriter prop_file_xml;
    push(source_file_xml, "source");
    write(source_file_xml, "id", uniqueId()); 
    pop(source_file_xml);
    
    push(prop_file_xml, "propagator");
    write(prop_file_xml, "id", uniqueId());
    pop(prop_file_xml);


    // Create the usual Chroma Record XMLs
    // From the structures we bound earlier
    XMLBufferWriter prop_info_out;
    XMLBufferWriter source_info_out;
    write(prop_info_out, "Propagator", prop_h_info);
    write(source_info_out, "MakeSource", make_source_header);
      
    // Set the source and prop XMLs
    TheNamedObjMap::Instance().get(params.prop_id).setFileXML(prop_file_xml);      
    TheNamedObjMap::Instance().get(params.prop_id).setRecordXML(prop_info_out);      
    
    TheNamedObjMap::Instance().get(params.source_id).setFileXML(source_file_xml);      
    TheNamedObjMap::Instance().get(params.source_id).setRecordXML(source_info_out);      
      
        


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


    swatch.stop();
    QDPIO::cout << InlineUSQCDReadDDPairsPropEnv::name << ": total time = " << swatch.getTimeInSeconds() << " secs"<< endl;
    QDPIO::cout << InlineUSQCDReadDDPairsPropEnv::name << ": ran successfully" << endl;

    pop(xml_out);  // qio_read_named_obj

    END_CODE();
}
};

