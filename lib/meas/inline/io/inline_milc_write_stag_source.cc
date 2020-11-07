/*! \file
 * \brief Inline task to write an object from a named buffer
 *
 * Named object writing
 */


#include "chromabase.h"
#include "handle.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/inline_milc_write_stag_source.h"
#include "meas/inline/io/named_objmap.h"
#include "io/enum_io/enum_io.h"
#include "util/ferm/transf.h"
#include "qio.h"


using namespace QDP;

namespace Chroma 
{

  namespace InlineMILCWriteStagSourceEnv
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
    	  return new InlineMILCWriteStagSource(InlineMILCWriteStagSourceParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "MILC_WRITE_STAGGERED_SOURCE";

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


  // Some very private things
  namespace { 

    // The tag in the record XML for sinks. This needs to be decided on.
    // For now I go with the latest version of QIO which says "usqcdInfo"
    // as opposed to "usqcdPropInfo" in the manual. The reader can eat both
    // supposedly, so once they decide which one we can change to that.

    const std::string sinkRecordName="usqcdPropInfo";
  }


  //! Object buffer
  void write(XMLWriter& xml, const std::string& path, const InlineMILCWriteStagSourceParams& p)
  {
    push(xml, path);
    write(xml, "Frequency", p.frequency);
    push(xml, "Param");
    write(xml, "OutputFile", p.output_file_name);
    write(xml, "OutputVolfmt", p.qio_volfmt);
    bool parallel_io_choice = false;
      if( p.parallel_io == QDPIO_PARALLEL) {
      	parallel_io_choice = true;
      }
      write(xml, "parallel_io", parallel_io_choice);
    write(xml, "Precision", p.precision);

    pop(xml); // Param
    
    push(xml, "NamedObject");
    write(xml, "source_id", p.source_id); // ID Of prop to create
    write(xml, "t_slice", p.t_slice);
    pop(xml);

    if ( p.xml_file != "") { 
      write(xml, "xml_file", p.xml_file);
    }
    pop(xml);
  }

  // Param stuff
  InlineMILCWriteStagSourceParams::InlineMILCWriteStagSourceParams() { frequency = 0; xml_file=""; }

  InlineMILCWriteStagSourceParams::InlineMILCWriteStagSourceParams(XMLReader& xml_in, const std::string& path)
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;


      read(paramtop, "Param/OutputFile", output_file_name);
      read(paramtop, "Param/OutputVolfmt", qio_volfmt);
      bool parallel_io_choice = Layout::isIOGridDefined() && ( Layout::numIONodeGrid() > 1);

      if( paramtop.count("Param/parallel_io") == 1) {
    	  read(paramtop, "Param/parallel_io", parallel_io_choice);
      }

      read(paramtop, "Param/t_slice", t_slice);

      if( parallel_io_choice ) {
    	parallel_io = QDPIO_PARALLEL;
      }
      else {
    	parallel_io = QDPIO_SERIAL;
      }
      read(paramtop, "Param/Precision", precision);


      read(paramtop, "./NamedObject/source_id", source_id);

      if( paramtop.count("xml_file") == 1 ) { 
	read(paramtop, "./xml_file", xml_file);
      }
      else { 
	xml_file="";
      }
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << __func__ << ": caught Exception reading XML: " << e << std::endl;
      QDP_abort(1);
    }
  }




  void 
  InlineMILCWriteStagSource::operator()(unsigned long update_no,
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


  void InlineMILCWriteStagSource::func(unsigned long update_no, XMLWriter& xml_out)
  {
    START_CODE();
    if ( Nc != 3 ){    /* Code is specific to Ns=4 and Nc=3. */
      QDPIO::cerr<<" code only works for Nc=3 \n";
      QDP_abort(111) ;
    }
#if QDP_NC == 3


    push(xml_out, "MILCWriteStaggSource");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlineMILCWriteStagSourceEnv::name << ":" << std::endl;
    StopWatch swatch;
    swatch.start();

    QDPIO::cout << "Attempt to write object name = " << params.source_id << std::endl;
    QDPIO::cout << "On timeslice: " << params.t_slice << std::endl;
    QDPIO::cout << "Output file = " << params.output_file_name << std::endl;

    /* ================== Routine to write the prop here =============== */
    
    /* ---------------------------------------------
     *  Step 1: Get hold of the prop and gauge field
     * --------------------------------------------- */
    
    /* Local references are used to try and see if the gets fail.
       If not, a wider scope reference will be used to bind the prop
       and u. This is to avoid a very very long try {} catch block */
    try {
      LatticeStaggeredPropagator& trial_source=TheNamedObjMap::Instance().getData<LatticeStaggeredPropagator>(params.source_id);
          }
    catch(...) { 
      QDPIO::cout << "Could not get the source from the named ObjectMap. Missing ID"<< params.source_id << std::endl;
      QDP_abort(1);
    }

    // OK If we're here, the prop and gauge field are in the store
    // so bind them 
    const LatticeStaggeredPropagator& theSrc=
      TheNamedObjMap::Instance().getData<LatticeStaggeredPropagator>(params.source_id);


    XMLBufferWriter file_xml;
    push(file_xml, "ChromaMILCStaggeredSource");
    write(file_xml, "note", "Staggered Source written from Chroma. Metadata irrelevant");
    pop(file_xml);
    
    // Open the QIO output file
    QDPFileWriter qio_out(file_xml, params.output_file_name,
    		params.qio_volfmt,
			params.parallel_io,
			QDPIO_CREATE);


    multi1d<int> lower_left(Nd);
    multi1d<int> upper_right(Nd);
    
    for(int i=0; i < Nd-1; ++i) {
    	lower_left[i] = 0;
    	upper_right[i] = Layout::lattSize()[i]-1;
    }
    lower_left[Nd-1] = params.t_slice;
    upper_right[Nd-1] = params.t_slice;
    
    for(int col=0; col < 3; ++col) {
    	XMLBufferWriter record_xml;
    	push(record_xml, "usqcdKSPropInfo");
    	write(record_xml, "version", "1.0");
    	write(record_xml, "color", col);
    	write(record_xml, "info", " ");
    	pop(record_xml);

    	// Extract the component
    	LatticeStaggeredFermion tmpFerm;
    	PropToFerm(theSrc,tmpFerm,col);

    	// Now depending on precision, cast and write
    	if ( params.precision == "single" ) {
    	  LatticeStaggeredFermionF castFerm=(tmpFerm) ;           // Single Precision cast
    	  LatticeColorVectorF fermOut;
    	  FermToCv(castFerm,fermOut);

    	  // Write
    	  QDPIO::cout << "About to write single prec. source component: " << col << std::endl;
    	  QDPIO::cout << "TypeSize is: " << sizeof(PColorVector<RComplex<REAL32>,3>) << std::endl;
    	  write(qio_out,record_xml, fermOut, lower_left, upper_right);

    	}
    	else {
    	  LatticeStaggeredFermionD castFerm(tmpFerm);             // Double precision cast
    	  LatticeColorVector fermOut;
    	  FermToCv(castFerm,fermOut);

    	  /* Write */
    	  QDPIO::cout << "About to write double prec. source component." << col << std::endl;
    	  write(qio_out,record_xml,fermOut, lower_left, upper_right);
    	}
    } // End loop over color
    qio_out.close();
    pop(xml_out);
#endif
    END_CODE();
}
}

