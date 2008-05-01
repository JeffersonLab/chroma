// $Id: inline_usqcd_write_ddpairs_prop.cc,v 1.3 2008-05-01 19:32:56 bjoo Exp $
/*! \file
 * \brief Inline task to write an object from a named buffer
 *
 * Named object writing
 */

#include "chromabase.h"
#include "handle.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/inline_usqcd_write_ddpairs_prop.h"
#include "meas/inline/io/named_objmap.h"
#include "io/enum_io/enum_io.h"
#include "io/qprop_io.h"
#include "meas/sources/source_const_factory.h"
#include "util/ferm/transf.h"
#include "qio.h"
#include "util/ft/sftmom.h"

using namespace QDP;

namespace Chroma 
{

  namespace InlineUSQCDWriteDDPairsPropEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineUSQCDWriteDDPairsProp(InlineUSQCDWriteDDPairsPropParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "USQCD_WRITE_DD_PAIRS_PROP";

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

    const std::string sinkRecordName="usqcdInfo";

    // Functions for QIO
    //-----------------------------------------
    static int get_node_number(const int coord[])
    {
      multi1d<int> crd(Nd);
      crd = coord;   // an array copy
      int node = Layout::nodeNumber(crd);
      return node;
    }
    
    static int get_node_index(const int coord[])
    {
      multi1d<int> crd(Nd);
      crd = coord;   // an array copy
      int linear = Layout::linearSiteIndex(crd);
      return linear;
    }
    
    static void get_coords(int coord[], int node, int linear)
    {
      multi1d<int> crd = Layout::siteCoords(node, linear);
      for(int i=0; i < Nd; ++i)
	coord[i] = crd[i];
    }
    
    static int get_sites_on_node(int node) 
    {
      return Layout::sitesOnNode();
    }

    typedef PSpinVector< PColorVector< RComplex<REAL32>, 3>, 4> FermF;
    typedef PSpinVector< PColorVector< RComplex<REAL64>, 3>, 4> FermD;

  }


  //! Object buffer
  void write(XMLWriter& xml, const string& path, const InlineUSQCDWriteDDPairsPropParams& p)
  {
    push(xml, path);
    write(xml, "Frequency", p.frequency);
    push(xml, "Param");
    write(xml, "OutputFile", p.output_file_name);
    write(xml, "OutputVolfmt", p.qio_volfmt);
    write(xml, "Precision", p.precision);

    pop(xml); // Param
    
    push(xml, "NamedObject");
    write(xml, "prop_id", p.prop_id); // ID Of prop to create
    write(xml, "gauge_id", p.gauge_id);
    pop(xml);

    if ( p.xml_file != "") { 
      write(xml, "xml_file", p.xml_file);
    }
    pop(xml);
  }

  // Param stuff
  InlineUSQCDWriteDDPairsPropParams::InlineUSQCDWriteDDPairsPropParams() { frequency = 0; xml_file=""; }

  InlineUSQCDWriteDDPairsPropParams::InlineUSQCDWriteDDPairsPropParams(XMLReader& xml_in, const std::string& path) 
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
      read(paramtop, "Param/Precision", precision);


      read(paramtop, "./NamedObject/prop_id", prop_id);
      read(paramtop, "./NamedObject/gauge_id", gauge_id);

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
  InlineUSQCDWriteDDPairsProp::operator()(unsigned long update_no,
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


  void InlineUSQCDWriteDDPairsProp::func(unsigned long update_no, XMLWriter& xml_out)
  {
    START_CODE();

    push(xml_out, "USQCDWriteDDPairsProp");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlineUSQCDWriteDDPairsPropEnv::name << ":" << endl;
    StopWatch swatch;
    swatch.start();

    QDPIO::cout << "Attempt to write object name = " << params.prop_id << endl;
    QDPIO::cout << "Output file = " << params.output_file_name << endl;

    /* ================== Routine to write the prop here =============== */
    
    /* ---------------------------------------------
     *  Step 1: Get hold of the prop and gauge field
     * --------------------------------------------- */
    
    /* Local references are used to try and see if the gets fail.
       If not, a wider scope reference will be used to bind the prop
       and u. This is to avoid a very very long try {} catch block */
    try {
      LatticePropagator& trial_prop=TheNamedObjMap::Instance().getData<LatticePropagator>(params.prop_id);
      
      const multi1d<LatticeColorMatrix>& u_trial = 
	TheNamedObjMap::Instance().getData<multi1d<LatticeColorMatrix> >(params.gauge_id);
    }
    catch(...) { 
      QDPIO::cout << "Could not get the prop from the named ObjectMap. Missing ID"<< params.prop_id << endl;
    }

    // OK If we're here, the prop and gauge field are in the store
    // so bind them 
    const LatticePropagator& theProp=
      TheNamedObjMap::Instance().getData<LatticePropagator>(params.prop_id);

    const multi1d<LatticeColorMatrix>& u = 
      TheNamedObjMap::Instance().getData<multi1d<LatticeColorMatrix> >(params.gauge_id);

    
    // Now get the Record XML for the prop where all the important stuff is
    XMLReader thePropHeader;
    TheNamedObjMap::Instance().get(params.prop_id).getRecordXML(thePropHeader);


    // Grab the header for the source
    XMLReader theSourceHeader(thePropHeader,"//PropSource");

    /* -------------------------------------------------------------
     * Make the source here. Use the XML for the source, to return the
     * appropriate source creation function. We apply this to - u - 
     * to get the appropriate source. WARNING WARNING WARNING: If the
     * user uses the wrong gauge field this can all go horribly pear
     * shaped. 
     *------------------------------------------------------------------*/

    LatticePropagator theSource;  // the source in floating precision
                                  // We'll cast to single/double later
    try
    {
      // Read the factory key for the source
      std::string source_type;
      read(theSourceHeader, "Source/SourceType", source_type);

      QDPIO::cout << "Source = " << source_type << endl;

      // Get the constructor function-object from the factory
      Handle< QuarkSourceConstruction<LatticePropagator> >
	sourceConstruction(ThePropSourceConstructionFactory::Instance().createObject(source_type,theSourceHeader,"Source")); 

      // Apply the source construction 
      theSource = (*sourceConstruction)(u);   // And we are done.

    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << InlineUSQCDWriteDDPairsPropEnv::name << ": Caught Exception creating source: " << e << endl;
      QDP_abort(1);
    }

    /* -------------------------------------------------------------------
     * OK. To output the source timeslice with QIO, we need to get the 
     * timeslice information out of the source header. The way we'll 
     * do this is just to make a source construction data object from 
     * th XML header -- ie data bind the header -- and interrogate the
     * structure 
     *------------------------------------------------------------------ */

    int t_slice;
    int j_decay;
    PropSourceConst_t source_info;
    {
      // Data bind the header info
      read(thePropHeader, "//PropSource", source_info);

      QDPIO::cout << "Source timeslice is " << source_info.t_source << endl;

      // Get timeslice and decay direction
      t_slice=source_info.t_source;
      j_decay=source_info.j_decay;
    }
    

    // Bind our Prop header to the relevant QPropIO datatype for ease of writing
    Propagator_t prop_header;
    read(thePropHeader, "//Propagator", prop_header);
    /* ------------------------------------------------------------------
     * Right. We're almost ready to write now. We'll need a couple of 
     * Stock XML snippets, like the File XML with our 'record' XML header
     * info being jammed into the <info> tag. Yuck! 
     * ------------------------------------------------------------------ */

    XMLBufferWriter file_xml;
    push(file_xml, "usqcdPropFile");
    write(file_xml, "version", "1.0");
    write(file_xml, "type", "USQCD_DiracFermion_Source_Sink_Pairs");
    push(file_xml, "info");
    write(file_xml, "Propagator", prop_header);
    pop(file_xml); // Info
    pop(file_xml); // usqcdPropFile

    /* Echo back for human checking */
    QDPIO::cout << "File XML is:" << endl;
    QDPIO::cout << file_xml.str();

    /* --------------------------------------------------------------------
     * For each source, we'll also need a source record XML. It appears to
     * be content free, but I can always jam the source header file into
     * its <info> tag
     * ------------------------------------------------------------------- */
    XMLBufferWriter source_xml;
    push(source_xml, "usqcdSourceInfo");
    write(source_xml, "version", "1.0");
    push(source_xml, "info");
    
    // Jam in our source info by writing (marshalling?) the PropSourceConst_t
    //    write(source_xml, "PropSource", source_info);
    pop(source_xml);
    pop(source_xml);

    // Display for the user...
    QDPIO::cout << "Source XML is:" << endl;
    QDPIO::cout << source_xml.str();

    /* --------------------------------------------------------------------
     * OK. Here is where we get in a little trouble. I can't very well use
     * QDP here, because it is not geared up to do 'hypercube' records.
     * Maybe I could change that later(?). So I have to go to the
     * QIO calls directly. What a nigthmare... 
     * ------------------------------------------------------------------ */

    // QIO Data structures
    QIO_Layout layout;
    int latsize[Nd];
    
    // Get the lattice size
    for(int m=0; m < Nd; ++m)
      latsize[m] = Layout::lattSize()[m];
    
    layout.node_number = &get_node_number;
    layout.node_index  = &get_node_index;
    layout.get_coords  = &get_coords;
    layout.num_sites = &get_sites_on_node;
    layout.latsize = latsize;
    layout.latdim = Nd; 
    layout.volume = Layout::vol(); 
    layout.sites_on_node = Layout::sitesOnNode(); 
    layout.this_node = Layout::nodeNumber(); 
    layout.number_of_nodes = Layout::numNodes(); 
    
    // Create the "xml string for the file"
    QIO_String* xml_c = QIO_string_create();
    
    // Put the file XML in it
    QIO_string_set(xml_c, file_xml.str().c_str());
    
    // Wrappers over simple ints
    int volfmt;
    switch(params.qio_volfmt) {
      
      case QDPIO_SINGLEFILE:
	volfmt = QIO_SINGLEFILE;
	break;
	
    case QDPIO_MULTIFILE:
      volfmt = QIO_MULTIFILE;
      break;
      
    case QDPIO_PARTFILE:
      volfmt = QIO_PARTFILE;
      break;
      
    default: 
      QDPIO::cerr << "Unknown value for qio_volfmt " << params.qio_volfmt << endl;
      QDP_abort(1);
      return;
    }    
    
    // QIO_Oflag
    QIO_Oflag oflag;
    oflag.serpar = QIO_SERIAL;
    oflag.mode = QIO_TRUNC; 
    oflag.ildgstyle = QIO_ILDGNO;
    oflag.ildgLFN=NULL;
    
    // This is recent, but allowed to be null
    QIO_Filesystem fs; 
    fs.my_io_node = NULL;
    fs.master_io_node = NULL;
    
    QIO_Writer* qio_out;
    
    // Open the qio_out writer
    if ((qio_out = QIO_open_write(xml_c, params.output_file_name.c_str(),
				  volfmt, 
				  &layout, 
				  &fs, &oflag)) == NULL ) {
      QDPIO::cerr << "QIO Error failed to open file " << params.output_file_name << endl;
      QDP_abort(1);  // just bail. Not sure I want this. This is not stream semantics
    }
    QIO_string_destroy(xml_c);
    
    /* Set the QIO precision and datatype. I hope the typesize is right...*/
    std::string usqcd_datatype;
    std::string usqcd_precision;
    int typesize;
    int wordsize;
    if( params.precision == "single" ) {
      usqcd_datatype = "USQCD_F3_DiracFermion";
      usqcd_precision ="F";
      typesize=(int)sizeof(FermF);
      wordsize=sizeof(WordType<FermF>::Type_t);
    }
    else if( params.precision == "double" ) { 
      usqcd_datatype = "USQCD_D3_DiracFermion";
      usqcd_precision= "D";
      typesize=(int)sizeof(FermD);
      wordsize=sizeof(WordType<FermD>::Type_t);
    }
    else { 
      QDPIO::cerr << "Invalid value for precision: "<< params.precision << endl;
      QDPIO::cerr << "The value must be either \"single\" or \"double\" " << endl;
      QDP_abort(1);
    }


    
    /* Set up the hypercube info */
    int lower[Nd];
    int upper[Nd];
    
    /* This should set the lower corner to be the origin and the max
       corner to be the largest accessible */
    for(int i=0; i < Nd; i++) { 
      lower[i] = 0;
      upper[i] = latsize[i]-1;
    } 

    /* We overwrite the time direction coordinates of the slice
       with the source timeslice */
    upper[j_decay] = t_slice;
    lower[j_decay] = t_slice;
    
    
    // Prepare the source record info. This won't change
    QIO_RecordInfo *source_record_info = QIO_create_record_info(QIO_HYPER, 
								lower,
								upper,
								Nd,
								(char *)(usqcd_datatype.c_str()),
								(char *)(usqcd_precision.c_str()),
								3, 4, 
								typesize,
								1);

    if( source_record_info == NULL ){ 
      QDPIO::cerr << "Failed to create record info" << endl;
      QDP_abort(1);
    }
    
    /* Prepare the source record XML. This won't change. */
    QIO_String* source_string = QIO_string_create();
    QIO_string_set(source_string, source_xml.str().c_str());
    if( source_string == NULL ){ 
      QDPIO::cerr << "Failed to set source file_xml" << endl;
      QDP_abort(1);
    }


    /* Make the record info for the propagator spinor components. */
    QIO_RecordInfo* spinor_record_info = QIO_create_record_info(QIO_FIELD,
							       NULL,
							       NULL,
							       Nd,
							       (char *)(usqcd_datatype.c_str()),
							       (char *)(usqcd_precision.c_str()),
							       3, 4, 
							       typesize,
							       1);

    if( spinor_record_info == NULL ){ 
      QDPIO::cerr << "Failed to create record info" << endl;
      QDP_abort(1);
    }

   
    // OK We're done with the headers and record infos.

    // Next we have to loop through the spin-color components and write
    for(int spin=0; spin < Ns; spin++) { 
      for(int color=0; color < Nc; color++) { 

	/* ---------------------------------------------------------------------------
	 * Make some XML for the spinor -- I couldn't pre-prepare this as it contains
	 * the spin color information
	 * --------------------------------------------------------------------------- */

	XMLBufferWriter spinor_xml;
	push(spinor_xml, sinkRecordName);
	write(spinor_xml, "version", "1.0");
	write(spinor_xml, "spin", spin);
	write(spinor_xml, "color", color);
	push(spinor_xml, "info");
	//	spinor_xml << thePropHeader ;
	pop(spinor_xml); // Info
	pop(spinor_xml); // usqcdPropInfo

	QIO_String* spinor_xml_string = QIO_string_create();
	QIO_string_set(spinor_xml_string, spinor_xml.str().c_str());

	/* -------------------------------------------------------------------------
	 * Write the source first.
	 * ------------------------------------------------------------------------- */

	//
	// Bring the source into the temporary fermion. 
	//
	LatticeFermion tmpFerm=zero;
	PropToFerm(theSource, tmpFerm, color, spin);
	
	// Now depending on precision, cast and write
	if ( params.precision == "single" ) { 
	  LatticeDiracFermionF3 ferm_out = tmpFerm ;           // Single Precision cast
    
	  {
	    SftMom phases(0, true, Nd-1);
	    multi1d<Double> source_corr = sumMulti(localNorm2(ferm_out),
						   phases.getSet());
	    write(xml_out, "output_source_corr", source_corr);
	  }

	  // Write 
	  QDPIO::cout << "About to write single prec. source record....";
	  int err_code = QIO_write(qio_out,source_record_info,source_string,
				   (QDPOLatticeFactoryGet<FermF>),
				   typesize,
				   wordsize,
				   (void *)(ferm_out.getF()));

	  // Error check
	  if (err_code != QIO_SUCCESS) { 
	    QDPIO::cout << "Failed!" << endl;
	    QDP_abort(1);
	  }
	  else {
	    QDPIO::cout << "Successful!" << endl;
	  }
	}
	else {
	  LatticeDiracFermionD3 ferm_out = tmpFerm;             // Double precision cast
    
	  /* Write */
	  QDPIO::cout << "About to write double prec. source record....";
	  int err_code = QIO_write(qio_out,source_record_info,source_string,
				   (QDPOLatticeFactoryGet<FermD>),
				   typesize,
				   wordsize,
				   (void *)(ferm_out.getF()));

	  /* Error check */
	  if (err_code != QIO_SUCCESS) { 
	    QDPIO::cerr << "Failed!" << endl;
	    QDP_abort(1);
	  }
	  else { 
	    QDPIO::cout << "Successful!" << endl;
	  }
	}

	/* --------------------------------------------------------------------------------
	 * Done with source. Now do the prop component
	 * -------------------------------------------------------------------------------- */

	// Bring the prop component into the tmpFerm
	PropToFerm(theProp, tmpFerm, color, spin);

	/* Now depending on the precision, cast and write */
	if ( params.precision == "single" ) { 
	  LatticeDiracFermionF3 ferm_out = tmpFerm ;  // Single precision Cast 
    
	  /* Write */
	  QDPIO::cout << "About to write single prec. prop component record....";
	  int err_code = QIO_write(qio_out,spinor_record_info, spinor_xml_string,
				   (QDPOLatticeFactoryGet<FermF>),
				   typesize,
				   wordsize,
				   (void *)(ferm_out.getF()));

	  /* Error check */
	  if (err_code != QIO_SUCCESS) { 
	    QDPIO::cout << "Failed!" << endl;
	    QDP_abort(1);
	  }
	  else {
	    QDPIO::cout << "Successful!" << endl;
	  }
	}
	else {
	  LatticeDiracFermionD3 ferm_out = tmpFerm;  // Double precision cast

	  /* Write */
	  QDPIO::cout << "About to write double prec. prop component record....";
	  int err_code = QIO_write(qio_out, spinor_record_info, spinor_xml_string,
				   (QDPOLatticeFactoryGet<FermD>),
				   typesize,
				   wordsize,
				   (void *)(ferm_out.getF()));

	  /* Error check */
	  if (err_code != QIO_SUCCESS) { 
	    QDPIO::cerr << "Failed!" << endl;
	    QDP_abort(1);
	  }
	  else { 
	    QDPIO::cout << "Successful!" << endl;
	  }
	}

	/* Destroy QIO string object. My own XML buffers will clean themselves up magically 
	   so "Nyaaah !!!!!" */
	QIO_string_destroy(spinor_xml_string);

      } // color
    } // spin

    
    /* Destroy the record infos */
    QIO_destroy_record_info(spinor_record_info);	
    QIO_destroy_record_info(source_record_info);

    /* Destroy the source string */
    QIO_string_destroy(source_string);
    
    
    /* Close file at the end. */
    QIO_close_write(qio_out);
    
    /* And we are done. Thank heavens! That was singularly unpleasant. 446 odd lines
       that should have taken much fewer */
    swatch.stop();

    QDPIO::cout << InlineUSQCDWriteDDPairsPropEnv::name << ": total time = " << swatch.getTimeInSeconds() <<" secs" << endl;
    QDPIO::cout << InlineUSQCDWriteDDPairsPropEnv::name << ": ran successfully" << endl;


    pop(xml_out);  // qio_write_named_obj

    END_CODE();
}
};
