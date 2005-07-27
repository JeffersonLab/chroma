// $Id: inline_building_blocks_w.cc,v 1.3 2005-07-27 16:23:52 edwards Exp $
/*! \file
 * \brief Inline construction of BuildingBlocks
 *
 * Building Blocks on forward and sequential props
 */

#include "meas/inline/hadron/inline_building_blocks_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/BuildingBlocks_w.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"

namespace Chroma 
{ 
  namespace InlineBuildingBlocksEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineBuildingBlocks(InlineBuildingBlocksParams(xml_in, path));
    }

    const std::string name = "BUILDING_BLOCKS";
    const bool registered = TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
  };


  //! Param input
  void read(XMLReader& xml, const string& path, InlineBuildingBlocksParams::Param_t& input)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    switch (version) 
    {
    case 1:
      /**************************************************************************/
      input.use_sink_offset = false;
      break;

    case 2:
      /**************************************************************************/
      read(paramtop, "use_sink_offset", input.use_sink_offset);
      break;

    default :
      /**************************************************************************/
      
      QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
      QDP_abort(1);
    }
    
    read(paramtop, "links_max", input.links_max);
    read(paramtop, "mom2_max", input.mom2_max);
    read(paramtop, "nrow", input.nrow);
  }

  //! Param write
  void write(XMLWriter& xml, const string& path, const InlineBuildingBlocksParams::Param_t& input)
  {
    push(xml, path);

    int version = 1;
    write(xml, "version", version);
    write(xml, "links_max", input.links_max);
    write(xml, "mom2_max", input.mom2_max);
    write(xml, "nrow", input.nrow);

    pop(xml);
  }

  //! Propagator input
  void read(XMLReader& xml, const string& path, InlineBuildingBlocksParams::Prop_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "BkwdPropFileName", input.BkwdPropFileName);
    read(inputtop, "BkwdPropG5Format", input.BkwdPropG5Format);
    read(inputtop, "GammaInsertion", input.GammaInsertion);
    read(inputtop, "BBFileNamePattern", input.BBFileNamePattern);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlineBuildingBlocksParams::Prop_t& input)
  {
    push(xml, path);

    write(xml, "BkwdPropFileName", input.BkwdPropFileName);
    write(xml, "BkwdPropG5Format", input.BkwdPropG5Format);
    write(xml, "GammaInsertion", input.GammaInsertion);
    write(xml, "BBFileNamePattern", input.BBFileNamePattern);

    pop(xml);
  }

  //! BB parameters
  void read(XMLReader& xml, const string& path, InlineBuildingBlocksParams::BB_out_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "OutFileName", input.OutFileName);
    read(inputtop, "FrwdPropFileName", input.FrwdPropFileName);
    read(inputtop, "BkwdProps", input.BkwdProps);
  }

  //! BB parameters
  void write(XMLWriter& xml, const string& path, const InlineBuildingBlocksParams::BB_out_t& input)
  {
    push(xml, path);

    write(xml, "OutFileName", input.OutFileName);
    write(xml, "FrwdPropFileName", input.FrwdPropFileName);
    write(xml, "BkwdProps", input.BkwdProps);

    pop(xml);
  }


  // Param stuff
  InlineBuildingBlocksParams::InlineBuildingBlocksParams() { frequency = 0; }

  InlineBuildingBlocksParams::InlineBuildingBlocksParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;

      // Read program parameters
      read(paramtop, "Param", param);

      // Read in the building block setup
      read(paramtop, "BuildingBlocks", bb);

      // Possible alternate XML file pattern
      if (paramtop.count("xml_file") != 0) 
      {
	read(paramtop, "xml_file", xml_file);
      }
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  void
  InlineBuildingBlocksParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    Chroma::write(xml_out, "Param", param);
    Chroma::write(xml_out, "BuildingBlocks", bb);
    QDP::write(xml_out, "xml_file", xml_file);

    pop(xml_out);
  }


  //###################################################################################//
  // Accept All Link Patterns                                                          //
  //###################################################################################//

  void AllLinkPatterns( bool &                          DoThisPattern,
			bool &                          DoFurtherPatterns,
			multi1d< unsigned short int > & LinkPattern )
  {
    DoThisPattern     = true;
    DoFurtherPatterns = true;
    
    return;
  }


  // Function call
  void 
  InlineBuildingBlocks::operator()(const multi1d<LatticeColorMatrix>& u,
				   XMLBufferWriter& gauge_xml,
				   unsigned long update_no,
				   XMLWriter& xml_out) 
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "ExampleBuildingBlocks");
      write(xml_out, "update_no", update_no);
      write(xml_out, "xml_file", xml_file);
      pop(xml_out);

      XMLFileWriter xml(xml_file);
      func(u, gauge_xml, update_no, xml);
    }
    else
    {
      func(u, gauge_xml, update_no, xml_out);
    }
  }


  // Function call
  void 
  InlineBuildingBlocks::func(const multi1d<LatticeColorMatrix>& U,
			     XMLBufferWriter& gauge_xml,
			     unsigned long update_no,
			     XMLWriter& XmlOut) 
  {
    START_CODE();

    push(XmlOut, "ExampleBuildingBlocks");
    write(XmlOut, "update_no", update_no);

    QDPIO::cout << " ExampleBuildingBlocks" << endl;
    QDPIO::cout << "     volume: " << params.param.nrow[0];
    for (int i=1; i<Nd; ++i) {
      QDPIO::cout << " x " << params.param.nrow[i];
    }
    QDPIO::cout << endl;

    //#################################################################################//
    // Echo Arguments                                                                  //
    //#################################################################################//

    // will capture all would-be standard output
    TextWriter Out( params.bb.OutFileName );

    Out <<                                                                      "\n";
    Out << "  NX                                  = " << params.param.nrow[0] << "\n";
    Out << "  NY                                  = " << params.param.nrow[1] << "\n";
    Out << "  NZ                                  = " << params.param.nrow[2] << "\n";
    Out << "  NT                                  = " << params.param.nrow[3] << "\n";
    Out << "  Forward Propagator                  = " << params.bb.FrwdPropFileName          << "\n";
    Out <<                                                                                     "\n";
    for(int loop=0; loop < params.bb.BkwdProps.size(); ++loop)
    {
      const InlineBuildingBlocksParams::Prop_t& prop = params.bb.BkwdProps[loop];
      Out << "  Backward Propagator                 = " << prop.BkwdPropFileName            << "\n";
      Out << "  Backward Propagator Gamma5 Format   = " << prop.BkwdPropG5Format            << "\n";
      Out << "  Gamma Insertion                     = " << prop.GammaInsertion              << "\n";
      Out << "  Building Blocks                     = " << prop.BBFileNamePattern           << "\n";
      Out <<                                                                                   "\n";
    }

    Out << "  Maximum Number of Links             = " << params.param.links_max              << "\n";
    Out << "  Maximum Spatial Momentum Squared    = " << params.param.mom2_max               << "\n"; 

    Out << "  Text Output File Name               = " << params.bb.OutFileName               << "\n";
    Out <<                                                                                     "\n";
    Out.flush();

    //#################################################################################//
    // Record the CVS Info for BuildingBlocks                                          //
    //#################################################################################//

    CVSBuildingBlocks( Out );
    Out << "CVSBuildingBlocks_hh = " << CVSBuildingBlocks_hh << "\n";
    Out <<                                                      "\n";

    //#################################################################################//
    // XML output
    //#################################################################################//

    proginfo(XmlOut);    // Print out basic program info

    push(XmlOut, "Output_version");
    write(XmlOut, "out_version", 2);
    pop(XmlOut);

    //###############################################################################//
    // Read Gauge Field                                                              //
    //###############################################################################//

    Out << "Attempt to initialize the gauge field" << "\n";  Out.flush();

    // Write out the input
    params.write(XmlOut, "Input");

    // Write out the config info
    write(XmlOut, "Config_info", gauge_xml);

    // check that the gauge field seems normal
    Double ave_plaq, ave_spacelike_plaq, ave_timelike_plaq, ave_link_trace;
    MesPlq( U, ave_plaq, ave_spacelike_plaq, ave_timelike_plaq, ave_link_trace );
    Out << "basic gauge field observables"                         << "\n";
    Out << "average plaquette            = " << ave_plaq           << "\n";
    Out << "average space-like plaquette = " << ave_spacelike_plaq << "\n";
    Out << "average time-like plaquette  = " << ave_timelike_plaq  << "\n";
    Out << "average link trace           = " << ave_link_trace     << "\n";

    push(XmlOut, "Observables");
    write(XmlOut, "ave_plaq", ave_plaq);
    write(XmlOut, "ave_spacelike_plaq", ave_spacelike_plaq);
    write(XmlOut, "ave_timelike_plaq", ave_timelike_plaq);
    write(XmlOut, "ave_link_trace", ave_link_trace);
    pop(XmlOut);

    //#################################################################################//
    // Read Forward Propagator                                                         //
    //#################################################################################//

    XMLReader FrwdPropXML, FrwdPropRecordXML;
    LatticePropagator F;
    ChromaProp_t prop_header;
    PropSource_t source_header;
    Out << "reading forward propagator " << params.bb.FrwdPropFileName << " ... " << "\n";  Out.flush();
    {
      Out << "assuming chroma format for forward propagator" << "\n";  Out.flush();
      readQprop( FrwdPropXML, FrwdPropRecordXML, F, params.bb.FrwdPropFileName, QDPIO_SERIAL );
   
      // Try to invert this record XML into a ChromaProp struct
      try
      {
	read(FrwdPropRecordXML, "/Propagator/ForwardProp", prop_header);
	read(FrwdPropRecordXML, "/Propagator/PropSource", source_header);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << "Error extracting forward_prop header: " << e << endl;
	QDP_abort(1);
      }
    }
    Out << "finished reading forward propagator " << params.bb.FrwdPropFileName << "\n";  Out.flush();


    // Derived from input prop
    int  j_decay = source_header.j_decay;

    // Sanity check - write out the norm2 of the forward prop in the j_decay direction
    // Use this for any possible verification
    {
      SftMom phases( 0, true, j_decay );

      multi1d<Double> FrwdPropCheck = sumMulti( localNorm2( F ), phases.getSet() );

      Out << "forward propagator check = " << FrwdPropCheck[0] << "\n";  Out.flush();

      // Write out the forward propagator header
      push(XmlOut, "ForwardProp");
      write(XmlOut, "FrwdPropFileName", params.bb.FrwdPropFileName);
      write(XmlOut, "FrwdPropXML", FrwdPropXML);
      write(XmlOut, "FrwdPropRecordXML", FrwdPropRecordXML);
      write(XmlOut, "FrwdPropCheck", FrwdPropCheck);
      pop(XmlOut);
    }

    //#################################################################################//
    // Read Backward (or Sequential) Propagators                                       //
    //#################################################################################//

    const int NF = params.bb.BkwdProps.size();

    XMLArrayWriter  XmlSeqSrc(XmlOut, NF);
    push(XmlSeqSrc, "SequentialSource");

    for(int loop = 0; loop < NF; ++loop)
    {
      push(XmlSeqSrc);
      write(XmlSeqSrc, "loop_ctr", loop);

      Out << "Loop = " << loop << "\n";  Out.flush();
      QDPIO::cout << "Loop = " << loop << endl;

      XMLReader BkwdPropXML, BkwdPropRecordXML;
      multi1d< LatticePropagator > B( 1 );
      SeqSource_t seqsource_header;
      Out << "reading backward u propagator " << params.bb.BkwdProps[loop].BkwdPropFileName << " ... " << "\n";  Out.flush();
      {
	Out << "assuming chroma format for backward u propagator" << "\n";  Out.flush();
	readQprop( BkwdPropXML, BkwdPropRecordXML, B[0], params.bb.BkwdProps[loop].BkwdPropFileName, QDPIO_SERIAL );

	// Try to invert this record XML into a ChromaProp struct
	// Also pull out the id of this source
	// NEED SECURITY HERE - need a way to cross check props. Use the ID.
	try
	{
	  read(BkwdPropRecordXML, "/SequentialProp/SeqSource", seqsource_header);
	}
	catch (const string& e) 
	{
	  QDPIO::cerr << "Error extracting seqprop header: " << e << endl;
	  QDP_abort(1);
	}
      }
      Out << "finished reading backward u propagator " << params.bb.BkwdProps[loop].BkwdPropFileName << " ... " << "\n";  Out.flush();
      
      // Sanity check - write out the norm2 of the forward prop in the j_decay direction
      // Use this for any possible verification
      {
	SftMom phases( 0, true, j_decay );

	multi1d<Double> BkwdPropCheck = sumMulti( localNorm2( B[0] ), phases.getSet() );

	Out << "backward u propagator check = " << BkwdPropCheck[0] << "\n";  Out.flush();
      
	// Write out the forward propagator header
	push(XmlSeqSrc, "BackwardProp");
	write(XmlSeqSrc, "BkwdPropFileName", params.bb.BkwdProps[loop].BkwdPropFileName);
	write(XmlSeqSrc, "BkwdPropG5Format", params.bb.BkwdProps[loop].BkwdPropG5Format);
	write(XmlSeqSrc, "SequentialSourceType", seqsource_header.seq_src);
	write(XmlSeqSrc, "BkwdPropXML", BkwdPropXML);
	write(XmlSeqSrc, "BkwdPropRecordXML", BkwdPropRecordXML);
	write(XmlSeqSrc, "BkwdPropCheck", BkwdPropCheck);
	pop(XmlSeqSrc);
      }

      //#################################################################################//
      // Additional Gamma Matrix Insertions                                              //
      //#################################################################################//
    
      multi1d<int> GammaInsertions(1);
      GammaInsertions[0] = params.bb.BkwdProps[loop].GammaInsertion;

      if (GammaInsertions[0] < 0 || GammaInsertions[0] >= Ns*Ns)
      {
	QDPIO::cerr << "InlineBuildingBlocks: Gamma insertion out of bounds: " << GammaInsertions[0] << endl;
	QDP_abort(1);
      }

      //#################################################################################//
      // Some output                                                                     //
      //#################################################################################//
    
      QDPIO::cout << "Seqsource name  = " << seqsource_header.seq_src << endl;
      QDPIO::cout << "Gamma insertion = " << params.bb.BkwdProps[loop].GammaInsertion << endl;

      write(XmlSeqSrc, "seq_src", seqsource_header.seq_src);
      write(XmlSeqSrc, "gamma_insertion", GammaInsertions[0]);
      write(XmlSeqSrc, "t_source", source_header.t_source);
      write(XmlSeqSrc, "t_sink", seqsource_header.t_sink);
      write(XmlSeqSrc, "sink_mom", seqsource_header.sink_mom);
	
      //#################################################################################//
      // Convert Backward Propagator Format                                              //
      //#################################################################################//
    
      if( params.bb.BkwdProps[loop].BkwdPropG5Format == "G5_B" )
      {
	LatticePropagator Bu = B[0];
	B[0] = Gamma( 15 ) * Bu;
      }
      else if( params.bb.BkwdProps[loop].BkwdPropG5Format == "B_G5" )
      {
	LatticePropagator Bu = B[0];
	B[0] = Bu * Gamma( 15 );
      }
      else if( params.bb.BkwdProps[loop].BkwdPropG5Format == "G5_B_G5" )
      {
	LatticePropagator Bu = B[0];
	B[0] = Gamma( 15 ) * Bu * Gamma( 15 );
      }

      //#################################################################################//
      // Set Momenta                                                                     //
      //#################################################################################//

      multi1d< int > SnkMom( Nd - 1 );
      if (params.param.use_sink_offset)
	SnkMom = seqsource_header.sink_mom;
      else
	SnkMom = 0;

      SftMom Phases( params.param.mom2_max, SnkMom, false, j_decay );

      //#################################################################################//
      // Construct File Names                                                            //
      //#################################################################################//

      int NQ = Phases.numMom();

      multi2d< string > Files( 1, NQ );

      const int BBFileNameLength = params.bb.BkwdProps[loop].BBFileNamePattern.length() + 3 * 3 + 1;

      for( int q = 0; q < NQ; q ++ )
      {
	multi1d< int > Q = Phases.numToMom( q );

	char XSign = '+';
	char YSign = '+';
	char ZSign = '+';

	if( Q[0] < 0 )
	{
	  XSign = '-';
	}
	if( Q[1] < 0 )
	{
	  YSign = '-';
	}
	if( Q[2] < 0 )
	{
	  ZSign = '-';
	}

	char* bbf = new char[BBFileNameLength + 1];
	sprintf( bbf, params.bb.BkwdProps[loop].BBFileNamePattern.c_str(), ZSign, abs(Q[2]), YSign, abs(Q[1]), XSign, abs(Q[0]) );

	Files(0,q) = bbf;

	delete[] bbf;
      }

      //#################################################################################//
      // Construct Building Blocks                                                       //
      //#################################################################################//
    
      Out << "calculating building blocks" << "\n";  Out.flush();
      QDPIO::cout << "calculating building blocks" << endl;

      multi1d< int > Flavors( 1 );
      Flavors[0] = loop;                    // Dru puts a Flavor into the BB

      const signed short int T1 = 0;
      const signed short int T2 = params.param.nrow[j_decay] - 1;
      const signed short int DecayDir = j_decay;

      BuildingBlocks( B, F, U, GammaInsertions, Flavors,
		      params.param.links_max, AllLinkPatterns, Phases, Files, T1, T2,
		      seqsource_header.seq_src, seqsource_header.sink_mom, DecayDir);

      Out << "finished calculating building blocks for loop = " << loop << "\n";  Out.flush();
      QDPIO::cout << "finished calculating building blocks for loop = " << loop << endl;

      pop(XmlSeqSrc);   // elem
    } // end loop over sequential sources

    pop(XmlSeqSrc);  // SequentialSource

    pop(XmlOut);   // ExampleBuildingBlocks

    Out << "\n" << "FINISHED" << "\n" << "\n";
    Out.close();

    END_CODE();
  } 

};
