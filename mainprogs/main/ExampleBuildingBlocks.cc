//###################################################################################//
//###################################################################################//
//                                                                                   //
// ExampleBuildingBlocks.cc                                                          //
//                                                                                   //
//###################################################################################//
//###################################################################################//
//                                                                                   //
// description:                                                                      //
//                                                                                   //
// Executing ExampleBuildingBlocks without any arguments (or an incorrect number of  //
// arguments) will result in the required argument list being printed to stdout.     //
//                                                                                   //
// Note that the file name patterns for the u and d building blocks must be a string //
// of the form "...%c%i...%c%i...%c%i..." where the first %c%i corresponds to qz,    //
// the second to qy, and the third to qx.                                            //
//                                                                                   //
// history:                                                                          //
//                                                                                   //
// There were at least four versions of "MIT" code.  Andrew Pochinsky has a c        //
// version. Dmitri Dolgov has a c++ version.  Dru B. Renner has c and c++ versions.  //
// All were independent and checked against one another.  Of course, all were        //
// developed under the guidance of John W. Negele.  The code here is just the        //
// "Building Blocks" portion of the MIT code.                                        //
//                                                                                   //
// authors:                                                                          //
//                                                                                   //
// Dru B. Renner, dru@mit.edu, 2002 - port of Building Blocks (MIT) code to qdp++    //
//                                                                                   //
// There are others who have contributed since the code has been migrated to qdp++.  //
// The cvs log entries indicate these other authors.                                 //
//                                                                                   //
//###################################################################################//
//###################################################################################//

#include <iostream>
#include <string.h>
#include <assert.h>
#include "chroma.h"

using namespace Chroma;

//###################################################################################//
// CVS Header                                                                        //
//###################################################################################//

static const char* const CVSExampleBuildingBlocks_hh =
"$Header: /home/bjoo/fromJLAB/cvsroot/chroma_base/mainprogs/main/ExampleBuildingBlocks.cc,v 1.20 2005-03-12 18:42:34 edwards Exp $";


/*
 * Input 
 */
//! Parameters for running program
struct Param_t
{
  int mom2_max;            // (mom)^2 <= mom2_max
  int links_max;           // maximum number of links
  multi1d<int> nrow;       // lattice size
};

//! Propagators
struct Prop_t
{
  string   BkwdPropFileName;   // backward propagator
  string   BkwdPropG5Format;   // backward propagators Gamma5 Format
  int      GammaInsertion;     // second gamma insertion
  string   BBFileNamePattern;  // BB output file name pattern
};

//! BB output
struct BB_out_t
{
  string            OutFileName;
  string            FrwdPropFileName;   // input forward prop
  multi1d<Prop_t>   BkwdProps;
};

//! Mega-structure of all input
struct BB_input_t
{
  Param_t          param;
  Cfg_t            cfg;
  BB_out_t         bb;
};


//! Propagator parameters
void read(XMLReader& xml, const string& path, Prop_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "BkwdPropFileName", input.BkwdPropFileName);
  read(inputtop, "BkwdPropG5Format", input.BkwdPropG5Format);
  read(inputtop, "GammaInsertion", input.GammaInsertion);
  read(inputtop, "BBFileNamePattern", input.BBFileNamePattern);
}



//! BB parameters
void read(XMLReader& xml, const string& path, BB_out_t& input)
{
  XMLReader inputtop(xml, path);

  read(inputtop, "OutFileName", input.OutFileName);
  read(inputtop, "FrwdPropFileName", input.FrwdPropFileName);
  read(inputtop, "BkwdProps", input.BkwdProps);
}


// Reader for input parameters
void read(XMLReader& xml, const string& path, Param_t& param)
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

  read(paramtop, "links_max", param.links_max);
  read(paramtop, "mom2_max", param.mom2_max);
  read(paramtop, "nrow", param.nrow);
}


// Reader for input parameters
void read(XMLReader& xml, const string& path, BB_input_t& input)
{
  XMLReader inputtop(xml, path);

  // Read all the input groups
  try
  {
    // Read program parameters
    read(inputtop, "Param", input.param);

    // Read in the gauge configuration info
    read(inputtop, "Cfg", input.cfg);

    // Read in the building block setup
    read(inputtop, "BuildingBlocks", input.bb);
  }
  catch (const string& e) 
  {
    QDPIO::cerr << "Error reading ExampleBuildingBlocks data: " << e << endl;
    QDP_abort(1);
  }
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

//###################################################################################//
// Main Function                                                                     //
//###################################################################################//

int main( int argc, char** argv )
{
  //#################################################################################//
  // Initialize Chroma                                                               //
  //#################################################################################//

  Chroma::initialize(&argc, &argv);

  START_CODE();

  //#################################################################################//
  // Read input parameters                                                           //
  //#################################################################################//

  // Input parameter structure
  BB_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader XmlIn(Chroma::getXMLInputFileName());

  // Read data
  read(XmlIn, "/ExampleBuildingBlocks", input);

  QDPIO::cout << " ExampleBuildingBlocks" << endl;
  QDPIO::cout << "     volume: " << input.param.nrow[0];
  for (int i=1; i<Nd; ++i) {
    QDPIO::cout << " x " << input.param.nrow[i];
  }
  QDPIO::cout << endl;

  //#################################################################################//
  // Echo Arguments                                                                  //
  //#################################################################################//

  // will capture all would-be standard output
  TextWriter Out( input.bb.OutFileName );

  Out <<                                                                      "\n";
  Out << "  NX                                  = " << input.param.nrow[0] << "\n";
  Out << "  NY                                  = " << input.param.nrow[1] << "\n";
  Out << "  NZ                                  = " << input.param.nrow[2] << "\n";
  Out << "  NT                                  = " << input.param.nrow[3] << "\n";
  Out << "  Gauge Field                         = " << input.cfg.cfg_file  << "\n";
  Out << "  Forward Propagator                  = " << input.bb.FrwdPropFileName          << "\n";
  Out <<                                                                                     "\n";
  for(int loop=0; loop < input.bb.BkwdProps.size(); ++loop)
  {
    const Prop_t& prop = input.bb.BkwdProps[loop];
    Out << "  Backward Propagator                 = " << prop.BkwdPropFileName            << "\n";
    Out << "  Backward Propagator Gamma5 Format   = " << prop.BkwdPropG5Format            << "\n";
    Out << "  Gamma Insertion                     = " << prop.GammaInsertion              << "\n";
    Out << "  Building Blocks                     = " << prop.BBFileNamePattern           << "\n";
    Out <<                                                                                   "\n";
  }

  Out << "  Maximum Number of Links             = " << input.param.links_max              << "\n";
  Out << "  Maximum Spatial Momentum Squared    = " << input.param.mom2_max               << "\n"; 

  Out << "  Text Output File Name               = " << input.bb.OutFileName               << "\n";
  Out << "  XML Output File Name                = " << Chroma::getXMLOutputFileName()     << "\n";
  Out <<                                                                                     "\n";
  Out.flush();

  //#################################################################################//
  // Record the CVS Info for BuildingBlocks                                          //
  //#################################################################################//

  CVSBuildingBlocks( Out );
  Out << "CVSExampleBuildingBlocks_hh = " << CVSExampleBuildingBlocks_hh << "\n";
  Out <<                                                                    "\n";

  //#################################################################################//
  // Set Lattice Geometry                                                            //
  //#################################################################################//

  Layout::setLattSize( input.param.nrow );
  Layout::create();

  //#################################################################################//
  // XML output
  //#################################################################################//

  // capture XML output
  XMLFileWriter& XmlOut = Chroma::getXMLOutputInstance();
  push(XmlOut, "ExampleBuildingBlocks");

  proginfo(XmlOut);    // Print out basic program info

  push(XmlOut, "Output_version");
  write(XmlOut, "out_version", 2);
  pop(XmlOut);

  //###############################################################################//
  // Read Gauge Field                                                              //
  //###############################################################################//

  Out << "Attempt to initialize the gauge field" << "\n";  Out.flush();

  multi1d<LatticeColorMatrix> U(Nd);
  {
    XMLReader GaugeFileXml, GaugeXml;

    Out << "reading gauge field " << input.cfg.cfg_file << " ... " << "\n";  Out.flush();
    gaugeStartup(GaugeFileXml, GaugeXml, U, input.cfg);
    Out << "finished reading gauge field " << input.cfg.cfg_file << "\n";  Out.flush();

    // Write out the input
    write(XmlOut, "Input", XmlIn);

    // Write out the config info
    write(XmlOut, "Config_info", GaugeXml);
  }

  // check that the gauge field seems normal
  unitarityCheck( U );
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
  XmlOut.flush();

  //#################################################################################//
  // Read Forward Propagator                                                         //
  //#################################################################################//

  XMLReader FrwdPropXML, FrwdPropRecordXML;
  LatticePropagator F;
  ChromaProp_t prop_header;
  PropSource_t source_header;
  Out << "reading forward propagator " << input.bb.FrwdPropFileName << " ... " << "\n";  Out.flush();
  {
    Out << "assuming chroma format for forward propagator" << "\n";  Out.flush();
    readQprop( FrwdPropXML, FrwdPropRecordXML, F, input.bb.FrwdPropFileName, QDPIO_SERIAL );
   
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
  Out << "finished reading forward propagator " << input.bb.FrwdPropFileName << "\n";  Out.flush();


  // Derived from input prop
  int  j_decay = source_header.j_decay;
  multi1d<int> t_source = source_header.t_source;

  // Sanity check - write out the norm2 of the forward prop in the j_decay direction
  // Use this for any possible verification
  {
    SftMom phases( 0, true, j_decay );

    multi1d<Double> FrwdPropCheck = sumMulti( localNorm2( F ), phases.getSet() );

    Out << "forward propagator check = " << FrwdPropCheck[0] << "\n";  Out.flush();

    // Write out the forward propagator header
    push(XmlOut, "ForwardProp");
    write(XmlOut, "FrwdPropFileName", input.bb.FrwdPropFileName);
    write(XmlOut, "FrwdPropXML", FrwdPropXML);
    write(XmlOut, "FrwdPropRecordXML", FrwdPropRecordXML);
    write(XmlOut, "FrwdPropCheck", FrwdPropCheck);
    pop(XmlOut);
    XmlOut.flush();
  }

  //#################################################################################//
  // Read Backward (or Sequential) Propagators                                       //
  //#################################################################################//

  const int NF = input.bb.BkwdProps.size();

  multi1d< LatticePropagator > B( NF );
  multi1d< int > GammaInsertions( NF );

  XMLArrayWriter  XmlSeqSrc(XmlOut, NF);
  push(XmlSeqSrc, "SequentialSource");

  for(int loop = 0; loop < NF; ++loop)
  {
    push(XmlSeqSrc);
    write(XmlSeqSrc, "loop_ctr", loop);

    XMLReader BkwdPropXML, BkwdPropRecordXML;
    LatticePropagator Bu;
    SeqSource_t seqsource_header;
    Out << "reading backward u propagator " << input.bb.BkwdProps[loop].BkwdPropFileName << " ... " << "\n";  Out.flush();
    {
      Out << "assuming chroma format for backward u propagator" << "\n";  Out.flush();
      readQprop( BkwdPropXML, BkwdPropRecordXML, Bu, input.bb.BkwdProps[loop].BkwdPropFileName, QDPIO_SERIAL );

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
    Out << "finished reading backward u propagator " << input.bb.BkwdProps[loop].BkwdPropFileName << " ... " << "\n";  Out.flush();
      
    // Sanity check - write out the norm2 of the forward prop in the j_decay direction
    // Use this for any possible verification
    {
      SftMom phases( 0, true, j_decay );

      multi1d<Double> BkwdPropCheck = sumMulti( localNorm2( Bu ), phases.getSet() );

      Out << "backward u propagator check = " << BkwdPropCheck[0] << "\n";  Out.flush();
      
      // Write out the forward propagator header
      push(XmlSeqSrc, "BackwardProp");
      write(XmlSeqSrc, "BkwdPropFileName", input.bb.BkwdProps[loop].BkwdPropFileName);
      write(XmlSeqSrc, "BkwdPropG5Format", input.bb.BkwdProps[loop].BkwdPropG5Format);
      write(XmlSeqSrc, "SequentialSourceType", seqsource_header.seq_src);
      write(XmlSeqSrc, "BkwdPropXML", BkwdPropXML);
      write(XmlSeqSrc, "BkwdPropRecordXML", BkwdPropRecordXML);
      write(XmlSeqSrc, "BkwdPropCheck", BkwdPropCheck);
      pop(XmlSeqSrc);
      XmlOut.flush();
    }

    // Derived from input seqprop
    string seq_src = seqsource_header.seq_src;
    QDPIO::cout << "Seqsource name = " << seqsource_header.seq_src << endl;
    int           t_sink   = seqsource_header.t_sink;
    multi1d<int>  sink_mom = seqsource_header.sink_mom;

    write(XmlSeqSrc, "seq_src", seq_src);
    write(XmlSeqSrc, "t_source", t_source);
    write(XmlSeqSrc, "t_sink", t_sink);
    write(XmlSeqSrc, "sink_mom", sink_mom);
	
    //#################################################################################//
    // Additional Gamma Matrix Insertions                                              //
    //#################################################################################//
    
    GammaInsertions[loop] = input.bb.BkwdProps[loop].GammaInsertion;

    if (GammaInsertions[loop] < 0 || GammaInsertions[loop] >= Ns*Ns)
    {
      QDPIO::cerr << argv[0] << ": Gamma insertion out of bounds: " << GammaInsertions[loop] << endl;
      QDP_abort(1);
    }

    //#################################################################################//
    // Convert Backward Propagator Format                                              //
    //#################################################################################//
    
    if( input.bb.BkwdProps[loop].BkwdPropG5Format == "G5_B" )
    {
      B[loop] = Gamma( 15 ) * Bu;
    }
    else if( input.bb.BkwdProps[loop].BkwdPropG5Format == "B_G5" )
    {
      B[loop] = Bu * Gamma( 15 );
    }
    else if( input.bb.BkwdProps[loop].BkwdPropG5Format == "G5_B_G5" )
    {
      B[loop] = Gamma( 15 ) * Bu * Gamma( 15 );
    }
    else
    {
      B[loop] = Bu;
    }

    pop(XmlSeqSrc);   // elem
  } // end loop over sequential sources

  pop(XmlSeqSrc);  // SequentialSource

  //#################################################################################//
  // Set Momenta                                                                     //
  //#################################################################################//

  // WARNING: DRU HAD SnkMom HARDWIRED TO ZERO. I THINK IT SHOULD FLOAT
  // THIS CHANGES THE INSERTION MOMENTA!
  multi1d< int > SnkMom( Nd - 1 );
  SnkMom = 0;

//  SftMom Phases( input.param.mom2_max, input.param.sink_mom, false, j_decay );
  SftMom Phases( input.param.mom2_max, SnkMom, false, j_decay );

  //#################################################################################//
  // Construct File Names                                                            //
  //#################################################################################//

  int NQ = Phases.numMom();

  multi2d< string > Files( NF, NQ );

  for(int loop = 0; loop < NF; ++loop)
  {
    const int BBFileNameLength = input.bb.BkwdProps[loop].BBFileNamePattern.length() + 3 * 3 + 1;

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
      sprintf( bbf, input.bb.BkwdProps[loop].BBFileNamePattern.c_str(), ZSign, abs(Q[2]), YSign, abs(Q[1]), XSign, abs(Q[0]) );

      Files(loop,q) = bbf;

      delete[] bbf;
    }
  }

  //#################################################################################//
  // Construct Building Blocks                                                       //
  //#################################################################################//

  Out << "calculating building blocks ..." << "\n";  Out.flush();
  QDPIO::cout << "calculating building blocks ..." << endl;

  const signed short int T1 = 0;
  const signed short int T2 = input.param.nrow[j_decay] - 1;

  BuildingBlocks( B, F, U, GammaInsertions,
		  input.param.links_max, AllLinkPatterns, Phases, Files, T1, T2 );

  Out << "finished calculating building blocks" << "\n";  Out.flush();
  QDPIO::cout << "finished calculating building blocks" << endl;

  //#################################################################################//
  // Close QDP                                                                       //
  //#################################################################################//

  XmlIn.close();

  pop(XmlOut);   // ExampleBuildingBlocks
  XmlOut.close();

  Out << "\n" << "FINISHED" << "\n" << "\n";
  Out.close();

  END_CODE();

  // Time to bolt
  Chroma::finalize();

  exit( 0 );
}

//###################################################################################//
//###################################################################################//
