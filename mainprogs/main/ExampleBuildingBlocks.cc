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

using namespace QDP;

//###################################################################################//
// CVS Header                                                                        //
//###################################################################################//

static const char* const CVSExampleBuildingBlocks_hh =
  "$Header: /home/bjoo/fromJLAB/cvsroot/chroma_base/mainprogs/main/ExampleBuildingBlocks.cc,v 1.16 2004-11-06 20:04:09 edwards Exp $";

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
  // Initialize QDP                                                                  //
  //#################################################################################//

  QDP_initialize( & argc, & argv );

  START_CODE();

  //#################################################################################//
  // Check Arguments                                                                 //
  //#################################################################################//

  if( argc != 21 )
  {
    QDPIO::cout << ""                                                                                                  << endl;
    QDPIO::cout << "Arguments:"                                                                                        << endl;
    QDPIO::cout << ""                                                                                                  << endl;
    QDPIO::cout << "(1)  NX                                  - \"1\", \"2\", ..."                                      << endl;
    QDPIO::cout << "(2)  NY                                  - \"1\", \"2\", ..."                                      << endl;
    QDPIO::cout << "(3)  NZ                                  - \"1\", \"2\", ..."                                      << endl;
    QDPIO::cout << "(4)  NT                                  - \"1\", \"2\", ..."                                      << endl;
    QDPIO::cout << "(5)  Gauge Field                         - Input File Name"                                        << endl;
    QDPIO::cout << "(6)  Gauge Field Format                  - \"archive\" or \"szin\""                                << endl;
    QDPIO::cout << "(7)  Forward Propagator                  - Input File Name"                                        << endl;
    QDPIO::cout << "(8)  Forward Propagator Format           - \"chroma\" or \"szin\""                                 << endl;
    QDPIO::cout << "(9)  Proton Backward U Propagator        - Input File Name (\"NULL\" if u not available)"          << endl;
    QDPIO::cout << "(10) Backward U Propagator Format        - \"chroma\" or \"szin\""                                 << endl;
    QDPIO::cout << "(11) Backward U Propagator Gamma5 Format - \"B\", \"G5_B\", \"B_G5\", or \"G5_B_G5\""              << endl;
    QDPIO::cout << "(12) Proton Backward D Propagator        - Input File Name (\"NULL\" if d not available)"          << endl;
    QDPIO::cout << "(13) Backward D Propagator Format        - \"chroma\" or \"szin\""                                 << endl;
    QDPIO::cout << "(14) Backward D Propagator Gamma5 Format - \"B\", \"G5_B\", \"B_G5\", or \"G5_B_G5\""              << endl;
    QDPIO::cout << "(15) Maximum Number of Links             - \"0\", \"1\", ..."                                      << endl;
    QDPIO::cout << "(16) Maximum Spatial Momentum Squared    - \"0\", \"1\", ..."                                      << endl;
    QDPIO::cout << "(17) Proton U Building Blocks            - Output File Name Pattern (\"NULL\" if u not available)" << endl;
    QDPIO::cout << "(18) Proton D Building Blocks            - Output File Name Pattern (\"NULL\" if d not available)" << endl;
    QDPIO::cout << "(19) Text Output File Name               - Output File Name"                                       << endl;
    QDPIO::cout << "(20) XML Output File Name                - XML Output File Name"                                   << endl;
    QDPIO::cout << ""                                                                                                  << endl;

    exit( 1 );
  }

  //#################################################################################//
  // Capture Arguments                                                               //
  //#################################################################################//

  const unsigned short int NX            = atoi( argv[1] );
  const unsigned short int NY            = atoi( argv[2] );
  const unsigned short int NZ            = atoi( argv[3] );
  const unsigned short int NT            = atoi( argv[4] );
  const char* const GaugeFieldFileName   = argv[5];
  const char* const GaugeFieldFormat     = argv[6];
  const char* const FrwdPropFileName     = argv[7];
  const char* const FrwdPropFormat       = argv[8];
  const char* const BkwdUPropFileName    = argv[9];
  const char* const BkwdUPropFormat      = argv[10];
  const char* const BkwdUPropG5Format    = argv[11];
  const char* const BkwdDPropFileName    = argv[12];
  const char* const BkwdDPropFormat      = argv[13];
  const char* const BkwdDPropG5Format    = argv[14];
  const unsigned short int MaxNLinks     = atoi( argv[15] );
  const unsigned short int MaxQSquared   = atoi( argv[16] );
  const char* const UBBFileNamePattern   = argv[17];
  const char* const DBBFileNamePattern   = argv[18];
  const char* const OutFileName          = argv[19];
  const char* const XMLFileName          = argv[20];

  //#################################################################################//
  // Echo Arguments                                                                  //
  //#################################################################################//

  // will capture all would-be standard output
  TextWriter Out( OutFileName );

  Out <<                                                                        "\n";
  Out << "(1)  NX                                  = " << NX                 << "\n";
  Out << "(2)  NY                                  = " << NY                 << "\n";
  Out << "(3)  NZ                                  = " << NZ                 << "\n";
  Out << "(4)  NT                                  = " << NT                 << "\n";
  Out << "(5)  Gauge Field                         = " << GaugeFieldFileName << "\n";
  Out << "(6)  Gauge Field Format                  = " << GaugeFieldFormat   << "\n";
  Out << "(7)  Forward Propagator                  = " << FrwdPropFileName   << "\n";
  Out << "(8)  Forward Propagator Format           = " << FrwdPropFormat     << "\n";
  Out << "(9)  Proton Backward U Propagator        = " << BkwdUPropFileName  << "\n";
  Out << "(10) Backward U Propagator Format        = " << BkwdUPropFormat    << "\n";
  Out << "(11) Backward U Propagator Gamma5 Format = " << BkwdUPropG5Format  << "\n";
  Out << "(12) Proton Backward D Propagator        = " << BkwdDPropFileName  << "\n";
  Out << "(13) Backward D Propagator Format        = " << BkwdDPropFormat    << "\n";
  Out << "(14) Backward D Propagator Gamma5 Format = " << BkwdDPropG5Format  << "\n";
  Out << "(15) Maximum Number of Links             = " << MaxNLinks          << "\n";
  Out << "(16) Maximum Spatial Momentum Squared    = " << MaxQSquared        << "\n";
  Out << "(17) Proton U Building Blocks            = " << UBBFileNamePattern << "\n";
  Out << "(18) Proton D Building Blocks            = " << DBBFileNamePattern << "\n";
  Out << "(19) Text Output File Name               = " << OutFileName        << "\n";
  Out << "(20) XML Output File Name                = " << XMLFileName        << "\n";
  Out <<                                                                        "\n";
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

  multi1d< int > N ( Nd );

  N[0] = NX;
  N[1] = NY;
  N[2] = NZ;
  N[3] = NT;

  Layout::setLattSize( N );
  Layout::create();

  //#################################################################################//
  // XML output
  //#################################################################################//

  // capture XML output
  XMLFileWriter Xml( XMLFileName );
  push(Xml, "ExampleBuildingBlocks");

  proginfo(Xml);    // Print out basic program info

  push(Xml, "Input");
  write(Xml, "N", N);
  write(Xml, "GaugeFieldFileName", GaugeFieldFileName);
  write(Xml, "GaugeFieldFormat", GaugeFieldFormat);
  write(Xml, "FrwdPropFileName", FrwdPropFileName);
  write(Xml, "FrwdPropFormat", FrwdPropFormat);
  write(Xml, "BkwdUPropFileName", BkwdUPropFileName);
  write(Xml, "BkwdUPropFormat", BkwdUPropFormat);
  write(Xml, "BkwdUPropG5Format", BkwdUPropG5Format);
  write(Xml, "BkwdDPropFileName", BkwdDPropFileName);
  write(Xml, "BkwdDPropFormat", BkwdDPropFormat);
  write(Xml, "MaxNLinks", MaxNLinks);
  write(Xml, "MaxQSquared", MaxQSquared);
  write(Xml, "UBBFileNamePattern", UBBFileNamePattern);
  write(Xml, "DBBFileNamePattern", DBBFileNamePattern);
  write(Xml, "OutFileName", OutFileName);
  write(Xml, "XMLFileName", XMLFileName);
  pop(Xml);

  push(Xml, "Output_version");
  write(Xml, "out_version", 1);
  pop(Xml);

  //###############################################################################//
  // Read Gauge Field                                                              //
  //###############################################################################//

  XMLReader GaugeFieldXML;
  multi1d< LatticeColorMatrix >* U = NULL;
  U = new multi1d< LatticeColorMatrix >( Nd );
  Out << "reading gauge field " << GaugeFieldFileName << " ... " << "\n";  Out.flush();
  if( strcmp( GaugeFieldFormat, "archive" ) == 0 )
  {
    Out << "assuming archive format for gauge field" << "\n";  Out.flush();
    readArchiv( GaugeFieldXML, *U, GaugeFieldFileName );
  }
  if( strcmp( GaugeFieldFormat, "szin" ) == 0 )
  {
    Out << "assuming szin format for gauge field" << "\n";  Out.flush();
    readSzin( GaugeFieldXML, *U, GaugeFieldFileName );
  }
  Out << "finished reading gauge field " << GaugeFieldFileName << "\n";  Out.flush();

  // check that the gauge field seems normal
  unitarityCheck( *U );
  Double ave_plaq, ave_spacelike_plaq, ave_timelike_plaq, ave_link_trace;
  MesPlq( *U, ave_plaq, ave_spacelike_plaq, ave_timelike_plaq, ave_link_trace );
  Out << "basic gauge field observables"                         << "\n";
  Out << "average plaquette            = " << ave_plaq           << "\n";
  Out << "average space-like plaquette = " << ave_spacelike_plaq << "\n";
  Out << "average time-like plaquette  = " << ave_timelike_plaq  << "\n";
  Out << "average link trace           = " << ave_link_trace     << "\n";

  // Write out the config header
  write(Xml, "Config_info", GaugeFieldXML);

  //#################################################################################//
  // Read Forward Propagator                                                         //
  //#################################################################################//

  XMLReader FrwdPropXML, FrwdPropRecordXML;
  LatticePropagator* F = NULL;
  F = new LatticePropagator;
  Out << "reading forward propagator " << FrwdPropFileName << " ... " << "\n";  Out.flush();
  if( strcmp( FrwdPropFormat, "chroma" ) == 0 )
  {
    Out << "assuming chroma format for forward propagator" << "\n";  Out.flush();
    readQprop( FrwdPropXML, FrwdPropRecordXML, *F, FrwdPropFileName, QDPIO_SERIAL );
  }
  if( strcmp( FrwdPropFormat, "szin" ) == 0 )
  {
    Out << "assuming szin format for forward propagator" << "\n";  Out.flush();
    readSzinQprop( FrwdPropXML, *F, FrwdPropFileName );
  }
  Out << "finished reading forward propagator " << FrwdPropFileName << "\n";  Out.flush();

  {
    SftMom phases( 0, true, Nd-1 );

    multi1d<Double> FrwdPropCheck = sumMulti( localNorm2( *F ), phases.getSet() );

    Out << "forward propagator check = " << FrwdPropCheck[0] << "\n";  Out.flush();
  }

  // Write out the forward propagator header
  write(Xml, "FrwdPropXML", FrwdPropXML);
  write(Xml, "FrwdPropRecordXML", FrwdPropRecordXML);

  Xml.flush();

  //#################################################################################//
  // Read Backward (or Sequential) U Propagator                                      //
  //#################################################################################//

  int NF = 0;

  int HasU = 0;
  XMLReader BkwdUPropXML, BkwdUPropRecordXML;
  LatticePropagator* Bu = NULL;
  Out << "reading backward u propagator " << BkwdUPropFileName << " ... " << "\n";  Out.flush();
  if( strcmp( BkwdUPropFileName, "NULL" ) != 0 )
  {
    Bu = new LatticePropagator;

    if( strcmp( BkwdUPropFormat, "chroma" ) == 0 )
    {
      Out << "assuming chroma format for backward u propagator" << "\n";  Out.flush();
      readQprop( BkwdUPropXML, BkwdUPropRecordXML, *Bu, BkwdUPropFileName, QDPIO_SERIAL );
    }
    if( strcmp( BkwdUPropFormat, "szin" ) == 0 )
    {
      Out << "assuming szin format for backward u propagator" << "\n";  Out.flush();
      readSzinQprop( BkwdUPropXML, *Bu, BkwdUPropFileName );
    }

    HasU = 1;
    NF ++;
  }
  Out << "finished reading backward u propagator " << BkwdUPropFileName << " ... " << "\n";  Out.flush();

  if ( HasU == 1 )
  {
    SftMom phases( 0, true, Nd-1 );

    multi1d<Double> BkwdUPropCheck = sumMulti( localNorm2( *Bu ), phases.getSet() );

    Out << "backward u propagator check = " << BkwdUPropCheck[0] << "\n";  Out.flush();
  }

  // Write out the forward propagator header
  write(Xml, "BkwdUPropXML", BkwdUPropXML);
  write(Xml, "BkwdUPropRecordXML", BkwdUPropRecordXML);

  //#################################################################################//
  // Read Backward (or Sequential) D Propagator                                      //
  //#################################################################################//

  int HasD = 0;
  XMLReader BkwdDPropXML, BkwdDPropRecordXML;
  LatticePropagator* Bd = NULL;
  Out << "reading backward d propagator " << BkwdDPropFileName << " ... " << "\n";  Out.flush();
  if( strcmp( BkwdDPropFileName, "NULL" ) != 0 )
  {
    Bd = new LatticePropagator;

    if( strcmp( BkwdDPropFormat, "chroma" ) == 0 )
    {
      Out << "assuming chroma format for backward d propagator" << "\n";  Out.flush();
      readQprop( BkwdDPropXML, BkwdDPropRecordXML, *Bd, BkwdDPropFileName, QDPIO_SERIAL );
    }
    if( strcmp( BkwdDPropFormat, "szin" ) == 0 )
    {
      Out << "assuming szin format for backward d propagator" << "\n";  Out.flush();
      readSzinQprop( BkwdDPropXML, *Bd, BkwdDPropFileName );
    }

    HasD = 1;
    NF ++;
  }
  Out << "finished reading backward d propagator " << BkwdDPropFileName << " ... " << "\n";  Out.flush();
 
  if ( HasD == 1 )
  {
    SftMom phases( 0, true, Nd-1 );

    multi1d<Double> BkwdDPropCheck = sumMulti( localNorm2( *Bd ), phases.getSet() );

    Out << "backward d propagator check = " << BkwdDPropCheck[0] << "\n";  Out.flush();
  }

  // Write out the forward propagator header
  write(Xml, "BkwdDPropXML", BkwdDPropXML);
  write(Xml, "BkwdDPropRecordXML", BkwdDPropRecordXML);

  //#################################################################################//
  // Basic Information                                                               //
  //#################################################################################//

  Out << "NF = " << NF << "\n";
  if( HasU == 1 )
  {
    Out << "u flavor available" << "\n";
  }
  else
  {
    Out << "u flavor not available" << "\n";
  }
  if( HasD == 1 )
  {
    Out << "d flavor available" << "\n";
  }
  else
  {
    Out << "d flavor not available" << "\n";
  }
  Out.flush();

  if( NF != 2 )
  {
    Out << "NF must be 2 due to simple assumption in BuildingBlocks.cc when writing Flavor entry in footer" << "\n";
    QDP_finalize();
  }

  //#################################################################################//
  // Convert Backward Propagator Format                                              //
  //#################################################################################//

  if( HasU == 1 )
  {
  if( strcmp( BkwdUPropG5Format, "G5_B" ) == 0 )
  {
    *Bu = Gamma( 15 ) * *Bu;
  }
  if( strcmp( BkwdUPropG5Format, "B_G5" ) == 0 )
  {
    *Bu = *Bu * Gamma( 15 );
  }
  if( strcmp( BkwdUPropG5Format, "G5_B_G5" ) == 0 )
  {
    *Bu = Gamma( 15 ) * *Bu * Gamma( 15 );
  }
  }

  if( HasD == 1 )
  {
  if( strcmp( BkwdDPropG5Format, "G5_B" ) == 0 )
  {
    *Bd = Gamma( 15 ) * *Bd;
  }
  if( strcmp( BkwdDPropG5Format, "B_G5" ) == 0 )
  {
    *Bd = *Bd * Gamma( 15 );
  }
  if( strcmp( BkwdDPropG5Format, "G5_B_G5" ) == 0 )
  {
    *Bd = Gamma( 15 ) * *Bd * Gamma( 15 );
  }
  }

  multi1d< LatticePropagator >* B = NULL;
  B = new multi1d< LatticePropagator >( NF );
  int f = 0;
  if( HasU == 1 )
  {
    (*B)[ f ] = *Bu;
    assert( Bu != NULL );
    delete Bu;  Bu = NULL;
    f ++;
  }
  if( HasD == 1 )
  {
    (*B)[ f ] = *Bd;
    assert( Bd != NULL );
    delete Bd;  Bd = NULL;
    f ++;
  }

  //#################################################################################//
  // Set Momenta                                                                     //
  //#################################################################################//

  multi1d< int > SnkMom( Nd - 1 );
  SnkMom = 0;

  SftMom Phases( MaxQSquared, SnkMom, false, Nd - 1 );

  //#################################################################################//
  // Construct File Names                                                            //
  //#################################################################################//

  int NQ = Phases.numMom();

  multi1d< char* > UBBFileNames( NQ );
  multi1d< char* > DBBFileNames( NQ );

  const int UBBFileNameLength = strlen( UBBFileNamePattern ) + 3 * 3 + 1;
  const int DBBFileNameLength = strlen( DBBFileNamePattern ) + 3 * 3 + 1;

  for( int q = 0; q < NQ; q ++ )
  {
    UBBFileNames[ q ] = (char*) malloc( UBBFileNameLength * sizeof( char ) );
    DBBFileNames[ q ] = (char*) malloc( DBBFileNameLength * sizeof( char ) );
    assert( UBBFileNames[ q ] != NULL );
    assert( DBBFileNames[ q ] != NULL );

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

    sprintf( UBBFileNames[q], UBBFileNamePattern, ZSign, abs(Q[2]), YSign, abs(Q[1]), XSign, abs(Q[0]) );
    sprintf( DBBFileNames[q], DBBFileNamePattern, ZSign, abs(Q[2]), YSign, abs(Q[1]), XSign, abs(Q[0]) );
  }

  multi2d< const char* > Files( NF, NQ );

  for( int q = 0; q < NQ; q ++ )
  {
    f = 0;
    if( HasU == 1 )
    {
      Files[ f ][ q ] = UBBFileNames[ q ];
      f ++;
    }
    if( HasD == 1 )
    {
      Files[ f ][ q ] = DBBFileNames[ q ];
      f ++;
    }
  }

  //#################################################################################//
  // Construct Building Blocks                                                       //
  //#################################################################################//

  Out << "calculating building blocks ..." << "\n";  Out.flush();

  const signed short int T1 = 0;
  const signed short int T2 = NT - 1;

  BuildingBlocks( *B, *F, *U, MaxNLinks, AllLinkPatterns, Phases, Files, T1, T2 );

  Out << "finished calculating building blocks" << "\n";  Out.flush();

  //#################################################################################//
  // Close QDP                                                                       //
  //#################################################################################//

  delete B;
  delete F;
  delete U;

  pop(Xml);
  Xml.close();

  Out << "\n" << "FINISHED" << "\n" << "\n";
  Out.close();

  END_CODE();

  QDP_finalize();

  exit( 0 );
}

//###################################################################################//
//###################################################################################//
