//#####################################################################################
//#####################################################################################
//
// ExampleBuildingBlocks.cc
//
//#####################################################################################
//#####################################################################################
//
// description:
//
// Executing ExampleBuildingBlocks without any arguments (or an incorrect number of
// arguments) will result in the required argument list being printed to the stdout.
//
// Note that the file name patterns for the u and d building blocks must be a string
// of the form "...%c%i...%c%i...%c%i..." where the first %c%i corresponds to qz, the
// second to qy, and the third to qx.
//
// authors:
//
// Dru B. Renner, dru@mit.edu, 2002 - port of "MIT" code to qdp++
// Dru B. Renner, dru@mit.edu, July 2003 - analyze MILC configurations
//
//#####################################################################################
//#####################################################################################

static const char* const CVSExampleBuildingBlocks_hh =
  "$Header: /home/bjoo/fromJLAB/cvsroot/chroma_base/mainprogs/main/ExampleBuildingBlocks.cc,v 1.5 2003-10-01 02:59:08 edwards Exp $";

//#####################################################################################
//#####################################################################################

#include <iostream>
#include <string.h>
#include <assert.h>
#include "qdp.h"
#include "qdp_util.h"
#include "chroma.h"

using namespace QDP;

//#####################################################################################
// Accept All Link Patterns
//#####################################################################################

void AllLinkPatterns( bool &                          DoThisPattern,
                      bool &                          DoFurtherPatterns,
                      multi1d< unsigned short int > & LinkPattern )
{
  DoThisPattern     = true;
  DoFurtherPatterns = true;

  return;
}

//#####################################################################################
//#####################################################################################

int main( int argc, char** argv )
{
  //#####################################################################################
  // Check Arguments
  //#####################################################################################

  if( argc != 17 )
  {
    // cout will create multiple outputs for multi-node code
    cout << ""                                                                                               << endl;
    cout << "Arguments:"                                                                                     << endl;
    cout << ""                                                                                               << endl;
    cout << "(1)  NX                               - \"1\", \"2\", ..."                                      << endl;
    cout << "(2)  NY                               - \"1\", \"2\", ..."                                      << endl;
    cout << "(3)  NZ                               - \"1\", \"2\", ..."                                      << endl;
    cout << "(4)  NT                               - \"1\", \"2\", ..."                                      << endl;
    cout << "(5)  Gauge Field                      - Input File Name"                                        << endl;
    cout << "(6)  Gauge Field Format               - \"archive\" or \"szin\""                                << endl;
    cout << "(7)  Forward Propagator               - Input File Name"                                        << endl;
    cout << "(8)  Proton Sequential U Propagator   - Input File Name (\"NULL\" if u not available)"          << endl;
    cout << "(9)  Sequential U Propagator Format   - \"B\", \"G5_B\", \"B_G5\", or \"G5_B_G5\" "             << endl;
    cout << "(10) Proton Sequential D Propagator   - Input File Name (\"NULL\" if d not available)"          << endl;
    cout << "(11) Sequential D Propagator Format   - \"B\", \"G5_B\", \"B_G5\", or \"G5_B_G5\" "             << endl;
    cout << "(12) Maximum Number of Links          - \"0\", \"1\", ..."                                      << endl;
    cout << "(13) Maximum Spatial Momentum Squared - \"0\", \"1\", ..."                                      << endl;
    cout << "(14) Proton U Building Blocks         - Output File Name Pattern (\"NULL\" if u not available)" << endl;
    cout << "(15) Proton D Building Blocks         - Output File Name Pattern (\"NULL\" if d not available)" << endl;
    cout << "(16) Text Output File Name            - Output File Name"                                       << endl;
    cout << ""                                                                                               << endl;

    exit( 1 );
  }

  //#####################################################################################
  // Capture Arguments
  //#####################################################################################

  const unsigned short int NX            = atoi( argv[1] );
  const unsigned short int NY            = atoi( argv[2] );
  const unsigned short int NZ            = atoi( argv[3] );
  const unsigned short int NT            = atoi( argv[4] );
  const char* const GaugeFieldFileName   = argv[5];
  const char* const GaugeFieldFormat     = argv[6];
  const char* const FrwdPropFileName     = argv[7];
  const char* const SeqUPropFileName     = argv[8];
  const char* const SeqUPropFormat       = argv[9];
  const char* const SeqDPropFileName     = argv[10];
  const char* const SeqDPropFormat       = argv[11];
  const unsigned short int MaxNLinks     = atoi( argv[12] );
  const unsigned short int MaxQSquared   = atoi( argv[13] );
  const char* const UBBFileNamePattern   = argv[14];
  const char* const DBBFileNamePattern   = argv[15];
  const char* const OutFileName          = argv[16];

  //#####################################################################################
  // Initialize qdp
  //#####################################################################################

  QDP_initialize( & argc, & argv );

  //#####################################################################################
  // Echo Arguments
  //#####################################################################################

  // will capture all would-be standard output
  TextWriter Out( OutFileName );

  Out <<                                                                "\n";
  Out << "NX                               = " << NX                 << "\n";
  Out << "NY                               = " << NY                 << "\n";
  Out << "NZ                               = " << NZ                 << "\n";
  Out << "NT                               = " << NT                 << "\n";
  Out << "Gauge Field                      = " << GaugeFieldFileName << "\n";
  Out << "Gauge Field Format               = " << GaugeFieldFormat   << "\n";
  Out << "Forward Propagator               = " << FrwdPropFileName   << "\n";
  Out << "Proton Sequential U Propagator   = " << SeqUPropFileName   << "\n";
  Out << "Sequential U Propagator Format   = " << SeqUPropFormat     << "\n";
  Out << "Proton Sequential D Propagator   = " << SeqDPropFileName   << "\n";
  Out << "Sequential D Propagator Format   = " << SeqDPropFormat     << "\n";
  Out << "Maximum Number of Links          = " << MaxNLinks          << "\n";
  Out << "Maximum Spatial Momentum Squared = " << MaxQSquared        << "\n";
  Out << "Proton U Building Blocks         = " << UBBFileNamePattern << "\n";
  Out << "Proton D Building Blocks         = " << DBBFileNamePattern << "\n";
  Out << "Text Output File Name            = " << OutFileName        << "\n";
  Out <<                                                                "\n";
  Out.flush();

  //#####################################################################################
  // Record the CVS Info for BuildingBlocks
  //#####################################################################################

  CVSBuildingBlocks( Out );
  Out << "CVSExampleBuildingBlocks_hh = " << CVSExampleBuildingBlocks_hh << "\n";
  Out <<                                                                "\n";

  //#####################################################################################
  // Set Lattice Geometry
  //#####################################################################################

  multi1d< int > N ( Nd );

  N[0] = NX;
  N[1] = NY;
  N[2] = NZ;
  N[3] = NT;

  Layout::setLattSize( N );
  Layout::create();

  //#####################################################################################
  // Read Propagators
  //#####################################################################################

  // Forward Quark Propagator
  XMLReader ForwardPropXML;
  LatticePropagator* F = NULL;
  F = new LatticePropagator;
  Out << "reading forward propagator " << FrwdPropFileName << " ... ";  Out.flush();
  readSzinQprop( ForwardPropXML, *F, FrwdPropFileName );
  Out << "done." << "\n";  Out.flush();

  int NF = 0;

  // Corresponding Sequential or Backward u Type Quark Propagator
  int HasU = 0;
  LatticePropagator* Bu = NULL;
  if( strcmp( SeqUPropFileName, "NULL" ) != 0 )
  {
    Bu = new LatticePropagator;
    Out << "reading proton sequential u propagator " << SeqUPropFileName << " ... ";  Out.flush();
    readSzinQprop( ForwardPropXML, *Bu, SeqUPropFileName );
    Out << "done." << "\n";  Out.flush();
    HasU = 1;
    NF ++;
  }

  // Corresponding Sequential or Backward d Type Quark Propagator
  XMLReader BackwardPropXML;
  int HasD = 0;
  LatticePropagator* Bd = NULL;
  if( strcmp( SeqDPropFileName, "NULL" ) != 0 )
  {
    Bd = new LatticePropagator;
    Out << "reading proton sequential d propagator " << SeqDPropFileName << " ... ";  Out. flush();
    readSzinQprop( BackwardPropXML, *Bd, SeqDPropFileName );
    Out << "done." << "\n";  Out.flush();
    HasD = 1;
    NF ++;
  }

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

  //#####################################################################################
  // Convert Sequential Propagator Format
  //#####################################################################################

  if( HasU == 1 )
  {
  if( strcmp( SeqUPropFormat, "G5_B" ) == 0 )
  {
    *Bu = Gamma( 15 ) * *Bu;
  }
  if( strcmp( SeqUPropFormat, "B_G5" ) == 0 )
  {
    *Bu = *Bu * Gamma( 15 );
  }
  if( strcmp( SeqUPropFormat, "G5_B_G5" ) == 0 )
  {
    *Bu = Gamma( 15 ) * *Bu * Gamma( 15 );
  }
  }

  if( HasD == 1 )
  {
  if( strcmp( SeqDPropFormat, "G5_B" ) == 0 )
  {
    *Bd = Gamma( 15 ) * *Bd;
  }
  if( strcmp( SeqDPropFormat, "B_G5" ) == 0 )
  {
    *Bd = *Bd * Gamma( 15 );
  }
  if( strcmp( SeqDPropFormat, "G5_B_G5" ) == 0 )
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

  //#####################################################################################
  // Read Gauge Field
  //#####################################################################################

  XMLReader GaugeXML;
  multi1d< LatticeColorMatrix >* U = NULL;
  U = new multi1d< LatticeColorMatrix >( Nd );
  Out << "reading gauge field " << GaugeFieldFileName << " ... ";  Out.flush();
  if( strcmp( GaugeFieldFormat, "archive" ) == 0 )
  {
    readArchiv( GaugeXML, *U, GaugeFieldFileName );
  }
  if( strcmp( GaugeFieldFormat, "szin" ) == 0 )
  {
    readSzin( GaugeXML, *U, GaugeFieldFileName );
  }
  Out << "done." << "\n";  Out.flush();

  //#####################################################################################
  // Set Momenta
  //#####################################################################################

  multi1d< int > SnkMom( ND - 1 );
  SnkMom = 0;

  SftMom Phases( MaxQSquared, SnkMom, false, ND - 1 );

  //#####################################################################################
  // Construct File Names
  //#####################################################################################

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

  //#####################################################################################
  // Construct Building Blocks
  //#####################################################################################

  Out << "calculating building blocks ... ";  Out.flush();

  BuildingBlocks( *B, *F, *U, MaxNLinks, AllLinkPatterns, Phases, Files );

  Out << "done." << "\n";  Out.flush();

  //#####################################################################################
  // Close qdp
  //#####################################################################################

  delete B;
  delete F;
  delete U;

  Out << "\n" << "FINISHED" << "\n" << "\n";
  Out.close();

  QDP_finalize();

  exit(0);
}

//#####################################################################################
//#####################################################################################
