//#####################################################################################
//#####################################################################################
//
// BuildingBlocks.cc
//
//#####################################################################################
//#####################################################################################
//
// description:
//
// Read BuildingBlocks.hh.
//
// history:
//
// There were at least four versions of "MIT" code.  Andrew Pochinsky has a c version.
// Dmitri Dolgov has a c++ version.  Dru B. Renner has c and c++ versions.  All were
// independent and checked against one another.  Of course, all were developed under 
// the guidance of John W. Negele.  This code here is just the "Building Blocks" 
// portion of the MIT code.
//
// authors:
//
// Dru B. Renner, dru@mit.edu, 2002 - port of Building Blocks (MIT) code to qdp++
//
// There are others who have contributed since the code has been migrated to qdp++.
// The cvs log entries indicate these other authors.
//
//#####################################################################################
//#####################################################################################

static const char* const CVSBuildingBlocks_cc =
  "$Header: /home/bjoo/fromJLAB/cvsroot/chroma_base/lib/meas/hadron/BuildingBlocks_w.cc,v 1.7 2004-05-27 16:26:55 dru Exp $";

//#####################################################################################
//#####################################################################################

#include <iostream>
#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/BuildingBlocks_w.h"

using namespace QDP;

#define _DEBUG_BB_C_ 1

//#####################################################################################
// record the CVS info
//#####################################################################################

void CVSBuildingBlocks( TextWriter & Out )
{
  Out << "CVSBuildingBlocks_hh = " << CVSBuildingBlocks_hh << "\n";
  Out << "CVSBuildingBlocks_cc = " << CVSBuildingBlocks_cc << "\n";
}

//#####################################################################################
// backward forward trace
//#####################################################################################

static int GBB_NLinkPatterns = 0;

void BkwdFrwdTr( const LatticePropagator &             B,
                 const LatticePropagator &             F,
                 const SftMom &                        Phases,
                 const multi2d< BinaryWriter* > &      BinaryWriters,
                 const int                             f,
                 const multi1d< unsigned short int > & LinkDirs )
{
  const unsigned short int NLinks = LinkDirs.size();
  unsigned short int Link;
  const int T1 = 0;
  const int T2 = 7;
  const int NQ = Phases.numMom();

  //###################################################################################
  // add a tag to identify the link pattern
  //###################################################################################

  for( int q = 0; q < NQ; q ++ )
  {
    BinaryWriters[f][q]->write( NLinks );

    #if _DEBUG_BB_C_ == 1
    {
      QDPIO::cout << "DEBUG: " << __FILE__ << " " << __LINE__ << "\n";
      QDPIO::cout << "q = " << q << "\n";
      QDPIO::cout << "f = " << f << "\n";
      QDPIO::cout << "NLinks = " << NLinks << "\n";
    }
    #endif

    for( Link = 0; Link < NLinks; Link ++ )
    {
      BinaryWriters[f][q]->write( LinkDirs[ Link ] );

      #if _DEBUG_BB_C_ == 1
      {
        QDPIO::cout << "Link = " << Link << "\n";
        QDPIO::cout << "LinkDirs[ Link ] = " << LinkDirs[ Link ] << "\n";
      }
      #endif
    }

    if( f == 0 )  // counts number of link patterns per flavor - assumes that u quark is always present - potential problem
    {
      GBB_NLinkPatterns ++;
    }
  }

  for( int i = 0; i < Ns * Ns; i ++ )
  {
    LatticeComplex Trace = trace( adj( B ) * Gamma( i ) * F );  // assumes any Gamma5 matrices have already been absorbed

    multi2d< DComplex > Projections = Phases.sft( Trace );

    for( int q = 0; q < NQ; q ++ )
    {
      multi1d< DComplex > Projection = Projections[ q ];

      for( int t = T1; t <= T2; t ++ )
      {
        float r = toFloat( real( Projection[ t ] ) );
        float i = toFloat( imag( Projection[ t ] ) );

        BinaryWriters[f][q]->write( r );
        BinaryWriters[f][q]->write( i );

        #if _DEBUG_BB_C_ == 1
        {
          QDPIO::cout << "DEBUG: " << __FILE__ << " " << __LINE__ << "\n";
          QDPIO::cout << "q = " << q << "\n";
          QDPIO::cout << "f = " << f << "\n";
          QDPIO::cout << "t = " << t << "\n";
          QDPIO::cout << "r = " << r << "\n";
          QDPIO::cout << "i = " << i << "\n";
        }
        #endif
      }
    }
  }

  return;
}

//#####################################################################################
// accumulate link operators
//#####################################################################################

void AddLinks( const multi1d< LatticePropagator > &  B,
               const LatticePropagator &             F,
               const multi1d< LatticeColorMatrix > & U,
               const SftMom &                        Phases,
               multi1d< unsigned short int > &       LinkDirs,
               const unsigned short int              MaxNLinks,
               BBLinkPattern                         LinkPattern,
               const short int                       PreviousDir,
               const short int                       PreviousMu,
               const multi2d< BinaryWriter* > &      BinaryWriters )
{
  const unsigned short int NLinks = LinkDirs.size();

  if( NLinks == MaxNLinks )
  {
    return;
  }

  LatticePropagator F_mu;
  const int NF = B.size();
  multi1d< unsigned short int > NextLinkDirs( NLinks + 1 );
  int Link;

  for( Link = 0; Link < NLinks; Link ++ )
  {
    NextLinkDirs[ Link ] = LinkDirs[ Link ];
  }

  // add link in forward mu direction
  for( int mu = 0; mu < Nd; mu ++ )
  {
    // skip the double back
    if( ( PreviousDir != -1 ) || ( PreviousMu != mu ) )
    {
      bool DoThisPattern = true;
      bool DoFurtherPatterns = true;

      NextLinkDirs[ NLinks ] = mu;

      LinkPattern( DoThisPattern, DoFurtherPatterns, NextLinkDirs );

      if( DoFurtherPatterns == true )
      {
        // accumulate product of link fields
        F_mu = shift( adj( U[ mu ] ) * F, BACKWARD, mu );
      }

      if( DoThisPattern == true )
      {
        // form correlation functions
        for( int f = 0; f < NF; f ++ )
        {
          BkwdFrwdTr( B[ f ], F_mu, Phases, BinaryWriters, f, NextLinkDirs );
        }
      }

      if( DoFurtherPatterns == true )
      {
        // add another link
        AddLinks( B, F_mu, U, Phases, NextLinkDirs, MaxNLinks, LinkPattern, 1, mu, BinaryWriters );
      }
    }
  }

  // add link in backward mu direction
  for( int mu = 0; mu < Nd; mu ++ )
  {
    // skip the double back
    if( ( PreviousDir != 1 ) || ( PreviousMu != mu ) )
    {
      bool DoThisPattern = true;
      bool DoFurtherPatterns = true;

      NextLinkDirs[ NLinks ] = mu + Nd;

      LinkPattern( DoThisPattern, DoFurtherPatterns, NextLinkDirs );

      if( DoFurtherPatterns == true )
      {
        // accumulate product of link fields
        F_mu = U[ mu ] * shift( F, FORWARD, mu );
      }

      if( DoThisPattern == true )
      {
        // form correlation functions
        for( int f = 0; f < NF; f ++ )
        {
          BkwdFrwdTr( B[ f ], F_mu, Phases, BinaryWriters, f, NextLinkDirs );
        }
      }

      if( DoFurtherPatterns == true )
      {
        // add another link
        AddLinks( B, F_mu, U, Phases, NextLinkDirs, MaxNLinks, LinkPattern, -1, mu, BinaryWriters );
      }
    }
  }

  return;
}

//#####################################################################################
// construct building blocks
//#####################################################################################

void BuildingBlocks( const multi1d< LatticePropagator > &  B,
                     const LatticePropagator &             F,
                     const multi1d< LatticeColorMatrix > & U,
                     const unsigned short int              MaxNLinks,
                     const BBLinkPattern                   LinkPattern,
                     const SftMom &                        Phases,
	             const multi2d< const char* > &        BinaryDataFileNames )
{
  //###################################################################################
  // open building blocks data files
  //###################################################################################

  const int NF = B.size();
  const int NQ = Phases.numMom();
  multi2d< BinaryWriter* > BinaryWriters( NF, NQ );

  for( int f = 0; f < NF; f ++ )
  {
  for( int q = 0; q < NQ; q ++ )
  {
    BinaryWriters[ f ][ q ] = new BinaryWriter( BinaryDataFileNames[ f ][ q ] );
  }
  }

  //###################################################################################
  // calculate building blocks
  //###################################################################################

  const unsigned short int NLinks = 0;
  multi1d< unsigned short int > LinkDirs( 0 );

  for( int f = 0; f < NF; f ++ )
  {
    BkwdFrwdTr( B[ f ], F, Phases, BinaryWriters, f, LinkDirs );
  }

  AddLinks( B, F, U, Phases, LinkDirs, MaxNLinks, LinkPattern, 0, -1, BinaryWriters );

  //###################################################################################
  // add footer and close files
  //###################################################################################

  const unsigned short int Id = 0;  // indicates building blocks
  const unsigned short int Version = 1;  // building blocks version
  const unsigned short int Contraction = 0;  // 0 indicates connected diagram
  const unsigned short int NX = Layout::lattSize()[0];
  const unsigned short int NY = Layout::lattSize()[1];
  const unsigned short int NZ = Layout::lattSize()[2];
  const unsigned short int NT = Layout::lattSize()[3];
  const signed short int T1 = 0;
  const signed short int T2 = NT - 1;
  const unsigned short int NLinkPatterns = GBB_NLinkPatterns / NQ;

  for( int f = 0; f < NF; f ++ )
  {
    const signed short int Flavor = f;  // currently assumes u and d are given as f=0 and f=1

  for( int q = 0; q < NQ; q ++ )
  {
    multi1d< int > Q = Phases.numToMom( q );

    const signed short int QX = Q[0];
    const signed short int QY = Q[1];
    const signed short int QZ = Q[2];

    #if _DEBUG_BB_C_ == 1
    {
      QDPIO::cout << "DEBUG: " << __FILE__ << " " << __LINE__ << "\n";

      QDPIO::cout << "Id            = " << Id            << "\n";
      QDPIO::cout << "Version       = " << Version       << "\n";
      QDPIO::cout << "Flavor        = " << Flavor        << "\n";
      QDPIO::cout << "Contraction   = " << Contraction   << "\n";
      QDPIO::cout << "NX            = " << NX            << "\n";
      QDPIO::cout << "NY            = " << NY            << "\n";
      QDPIO::cout << "NZ            = " << NZ            << "\n";
      QDPIO::cout << "NT            = " << NT            << "\n";
      QDPIO::cout << "T1            = " << T1            << "\n";
      QDPIO::cout << "T2            = " << T2            << "\n";
      QDPIO::cout << "MaxNLinks     = " << MaxNLinks     << "\n";
      QDPIO::cout << "NLinkPatterns = " << NLinkPatterns << "\n";
      QDPIO::cout << "QX            = " << QX            << "\n";
      QDPIO::cout << "QY            = " << QY            << "\n";
      QDPIO::cout << "QZ            = " << QZ            << "\n";
    }
    #endif

    // possibly specific to this version
    BinaryWriters[f][q]->write( Flavor );
    BinaryWriters[f][q]->write( Contraction );
    BinaryWriters[f][q]->write( NX );
    BinaryWriters[f][q]->write( NY );
    BinaryWriters[f][q]->write( NZ );
    BinaryWriters[f][q]->write( NT );
    BinaryWriters[f][q]->write( T1 );
    BinaryWriters[f][q]->write( T2 );
    BinaryWriters[f][q]->write( MaxNLinks );
    BinaryWriters[f][q]->write( NLinkPatterns );
    BinaryWriters[f][q]->write( QX );
    BinaryWriters[f][q]->write( QY );
    BinaryWriters[f][q]->write( QZ );
    // generic to any building blocks file
    BinaryWriters[f][q]->write( Id );
    BinaryWriters[f][q]->write( Version );
    // close file
    BinaryWriters[f][q]->close();
    delete BinaryWriters[f][q];
  }
  }

  return;
}

//#####################################################################################
//#####################################################################################

#undef _DEBUG_BB_C_

//#####################################################################################
//#####################################################################################
