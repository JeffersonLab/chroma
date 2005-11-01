//###################################################################################//
//###################################################################################//
//                                                                                   //
// BuildingBlocks.cc                                                                 //
//                                                                                   //
//###################################################################################//
//###################################################################################//
//                                                                                   //
// description:                                                                      //
//                                                                                   //
// Read BuildingBlocks.hh.                                                           //
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
#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/BuildingBlocks_w.h"

namespace Chroma {


//###################################################################################//
// debug flag                                                                        //
//###################################################################################//

#define _DEBUG_BB_C_ 0

//###################################################################################//
// cvs header                                                                        //
//###################################################################################//

static const char* const CVSBuildingBlocks_cc =
  "$Header: /home/bjoo/fromJLAB/cvsroot/chroma_base/lib/meas/hadron/BuildingBlocks_w.cc,v 2.4 2005-11-01 04:17:29 edwards Exp $";

//###################################################################################//
// record the CVS info                                                               //
//###################################################################################//

void CVSBuildingBlocks( TextWriter & Out )
{
  Out << "CVSBuildingBlocks_hh = " << CVSBuildingBlocks_hh << "\n";
  Out << "CVSBuildingBlocks_cc = " << CVSBuildingBlocks_cc << "\n";
}

//###################################################################################//
// backward forward trace                                                            //
//###################################################################################//

void BkwdFrwdTr( const LatticePropagator &             B,
                 const LatticePropagator &             F,
		 int                                   GammaInsertion,
                 const SftMom &                        Phases,
                 const SftMom &                        PhasesCanonical,
                 multi2d< BinaryWriter > &             BinaryWriters,
		 multi1d< int > &                      GBB_NLinkPatterns,
		 multi2d< int > &                      GBB_NMomPerms,
                 const int                             f,
                 const multi1d< unsigned short int > & LinkDirs,
                 const signed short int T1, 
                 const signed short int T2 )
{
  StopWatch TotalTime;
  TotalTime.reset();
  TotalTime.start();

  StopWatch Timer;

  int TRCalls = 0;
  int FTCalls = 0;
  int GFGCalls = 0;
  int IPCalls = 0;
  double TRTime = 0.0;
  double FTTime = 0.0;
  double GFGTime = 0.0;
  double IPTime = 0.0;
  double IOTime = 0.0;

  const unsigned short int NLinks = LinkDirs.size();
  unsigned short int Link;
  const int NumQ = Phases.numMom();
  const int NumO = BinaryWriters.size1();

  //#################################################################################//
  // add a tag to identify the link pattern                                          //
  //#################################################################################//

  Timer.reset();
  Timer.start();

  for( int o = 0; o < NumO; o ++ )
  {
    BinaryWriters(f,o).write( NLinks );

    #if _DEBUG_BB_C_ == 1
    {
      QDPIO::cout << "DEBUG: " << __FILE__ << " " << __LINE__ << "\n";
      QDPIO::cout << "q = " << o << "\n";
      QDPIO::cout << "f = " << f << "\n";
      QDPIO::cout << "NLinks = " << NLinks << "\n";
    }
    #endif

    for( Link = 0; Link < NLinks; Link ++ )
    {
      BinaryWriters(f,o).write( LinkDirs[ Link ] );

      #if _DEBUG_BB_C_ == 1
      {
        QDPIO::cout << "Link = " << Link << "\n";
        QDPIO::cout << "LinkDirs[ Link ] = " << LinkDirs[ Link ] << "\n";
      }
      #endif
    }

     // counts number of link patterns per flavor
    GBB_NLinkPatterns[f] ++;
  }

  Timer.stop();
  IOTime += Timer.getTimeInSeconds();

  for( int i = 0; i < Ns * Ns; i ++ )
  {
    Timer.reset();
    Timer.start();

    // assumes any Gamma5 matrices have already been absorbed
    //LatticeComplex Trace = trace( adj( B ) * Gamma( i ) * F * Gamma( GammaInsertion ) );
    LatticePropagator GFG = Gamma(i) * F * Gamma( GammaInsertion );

    Timer.stop();
    GFGCalls += 1;
    GFGTime += Timer.getTimeInSeconds();
    Timer.reset();
    Timer.start();

    LatticeComplex Trace = localInnerProduct( B, GFG );

    Timer.stop();
    IPCalls += 1;
    IPTime += Timer.getTimeInSeconds();
    Timer.reset();
    Timer.start();

    multi2d< DComplex > Projections = Phases.sft( Trace );

    Timer.stop();
    FTTime += Timer.getTimeInSeconds();
    FTCalls += 1;

    Timer.reset();
    Timer.start();

    for( int q = 0; q < NumQ; q ++ )
    {
      multi1d< DComplex > Projection = Projections[ q ];
      multi1d< int > Q = Phases.numToMom( q );

      int o = PhasesCanonical.momToNum( Q );
      if (o == -1)
      {
	QDPIO::cerr << __func__ << ": internal error: failed to find index of ordered momentum" << endl;
	QDP_abort(1);
      }
      
      const signed short int QX = Q[0];
      const signed short int QY = Q[1];
      const signed short int QZ = Q[2];
      BinaryWriters(f,o).write( QX );
      BinaryWriters(f,o).write( QY );
      BinaryWriters(f,o).write( QZ );

      // counts number of momenta permutations per canonical ordering
      GBB_NMomPerms(f,o) ++;

      // Write correlator
      for( int t = T1; t <= T2; t ++ )
      {
        float r = toFloat( real( Projection[ t ] ) );
        float i = toFloat( imag( Projection[ t ] ) );

        BinaryWriters(f,o).write( r );
        BinaryWriters(f,o).write( i );

        #if _DEBUG_BB_C_ == 1
        {
          QDPIO::cout << "DEBUG: " << __FILE__ << " " << __LINE__ << "\n";
          QDPIO::cout << "q = " << q << "\n";
          QDPIO::cout << "o = " << o << "\n";
          QDPIO::cout << "f = " << f << "\n";
          QDPIO::cout << "t = " << t << "\n";
          QDPIO::cout << "r = " << r << "\n";
          QDPIO::cout << "i = " << i << "\n";
        }
        #endif
      }
    }

    Timer.stop();
    IOTime += Timer.getTimeInSeconds();
  }

  QDPIO::cout << __func__ << ":  io time = " << IOTime << " seconds" << endl;
  QDPIO::cout << __func__ << ": gfg time = " << GFGTime / (double) GFGCalls << " seconds" << endl;
  QDPIO::cout << __func__ << ":  ip time = " << IPTime / (double) IPCalls << " seconds" << endl;
  QDPIO::cout << __func__ << ":  ft time = " << FTTime / (double) FTCalls << " seconds" << endl;
  TotalTime.stop();
  QDPIO::cout << __func__ << ": total time = " << TotalTime.getTimeInSeconds() << " seconds" << endl;

  return;
}

//###################################################################################//
// accumulate link operators                                                         //
//###################################################################################//

void AddLinks( const multi1d< LatticePropagator > &  B,
               const LatticePropagator &             F,
               const multi1d< LatticeColorMatrix > & U,
	       const multi1d< int > &                GammaInsertions,
               const SftMom &                        Phases,
	       const SftMom &                        PhasesCanonical,
               multi1d< unsigned short int > &       LinkDirs,
               const unsigned short int              MaxNLinks,
               BBLinkPattern                         LinkPattern,
               const short int                       PreviousDir,
               const short int                       PreviousMu,
               multi2d< BinaryWriter > &             BinaryWriters,
	       multi1d< int > &                      GBB_NLinkPatterns,
	       multi2d< int > &                      GBB_NMomPerms,
               const signed short int T1, 
               const signed short int T2 )
{
  StopWatch Timer;
  int ShiftCalls = 0;
  double ShiftTime = 0.0;

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
        Timer.reset();
        Timer.start();

        // accumulate product of link fields
        F_mu = shift( adj( U[ mu ] ) * F, BACKWARD, mu );

        Timer.stop();
        ShiftTime += Timer.getTimeInSeconds();
        ShiftCalls += 1;
      }

      if( DoThisPattern == true )
      {
        // form correlation functions
        for( int f = 0; f < NF; f ++ )
        {
          BkwdFrwdTr( B[ f ], F_mu, GammaInsertions[ f ], Phases, PhasesCanonical,
		      BinaryWriters, GBB_NLinkPatterns, GBB_NMomPerms,
		      f, NextLinkDirs, T1, T2 );
        }
      }

      if( DoFurtherPatterns == true )
      {
        // add another link
        AddLinks( B, F_mu, U, GammaInsertions, 
		  Phases, PhasesCanonical,
		  NextLinkDirs, MaxNLinks, LinkPattern, 1, mu, 
		  BinaryWriters, GBB_NLinkPatterns, GBB_NMomPerms,
		  T1, T2 );
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
        Timer.reset();
        Timer.start();

        // accumulate product of link fields
        F_mu = U[ mu ] * shift( F, FORWARD, mu );

        Timer.stop();
        ShiftTime += Timer.getTimeInSeconds();
        ShiftCalls += 1;
      }

      if( DoThisPattern == true )
      {
        // form correlation functions
        for( int f = 0; f < NF; f ++ )
        {
          BkwdFrwdTr( B[ f ], F_mu, GammaInsertions[ f ], Phases, PhasesCanonical,
		      BinaryWriters, GBB_NLinkPatterns, GBB_NMomPerms,
		      f, NextLinkDirs, T1, T2 );
        }
      }

      if( DoFurtherPatterns == true )
      {
        // add another link
        AddLinks( B, F_mu, U, GammaInsertions, Phases, PhasesCanonical,
		  NextLinkDirs, MaxNLinks, LinkPattern, -1, mu, BinaryWriters, 
		  GBB_NLinkPatterns, GBB_NMomPerms,
		  T1, T2 );
      }
    }
  }

  QDPIO::cout << __func__ << ": shift time = "
	      << ShiftTime
	      << " seconds with shift calls = "
              << ShiftCalls 
	      << endl;

  return;
}

//###################################################################################//
// construct building blocks                                                         //
//###################################################################################//

void BuildingBlocks( const multi1d< LatticePropagator > &  B,
                     const LatticePropagator &             F,
                     const multi1d< LatticeColorMatrix > & U,
                     const multi1d< int > &                GammaInsertions,
                     const multi1d< int > &                Flavors,
                     const unsigned short int              MaxNLinks,
                     const BBLinkPattern                   LinkPattern,
                     const SftMom &                        Phases,
                     const SftMom &                        PhasesCanonical,
	             const multi2d< string > &             BinaryDataFileNames,
                     const signed short int                T1,
                     const signed short int                T2,
		     const std::string&                    SeqSourceType, 
		     const multi1d< int >&                 SnkMom, 
		     const signed short int                DecayDir)
{
  StopWatch TotalTime;
  TotalTime.reset();
  TotalTime.start();

  StopWatch Timer;

  //#################################################################################//
  // open building blocks data files                                                 //
  //#################################################################################//

  Timer.reset();
  Timer.start();

  if ( BinaryDataFileNames.size1() != PhasesCanonical.numMom() )
  {
    QDPIO::cerr << __func__ << ": number of datafile names incompatible with number of momenta"
		<< endl;
    QDP_abort(1);
  }

  const int NumF = B.size();
  const int NumO = BinaryDataFileNames.size1();
  multi2d< BinaryWriter > BinaryWriters( NumF, NumO );
  multi1d< int > GBB_NLinkPatterns( NumF );
  multi2d< int > GBB_NMomPerms( NumF, NumO );

  for( int f = 0; f < NumF; f ++ )
  {
    GBB_NLinkPatterns[f] = 0;
    for( int o = 0; o < NumO; o ++ )
    {
      BinaryWriters(f,o).open( BinaryDataFileNames(f,o) );
      GBB_NMomPerms(f,o) = 0;
    }
  }

  Timer.stop();
  QDPIO::cout << __func__ << ": time to open files = "
	      << Timer.getTimeInSeconds() 
	      << " seconds" << endl;

  //#################################################################################//
  // calculate building blocks                                                       //
  //#################################################################################//

  Timer.reset();
  Timer.start();

  QDPIO::cout << __func__ << ": start BkwdFrwdTr" << endl;

  const unsigned short int NLinks = 0;
  multi1d< unsigned short int > LinkDirs( 0 );

  for( int f = 0; f < NumF; f ++ )
  {
    BkwdFrwdTr( B[ f ], F, GammaInsertions[ f ], Phases, PhasesCanonical,
		BinaryWriters, GBB_NLinkPatterns, GBB_NMomPerms, f, LinkDirs, T1, T2 );
  }

  Timer.stop();
  QDPIO::cout << __func__ << ": total time for 0 links (single BkwdFrwdTr call) = "
	      << Timer.getTimeInSeconds() 
	      << " seconds" << endl;

  Timer.reset();
  Timer.start();

  QDPIO::cout << __func__ << ": start AddLinks" << endl;

  AddLinks( B, F, U, GammaInsertions, 
	    Phases, PhasesCanonical,
	    LinkDirs, MaxNLinks, LinkPattern, 0, -1, 
	    BinaryWriters, GBB_NLinkPatterns, GBB_NMomPerms,
	    T1, T2 );

  Timer.stop();
  QDPIO::cout << __func__ << ": total time for remaining links (outermost AddLinks call) = "
	      << Timer.getTimeInSeconds() 
	      << " seconds" << endl;

  //#################################################################################//
  // add footer and close files                                                      //
  //#################################################################################//

  Timer.reset();
  Timer.start();

  const unsigned short int Id = 0;  // indicates building blocks
  const unsigned short int Version = 3;  // building blocks version
  const unsigned short int Contraction = 0;  // 0 indicates connected diagram
  const unsigned short int NX = Layout::lattSize()[0];
  const unsigned short int NY = Layout::lattSize()[1];
  const unsigned short int NZ = Layout::lattSize()[2];
  const unsigned short int NT = Layout::lattSize()[3];
  //const signed short int T1 = 0;
  //const signed short int T2 = NT - 1;
  const signed short int   PX = SnkMom[0];
  const signed short int   PY = SnkMom[1];
  const signed short int   PZ = SnkMom[2];
  const signed short int   SeqSourceLen = 64;
  string SeqSource = SeqSourceType;
  SeqSource.resize(SeqSourceLen, 0);

  for( int f = 0; f < NumF; f ++ )
  {
    const signed short int Flavor = Flavors[f];  // currently assumes u and d are given as f=0 and f=1
    const signed short int GammaInsertion = GammaInsertions[f];
    const unsigned short int NLinkPatterns = GBB_NLinkPatterns[f] / NumO;

    for( int o = 0; o < NumO; o ++ )
    {
      const unsigned short int NMomPerms = GBB_NMomPerms(f,o) / (Ns * Ns * NLinkPatterns);

#if _DEBUG_BB_C_ == 1
      {
	multi1d< int > Q = PhasesCanonical.numToMom( o );
	const signed short int QX = Q[0];
	const signed short int QY = Q[1];
	const signed short int QZ = Q[2];

	QDPIO::cout << "DEBUG: " << __FILE__ << " " << __LINE__ << "\n";

	QDPIO::cout << "Id              = " << Id            << "\n";
	QDPIO::cout << "Version         = " << Version       << "\n";
	QDPIO::cout << "Flavor          = " << Flavor        << "\n";
	QDPIO::cout << "Contraction     = " << Contraction   << "\n";
	QDPIO::cout << "SeqSource       = " << SeqSource     << "\n";
	QDPIO::cout << "GammaInsertion  = " << GammaInsertion<< "\n";
	QDPIO::cout << "NX              = " << NX            << "\n";
	QDPIO::cout << "NY              = " << NY            << "\n";
	QDPIO::cout << "NZ              = " << NZ            << "\n";
	QDPIO::cout << "NT              = " << NT            << "\n";
	QDPIO::cout << "T1              = " << T1            << "\n";
	QDPIO::cout << "T2              = " << T2            << "\n";
	QDPIO::cout << "MaxNLinks       = " << MaxNLinks     << "\n";
	QDPIO::cout << "NLinkPatterns   = " << NLinkPatterns << "\n";
	QDPIO::cout << "NMomPerms       = " << NMomPerms     << "\n";
	QDPIO::cout << "Canonical Index = " << o             << "\n";
	QDPIO::cout << "QX              = " << QX            << "\n";
	QDPIO::cout << "QY              = " << QY            << "\n";
	QDPIO::cout << "QZ              = " << QZ            << "\n";
	QDPIO::cout << "PX              = " << PX            << "\n";
	QDPIO::cout << "PY              = " << PY            << "\n";
	QDPIO::cout << "PZ              = " << PZ            << "\n";
      }
#endif

      // possibly specific to this version
      BinaryWriters(f,o).write( Flavor );
      BinaryWriters(f,o).write( Contraction );
      BinaryWriters(f,o).writeArray( SeqSource.data(), 1, SeqSourceLen );
      BinaryWriters(f,o).write( GammaInsertion );
      BinaryWriters(f,o).write( NX );
      BinaryWriters(f,o).write( NY );
      BinaryWriters(f,o).write( NZ );
      BinaryWriters(f,o).write( NT );
      BinaryWriters(f,o).write( DecayDir );
      BinaryWriters(f,o).write( T1 );
      BinaryWriters(f,o).write( T2 );
      BinaryWriters(f,o).write( MaxNLinks );
      BinaryWriters(f,o).write( NLinkPatterns );
      BinaryWriters(f,o).write( NMomPerms );
      BinaryWriters(f,o).write( PX );
      BinaryWriters(f,o).write( PY );
      BinaryWriters(f,o).write( PZ );
      BinaryWriters(f,o).write( BinaryWriters(f,o).getChecksum() );
      // generic to any building blocks file
      BinaryWriters(f,o).write( Id );
      BinaryWriters(f,o).write( Version );
      // close file
      BinaryWriters(f,o).close();
    }
  }

  Timer.stop();
  QDPIO::cout << __func__ << ": time to write footer = "
	      << Timer.getTimeInSeconds() 
	      << " seconds" << endl;

  TotalTime.stop();
  QDPIO::cout << __func__ << ": total time = "
	      << TotalTime.getTimeInSeconds() 
	      << " seconds" << endl;

  return;
}

//###################################################################################//
//###################################################################################//

#undef _DEBUG_BB_C_

//###################################################################################//
//###################################################################################//

}  // end namespace Chroma
