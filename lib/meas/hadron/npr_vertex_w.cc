//  $Id: npr_vertex_w.cc,v 1.3 2006-11-03 20:18:55 hwlin Exp $
/*! \file
 *  \brief NPR vertex calculations
 */

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/npr_vertex_w.h"

using namespace QDP;

namespace Chroma 
{

  void BkwdFrwd(const LatticePropagator&  B,
		const LatticePropagator&  F,
		QDPFileWriter& qio_file,
		int& GBB_NLinkPatterns,
		const multi1d< int > & LinkDirs)
  {
    StopWatch TotalTime;
    TotalTime.reset();
    TotalTime.start();

    for( int i = 0; i < Ns * Ns; i ++ )
    {
      XMLBufferWriter record_xml;
      push(record_xml, "Vertex");

      write(record_xml, "linkDirs", LinkDirs);   // link pattern
      write(record_xml, "gamma", i);

      // counts number of link patterns
      GBB_NLinkPatterns++;

      // assumes any Gamma5 matrices have already been absorbed
      int G5 = Ns*Ns-1;
      
      // This is a beastly hack to make the thing compile as QIO isn't wanting
      // to write just one single Site's worth of DPropagator
      // So I take a lattice propagetor and set every site's data equal to the 
      // site that I am interested in. THis is just to make the code compile
      // and should be fixed at a later time
      DPropagator prop;
      {
	LatticePropagator tmp = (Gamma(G5) * adj(B) * Gamma(G5)) * Gamma(i) * F;
	prop = sum(tmp);   // The site's worth of data of interest
      }
      
      // !!!! FIXME
      // Hack set every site's worth of data equal to 'prop' 
      LatticePropagatorD dprop = zero;
      dprop = prop;
      write(qio_file, record_xml, dprop);
    }

    TotalTime.stop();
    QDPIO::cout << __func__ << ": total time = " << TotalTime.getTimeInSeconds() << " seconds" << endl;

    return;
  }

//###################################################################################//
// accumulate link operators                                                         //
//###################################################################################//

  void AddLinks(const LatticePropagator&  B,
		const LatticePropagator&  F,
		const multi1d< LatticeColorMatrix > & U,
		multi1d< int >&    LinkDirs,
		const int          MaxNLinks,
		BBLinkPattern      LinkPattern,
		const int          PreviousDir,
		const int          PreviousMu,
		QDPFileWriter&     qio_file,
		int&               GBB_NLinkPatterns)
  {
    StopWatch Timer;
    Timer.reset();
    Timer.start();

    const int NLinks = LinkDirs.size();

    if( NLinks == MaxNLinks )
    {
      return;
    }

    LatticePropagator F_mu;
    multi1d< int > NextLinkDirs( NLinks + 1 );

    for(int Link = 0; Link < NLinks; Link ++)
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
	  BkwdFrwd(B, F_mu, qio_file, GBB_NLinkPatterns, NextLinkDirs);
	}

	if( DoFurtherPatterns == true )
	{
	  // add another link
	  AddLinks(B, F_mu, U,
		   NextLinkDirs, MaxNLinks, LinkPattern, 1, mu, 
		   qio_file, GBB_NLinkPatterns);
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
	  BkwdFrwd(B, F_mu, qio_file, GBB_NLinkPatterns, NextLinkDirs);
	}

	if( DoFurtherPatterns == true )
	{
	  // add another link
	  AddLinks(B, F_mu, U, 
		   NextLinkDirs, MaxNLinks, LinkPattern, -1, mu, 
		   qio_file, GBB_NLinkPatterns);
	}
      }
    }

    Timer.stop();
    QDPIO::cout << __func__ << ": total time = " << Timer.getTimeInSeconds() << " seconds" << endl;

    return;
  }


  //! NPR vertices
  void NprVertex(const LatticePropagator &             F,
		 const multi1d< LatticeColorMatrix > & U,
		 const unsigned short int              MaxNLinks,
		 const BBLinkPattern                   LinkPattern,
		 QDPFileWriter& qio_file)
  {
    StopWatch TotalTime;
    TotalTime.reset();
    TotalTime.start();

    StopWatch Timer;

    int GBB_NLinkPatterns;

    //#################################################################################//
    // open building blocks data files                                                 //
    //#################################################################################//

    Timer.reset();
    Timer.start();

    //#################################################################################//
    // calculate building blocks                                                       //
    //#################################################################################//

    QDPIO::cout << __func__ << ": start BkwdFrwd" << endl;

    const int NLinks = 0;
    multi1d< int > LinkDirs( 0 );
    {
      LatticePropagator B = Gamma(15)*adj(F)*Gamma(15);

      BkwdFrwd(B, F, qio_file, GBB_NLinkPatterns, LinkDirs);
    }

    Timer.stop();
    QDPIO::cout << __func__ << ": total time for 0 links (single BkwdFrwdTr call) = "
		<< Timer.getTimeInSeconds() 
		<< " seconds" << endl;

    Timer.reset();
    Timer.start();

    QDPIO::cout << __func__ << ": start AddLinks" << endl;

    AddLinks(F, F, U, 
	     LinkDirs, MaxNLinks, LinkPattern, 0, -1, 
	     qio_file, GBB_NLinkPatterns);

    Timer.stop();
    QDPIO::cout << __func__ << ": total time for remaining links (outermost AddLinks call) = "
		<< Timer.getTimeInSeconds() 
		<< " seconds" << endl;

    TotalTime.stop();
    QDPIO::cout << __func__ << ": total time = "
		<< TotalTime.getTimeInSeconds() 
		<< " seconds" << endl;

    return;
  }

}  // end namespace Chroma
