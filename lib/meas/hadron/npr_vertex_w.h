// -*- C++ -*-
//  $Id: npr_vertex_w.h,v 1.3 2006-11-02 22:26:10 edwards Exp $
/*! \file
 *  \brief NPR vertex calculations
 */

#ifndef __npr_vertex_w_h__
#define __npr_vertex_w_h__

namespace Chroma 
{
  //! Used to Set Requested Link Patterns
  /*! \ingroup hadron */
  typedef void (*BBLinkPattern)(bool &                          DoThisPattern,
				bool &                          DoFurtherPatterns,
				multi1d< unsigned short int > & LinkPattern);

  //! NPR vertices
  /*! \ingroup hadron */
  void NprVertex(const LatticePropagator &             F,
		 const multi1d< LatticeColorMatrix > & U,
		 const unsigned short int              MaxNLinks,
		 const BBLinkPattern                   LinkPattern,
		 QDPFileWriter& qio_file);

}  // end namespace Chroma

#endif

//###################################################################################//
//###################################################################################//
