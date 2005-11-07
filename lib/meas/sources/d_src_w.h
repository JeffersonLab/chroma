// -*- C++ -*-
// $Id: d_src_w.h,v 2.1 2005-11-07 06:27:47 edwards Exp $
/*! \file
 *  \brief D-wave sources
 */

#ifndef __d_src_h__
#define __d_src_h__

namespace Chroma 
{

  //! D-wave source
  /*! @ingroup sources */
  void d_src(const multi1d<LatticeColorMatrix>& u, 
	     LatticeColorVector& chi,
	     int direction);

  //! D-wave source
  /*! @ingroup sources */
  void d_src(const multi1d<LatticeColorMatrix>& u, 
	     LatticePropagator& chi,
	     int direction);

  //! D-wave source
  /*! @ingroup sources */
  void d_src(const multi1d<LatticeColorMatrix>& u, 
	     LatticeFermion& chi,
	     int direction);

}  // end namespace Chroma
#endif
