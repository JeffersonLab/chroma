// -*- C++ -*-
// $Id: p_src_w.h,v 2.1 2005-11-07 06:27:47 edwards Exp $
/*! \file
 *  \brief P-wave sources
 */

#ifndef __p_src_h__
#define __p_src_h__

namespace Chroma 
{
  //! P-wave source
  /*! @ingroup sources */
  void p_src(const multi1d<LatticeColorMatrix>& u, 
	     LatticeColorVector& chi,
	     int direction);

  //! P-wave source
  /*! @ingroup sources */
  void p_src(const multi1d<LatticeColorMatrix>& u, 
	     LatticePropagator& chi,
	     int direction);

  //! P-wave source
  /*! @ingroup sources */
  void p_src(const multi1d<LatticeColorMatrix>& u, 
	     LatticeFermion& chi,
	     int direction);

}  // end namespace Chroma

#endif
