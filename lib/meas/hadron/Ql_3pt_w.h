// -*- C++ -*-
//$Id: Ql_3pt_w.h,v 1.1 2009-01-16 17:30:31 caubin Exp $
/*! \file
 *  \brief Static-Light 3pt function
 */

#ifndef __Ql_3pt_h__
#define __Ql_3pt_h__

#include "chromabase.h"
#include "util/ft/sftmom.h"

namespace Chroma 
{

  //! Heavy-light 3-pt function
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions! 
   *
   * The heavy quark is inserted in the infinitely heavy quark limit
   * by a Wilson-Line without spin indices.
   * We are effectively propagating a spin-1/2 light quark
   * This generates with two quark propagators all three-point
   * functions, with all 16 gamma matrix insertions.
   *
   * \param u                  gauge field (Read) 
   * \param quark_propagator1   quark propagator1 ( Read )
   * \param quark_propagator2   quark propagator2 ( Read )
   * \param src_coord          cartesian coordinates of the source ( Read )
   * \param snk_coord          cartesian coordinates of the sink ( Read )
   * \param phases             object holds list of momenta and Fourier phases ( Read )
   * \param xml                xml file object ( Read )
   * \param xml_group          group name for xml data ( Read )
   *
   */
  void QlQl(const multi1d<LatticeColorMatrix>& u, 
	     const LatticePropagator& quark_propagator1,
	     const LatticePropagator& quark_propagator2,
	     const multi1d<int>& src_coord, 
	     const multi1d<int>& snk_coord, 
	     const SftMom& phases,
	     XMLWriter& xml,
	     const string& xml_group);

}
#endif
