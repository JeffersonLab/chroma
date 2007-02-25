// -*- C++ -*-
// $Id: diquark_w.h,v 1.2 2007-02-25 22:39:28 edwards Exp $
/*! \file
 *  \brief Construct a diquark object
 */

#ifndef __diquark_w_h__
#define __diquark_w_h__

#include "chromabase.h"
#include "io/qprop_io.h"
#include "util/ft/sftmom.h"

namespace Chroma 
{

  //! Dense QQDiquark object
  struct QQDiquarkContract_t
  {
    multiNd<LatticeComplex>  comp;
  };


  //! Unpack a quark
  /*!
   * \ingroup hadron
   *
   * We need this fast, so at the expense of a lot of memory we will
   * expose all the color/spin indices of each propagator into a temporary
   */
  multi2d< multi2d<LatticeComplex> > unpackQuark(const LatticePropagator& quark_propagator);


  //! Construct a QQ diquark object
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions!
   *
   * In all baryons the colour components are contracted with the totally
   * antisymmetric 'tensor' eps(a,b,c) = antisym_tensor(a,b,c).
   *
   * \param diquark                  diquark object (in real space) ( Write )
   * \param quark_propagator_1       quark propagator ( Read )
   * \param quark_propagator_2       quark propagator ( Read )
   */

  void QQDiquark(QQDiquarkContract_t& diquark, 
		 const LatticePropagator& quark_propagator_1, 
		 const LatticePropagator& quark_propagator_2);

}  // end namespace Chroma

#endif
