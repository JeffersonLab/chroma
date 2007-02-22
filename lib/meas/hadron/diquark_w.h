// -*- C++ -*-
// $Id: diquark_w.h,v 1.1 2007-02-22 06:58:55 edwards Exp $
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
    //! Serialize generalized object
    multi1d<ComplexF> serialize();

    multiNd<LatticeComplex>  comp;
  };


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
