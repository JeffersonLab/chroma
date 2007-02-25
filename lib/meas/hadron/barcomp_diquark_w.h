// -*- C++ -*-
// $Id: barcomp_diquark_w.h,v 1.1 2007-02-25 22:39:03 edwards Exp $
/*! \file
 *  \brief Construct all components of a baryon propagator using a diquark
 */

#ifndef __barcomp_diquark_w_h__
#define __barcomp_diquark_w_h__

#include "meas/hadron/diquark_w.h"
#include "meas/hadron/barcomp_w.h"

namespace Chroma 
{

  //! Construct some components of a baryon propagator
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions!
   *
   * In all baryons the colour components are contracted with the totally
   * antisymmetric 'tensor' eps(a,b,c) = antisym_tensor(a,b,c).
   *
   * \param barprop                  baryon correlation function (in real space) ( Write )
   * \param diquark                  diquark ( Read )
   * \param quark_propagator_3       quark propagator ( Read )
   * \param spin_indices             holds list of source/sink spin indices ( Read )
   * \param phases                   object holds list of momenta ( Read )
   * \param t0                       coordinates of source in decay direction ( Read )
   * \param bc_spec                  boundary condition for spectroscopy ( Read )
   */

  void barcompDiquarkSparse(QQQSparse_t& barprop, 
			    const QQDiquarkContract_t& diquark,
			    const LatticePropagator& quark_propagator_3,
			    const multi1d<QQQSpinIndices_t> spin_indices,
			    const SftMom& phases,
			    int t0, int bc_spec);


  //! Construct all components of a baryon propagator
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions!
   *
   * In all baryons the colour components are contracted with the totally
   * antisymmetric 'tensor' eps(a,b,c) = antisym_tensor(a,b,c).
   *
   * \param barprop                  baryon correlation function (in real space) ( Write )
   * \param diquark                  diquark ( Read )
   * \param quark_propagator_3       quark propagator ( Read )
   * \param phases                   object holds list of momenta ( Read )
   * \param t0                       coordinates of source in decay direction ( Read )
   * \param bc_spec                  boundary condition for spectroscopy ( Read )
   */

  void barcompDiquarkDense(QQQDense_t& barprop,
			   const QQDiquarkContract_t& diquark,
			   const LatticePropagator& quark_propagator_3,
			   const SftMom& phases,
			   int t0, int bc_spec);

}  // end namespace Chroma

#endif
