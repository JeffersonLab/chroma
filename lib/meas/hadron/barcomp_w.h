// -*- C++ -*-
// $Id: barcomp_w.h,v 3.2 2007-02-02 05:22:47 edwards Exp $
/*! \file
 *  \brief Construct all components of a baryon propagator
 */

#ifndef __barcomp_w_h__
#define __barcomp_w_h__

#include "chromabase.h"
#include "io/qprop_io.h"
#include "util/ft/sftmom.h"

namespace Chroma 
{

  //! Sparse QQQ object
  struct QQQSparse_t
  {
    //! Serialize generalized object
    multi1d<ComplexF> serialize();

    //! QQQ spin element
    struct QQQElem_t
    {
      QQQSpinIndices_t  source_sink;
      multi1d<ComplexF> corr;
    };

    int                 length;
    multi1d<QQQElem_t>  corrs;
  };


  //! Dense QQQ object
  struct QQQDense_t
  {
    //! Serialize generalized object
    multi1d<ComplexF> serialize();

    int                 length;
    multiNd<Complex>    corrs;
  };


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
   * \param quark_propagator_1       quark propagator ( Read )
   * \param quark_propagator_2       quark propagator ( Read )
   * \param quark_propagator_3       quark propagator ( Read )
   * \param spin_indices             holds list of source/sink spin indices ( Read )
   * \param phases                   object holds list of momenta ( Read )
   * \param t0                       coordinates of source in decay direction ( Read )
   * \param bc_spec                  boundary condition for spectroscopy ( Read )
   */

  void barcompSparse(QQQSparse_t& barprop, 
		     const LatticePropagator& quark_propagator_1, 
		     const LatticePropagator& quark_propagator_2,
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
   * \param quark_propagator_1       quark propagator ( Read )
   * \param quark_propagator_2       quark propagator ( Read )
   * \param quark_propagator_3       quark propagator ( Read )
   * \param phases                   object holds list of momenta ( Read )
   * \param t0                       coordinates of source in decay direction ( Read )
   * \param bc_spec                  boundary condition for spectroscopy ( Read )
   */

  void barcomp(QQQDense_t& barprop, 
	       const LatticePropagator& quark_propagator_1, 
	       const LatticePropagator& quark_propagator_2,
	       const LatticePropagator& quark_propagator_3,
	       const SftMom& phases,
	       int t0, int bc_spec);

}  // end namespace Chroma

#endif
