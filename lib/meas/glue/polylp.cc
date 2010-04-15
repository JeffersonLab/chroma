// $Id: polylp.cc,v 3.1 2006-08-24 21:04:31 edwards Exp $
/*! \file
 *  \brief Calculate the global normalized sum of the Polyakov loop
 */

#include "chromabase.h"
#include "meas/glue/polylp.h"

namespace Chroma 
{

  //! Compute Polyakov loop
  /*!
   * \ingroup glue
   *
   * \param u          gauge field (Read)
   * \param poly_loop  Polyakov loop average in direction mu (Write) 
   * \param mu         direction of Polyakov loop (Read)
   */
  template<typename Q> 
  void polylp_t(const multi1d<Q>& u, DComplex& poly_loop, int mu)
  {
    START_CODE();
        
    // Initial Polyakov loop
    LatticeColorMatrix poly = u[mu];

    for(int n = 1; n < Layout::lattSize()[mu]; ++n)  // run over all links in mu dir
    {    
      LatticeColorMatrix tmp = shift(poly, FORWARD, mu);
      poly = u[mu] * tmp;
    }

    /* Take the trace and sum up */
    poly_loop = sum(trace(poly)) / Double(Nc*Layout::vol());

    END_CODE();
  }


  void polylp(const multi1d<LatticeColorMatrixF3>& u, DComplex& poly_loop, int mu) 
  {
     polylp_t( u, poly_loop, mu);
  }

  void polylp(const multi1d<LatticeColorMatrixD3>& u, DComplex& poly_loop, int mu)
  {
     polylp_t(u, poly_loop, mu);
  }

  //! Compute Polyakov loop
  /*!
   * \ingroup glue
   *
   * \param u          gauge field (Read)
   * \param poly_loop  Polyakov loop average (Write) 
   */
  template<typename Q>
  void polylp_t(const multi1d<Q>& u, multi1d<DComplex>& poly_loop)
  {
    START_CODE();
 
    poly_loop.resize(Nd);

    for(int mu = 0; mu < Nd; ++mu)
      polylp(u, poly_loop[mu], mu);

    END_CODE();
  }


  void polylp(const multi1d<LatticeColorMatrixF3>& u, multi1d<DComplex>& poly_loop)
  { 
      polylp_t(u,poly_loop);
  }

  void polylp(const multi1d<LatticeColorMatrixD3>& u, multi1d<DComplex>& poly_loop)
  {
      polylp_t(u,poly_loop);
  }
}  // end namespace Chroma
