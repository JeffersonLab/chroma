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
    Q poly = u[mu];

    for(int n = 1; n < Layout::lattSize()[mu]; ++n)  // run over all links in mu dir
    {    
      Q tmp = shift(poly, FORWARD, mu);
      poly = u[mu] * tmp;
    }

    /* Take the trace and sum up */
    poly_loop = sum(trace(poly)) / Double(Nc*Layout::vol());

    END_CODE();
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


  // ***** Enter boilerplate code hell *******
  // Isn't this what templates are for?

#if QDP_NC!=3
  void polylp(const multi1d<LatticeColorMatrix>& u, DComplex& poly_loop, int mu) 
  {
    polylp_t( u, poly_loop, mu);
  }
#endif
  void polylp(const multi1d<LatticeColorMatrixF3>& u, DComplex& poly_loop, int mu) 
  {
    if(Nc!=3) QDP_abort(1);
    polylp_t( u, poly_loop, mu);
  }
  void polylp(const multi1d<LatticeColorMatrixD3>& u, DComplex& poly_loop, int mu)
  {
    if(Nc!=3) QDP_abort(1);
    polylp_t(u, poly_loop, mu);
  }

#if QDP_NC!=3
  void polylp(const multi1d<LatticeColorMatrix>& u, multi1d<DComplex>& poly_loop)
  { 
      polylp_t(u,poly_loop);
  }
#endif
  void polylp(const multi1d<LatticeColorMatrixF3>& u, multi1d<DComplex>& poly_loop)
  { 
    if(Nc!=3) QDP_abort(1);
    polylp_t(u,poly_loop);
  }
  void polylp(const multi1d<LatticeColorMatrixD3>& u, multi1d<DComplex>& poly_loop)
  {
    if(Nc!=3) QDP_abort(1);
    polylp_t(u,poly_loop);
  }
}  // end namespace Chroma
