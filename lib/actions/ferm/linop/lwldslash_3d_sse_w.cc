// $Id: lwldslash_3d_sse_w.cc,v 3.1 2007-12-04 16:04:42 bjoo Exp $
/*! \file
 *  \brief Wilson Dslash linear operator
 */
#include "qdp_config.h"
#if QDP_NS == 4
#if QDP_ND == 4
#if QDP_NC == 3

#include "chromabase.h"
#include "actions/ferm/linop/lwldslash_3d_sse_w.h"
#include "sse_dslash_3d.h"
#include "sse_dslash_qdp_packer_3d.h"

namespace Chroma 
{ 

  //! Initialization routine
  void SSEWilsonDslash3D::init()
  {
    // Initialize internal structures for DSLASH
#if 0
    QDPIO::cout << "Calling init_sse_su3dslash()... " << endl;
#endif

    // Initialize using the total problem size
    init_sse_su3dslash_3d(Layout::lattSize().slice());
  }


  //! Empty constructor
  SSEWilsonDslash3D::SSEWilsonDslash3D()
  {
    init();
  }
  
  //! Full constructor
  SSEWilsonDslash3D::SSEWilsonDslash3D(Handle< FermState<T,P,Q> > state)
  { 
    init();
    create(state);
  }
  
  //! Full constructor with anisotropy
  SSEWilsonDslash3D::SSEWilsonDslash3D(Handle< FermState<T,P,Q> > state,
				   const AnisoParam_t& aniso_) 
  {
    init();
    create(state, aniso_);
  }

  //! Creation routine
  void SSEWilsonDslash3D::create(Handle< FermState<T,P,Q> > state)
  {
    AnisoParam_t foo;
    create(state, foo);
  }

  //! Creation routine with anisotropy
  void SSEWilsonDslash3D::create(Handle< FermState<T,P,Q> > state,
			       const AnisoParam_t& aniso_) 
  {
    START_CODE();

    // Save a copy of the aniso params original fields and with aniso folded in
    anisoParam = aniso_;

    // Save a copy of the fermbc
    fbc = state->getFermBC();

    // Sanity check
    if (fbc.operator->() == 0)
    {
      QDPIO::cerr << "SSEWilsonDslash3D: error: fbc is null" << endl;
      QDP_abort(1);
    }

    // Fold in anisotropy
    multi1d<LatticeColorMatrix> u = state->getLinks();
    Real ff = where(anisoParam.anisoP, anisoParam.nu / anisoParam.xi_0, Real(1));
  
    if (anisoParam.anisoP)
    {
      // Rescale the u fields by the anisotropy
      for(int mu=0; mu < u.size(); ++mu)
      {
	if (mu != anisoParam.t_dir)
	  u[mu] *= ff;
      }
    }

    // Pack the gauge fields
    packed_gauge.resize( 4 * Layout::sitesOnNode() );
    SSEDslash3D::qdp_pack_gauge_3d(u, packed_gauge);
  


    END_CODE();
  }


  SSEWilsonDslash3D::~SSEWilsonDslash3D() 
  {
    START_CODE();

#if 0
    QDPIO::cout << "Calling free_sse_su3dslash()... " << endl;
#endif

    free_sse_su3dslash_3d();

    END_CODE();
  }

  //! General Wilson-Dirac dslash
  /*! \ingroup linop
   * Wilson dslash
   *
   * Arguments:
   *
   *  \param chi	      Result				                (Write)
   *  \param psi	      Pseudofermion field				(Read)
   *  \param isign      D'^dag or D' ( MINUS | PLUS ) resp.		(Read)
   *  \param cb	      Checkerboard of OUTPUT vector			(Read) 
   */
  void
  SSEWilsonDslash3D::apply (LatticeFermion& chi, const LatticeFermion& psi, 
			  enum PlusMinus isign, int cb) const
  {
    START_CODE();

    /* Pass the right parities. 
     *
     * SZIN standard is that cb is cb of INPUT 
     * Chroma standard is that the cb is cb of OUTPUT
     *
     * Need to invert cb for SZIN style SSE call 
   
     *
     *
     * Pass all the fermion and all the gauge pieces
     *
     * NOTE: this breaks usage from SZIN. However, Chroma and SSE dslash can handle 
     * odd subgrid lattice sizes, whereas SZIN cannot. Thus, I must pass all fermion
     * cbs to support such flexibility.
     *
     */
    int source_cb = 1 - cb;
    int target_cb = cb;
    int cbsites = QDP::Layout::sitesOnNode()/2;

    sse_su3dslash_wilson_3d((SSEREAL *)&(packed_gauge[0]),
			    (SSEREAL *)&(psi.elem(0).elem(0).elem(0).real()),
			    (SSEREAL *)&(chi.elem(0).elem(0).elem(0).real()),
			    (int)isign, source_cb);
  

    getFermBC().modifyF(chi, QDP::rb[cb]);

    END_CODE();
  }

} // End Namespace Chroma

#endif
#endif
#endif
