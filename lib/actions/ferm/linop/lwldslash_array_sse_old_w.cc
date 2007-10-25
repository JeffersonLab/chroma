// $Id: lwldslash_array_sse_old_w.cc,v 3.1 2007-10-25 16:10:11 bjoo Exp $
/*! \file
 *  \brief Wilson Dslash linear operator array
 */

#include "chromabase.h"
#include "actions/ferm/linop/lwldslash_array_sse_old_w.h"
#include <sse_config.h>


// This is in C++ so it comes outside the extern "C" {}
extern void qdp_pack_gauge(const multi1d<LatticeColorMatrix>&_u, multi1d<Chroma::PrimitiveSU3Matrix>& u_tmp);

namespace Chroma 
{ 
  extern "C" 
  {
    void init_sse_su3dslash(const int* latt_size);
    void free_sse_su3dslash(void);
    void sse_su3dslash_wilson(SSEREAL* u, SSEREAL *psi, SSEREAL *res, int isign, int cb);
  }

  //! Initialization routine
  void SSEWilsonDslashArray::init()
  {
    START_CODE();

    // Initialize internal structures for DSLASH
#if 0
    QDPIO::cout << "Calling init_sse_su3dslash()... " << endl;
#endif

    // Initialize using the total problem size
    init_sse_su3dslash(Layout::lattSize().slice());

    END_CODE();
  }

  //! Empty constructor. Must use create later
  SSEWilsonDslashArray::SSEWilsonDslashArray() 
  {
    init();
  }

  //! Full constructor
  SSEWilsonDslashArray::SSEWilsonDslashArray(Handle< FermState<T,P,Q> > state,
					     int N5_,
					     const AnisoParam_t& aniso_)
  {
    init(); 
    create(state,N5_,aniso_);
  }


  //! Full constructor
  SSEWilsonDslashArray::SSEWilsonDslashArray(Handle< FermState<T,P,Q> > state,
					     int N5_)
  {
    init(); 
    create(state,N5_);
  }


  //! Creation routine
  void SSEWilsonDslashArray::create(Handle< FermState<T,P,Q> > state,
				    int N5_)
  {
    AnisoParam_t aniso;
    create(state, N5_, aniso);
  }


  //! Creation routine
  void SSEWilsonDslashArray::create(Handle< FermState<T,P,Q> > state,
				    int N5_, const AnisoParam_t& aniso)
  {
    START_CODE();

    N5 = N5_;
    anisoParam = aniso;

    // Save a copy of the fermbc
    fbc = state->getFermBC();

    if (fbc.operator->() == 0)
    {
      QDPIO::cerr << "SSEWilsonDslashArray: error: fbc is null" << endl;
      QDP_abort(1);
    }

    // Temporary copy - not kept
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
    packed_gauge.resize( Nd * Layout::sitesOnNode() );

#if 0
    QDPIO::cout << "Done " << endl << flush;

    QDPIO::cout << "Calling pack_gauge_field..." << flush;
#endif

    qdp_pack_gauge(u, packed_gauge);
  
#if 0
    QDPIO::cout << "Done" << endl << flush;
#endif
    
    END_CODE();
  }


  SSEWilsonDslashArray::~SSEWilsonDslashArray() 
  {
    START_CODE();

#if 0
    QDPIO::cout << "Calling free_sse_su3dslash()... " << endl;
#endif

    free_sse_su3dslash();

    END_CODE();
  }


  //! General Wilson-Dirac dslash
  /*! \ingroup linop
   * Wilson dslash
   *
   * Arguments:
   *
   *  \param chi      Result				                (Write)
   *  \param psi      Pseudofermion field				(Read)
   *  \param isign    D'^dag or D' ( MINUS | PLUS ) resp.		(Read)
   *  \param cb	      Checkerboard of OUTPUT vector			(Read) 
   */
  void 
  SSEWilsonDslashArray::apply (multi1d<LatticeFermion>& chi, 
			       const multi1d<LatticeFermion>& psi, 
			       enum PlusMinus isign, int cb) const
  {
    START_CODE();

    if( chi.size() != N5 ) chi.resize(N5);

    for(int n=0; n < N5; ++n)
      apply(chi[n], psi[n], isign, cb);

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
  SSEWilsonDslashArray::apply (LatticeFermion& chi, const LatticeFermion& psi, 
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
    sse_su3dslash_wilson((SSEREAL *)&(packed_gauge[0]),
			 (SSEREAL *)&(psi.elem(0).elem(0).elem(0).real()),
			 (SSEREAL *)&(chi.elem(0).elem(0).elem(0).real()),
			 (int)isign, (1-cb));
  
    getFermBC().modifyF(chi, QDP::rb[cb]);

    END_CODE();
  }

} // End Namespace Chroma

