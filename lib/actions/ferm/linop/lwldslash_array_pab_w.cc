// $Id: lwldslash_array_pab_w.cc,v 3.8 2008-01-21 20:18:50 edwards Exp $
/*! \file
 *  \brief Wilson Dslash linear operator array
 */

#include "chromabase.h"
#include "chroma_config.h"
// Need this for the refcount namespace
#include "actions/ferm/linop/lwldslash_w_pab.h"
#include "actions/ferm/linop/lwldslash_array_pab_w.h"

#include <wfm.h>

namespace Chroma 
{ 
  //! Empty constructor. Must use create later
  PABWilsonDslashArray::PABWilsonDslashArray()
  {
  }


  //! Full constructor
  PABWilsonDslashArray::PABWilsonDslashArray(Handle< FermState<T,P,Q> > state,
					     int N5_)
  {
    create(state,N5_);
  }


  //! Full constructor
  PABWilsonDslashArray::PABWilsonDslashArray(Handle< FermState<T,P,Q> > state,
					     int N5_,
					     const AnisoParam_t& aniso_)
  {
    create(state,N5_,aniso_);
  }


  //! Creation routine
  void PABWilsonDslashArray::create(Handle< FermState<T,P,Q> > state,
				    int N5_)
  {
    multi1d<Real> cf(Nd);
    cf = 1.0;
    create(state, N5_, cf);
  }


  //! Creation routine with anisotropy
  void PABWilsonDslashArray::create(Handle< FermState<T,P,Q> > state, int N5_,
				    const AnisoParam_t& anisoParam) 
  {
    START_CODE();

    create(state, N5_, makeFermCoeffs(anisoParam));

    END_CODE();
  }

  //! Creation routine
  void PABWilsonDslashArray::create(Handle< FermState<T,P,Q> > state, int N5_,
				    const multi1d<Real>& coeffs_)
  {
    START_CODE();

    N5 = N5_;
    coeffs = coeffs_;

    // Save a copy of the fermbc
    fbc = state->getFermBC();

    // Sanity check
    if (fbc.operator->() == 0)
    {
      QDPIO::cerr << "PABWilsonDslashArray: error: fbc is null" << endl;
      QDP_abort(1);
    }

    // Temporary copy - not kept
    multi1d<LatticeColorMatrix> u = state->getLinks();
  
    // Rescale the u fields by the anisotropy
    for(int mu=0; mu < u.size(); ++mu)
    {
      u[mu] *= coeffs[mu];
    }

    // Rearrange the gauge

    // Allocate the packed gauge
    bool got_fast = true;
    try { 
      packed_gauge=(PrimitiveSU3Matrix *)
         QDP::Allocator::theQDPAllocator::Instance().allocate(
			Nd*sizeof(PrimitiveSU3Matrix)*Layout::sitesOnNode(), 
			QDP::Allocator::FAST);
    }
    catch( std::bad_alloc ) { 
	got_fast = false;
    }

    if( got_fast == false ) { 
      try { 
      packed_gauge=(PrimitiveSU3Matrix *)
         QDP::Allocator::theQDPAllocator::Instance().allocate(
			Nd*sizeof(PrimitiveSU3Matrix)*Layout::sitesOnNode(), 
			QDP::Allocator::DEFAULT);
      }
      catch( std::bad_alloc ) { 
	QDPIO::cerr << "Could not allocate packed gauge in lwldslash_array_pab_w.cc line 67" << endl << flush;
        QDP_abort(1);
      }

    }




    // Rearrange gauge from  Dir,Site, Matrix 
    // to Site, Dir, Matrix
    // (and sites are ordered in cb order)
  
    int volume=Layout::sitesOnNode();
    for(int ix=0; ix < volume; ix++) { 
      for(int mu=0; mu < 4; mu++) { 
	packed_gauge[ 4*ix + mu ] = transpose(u[mu].elem(ix).elem());
      }
    }

    // Need this

    wil_cbsize=Layout::sitesOnNode()/2;
    if ( PABDslashEnv::refcount == 0 ) { 
      // I probably need the local size here...
      const multi1d<int>& lsize= Layout::subgridLattSize();
      const multi1d<int>& machsize = Layout::logicalSize();
      
      // Set up the wilson thingie
      for(int i = 0; i < Nd; i++) { 
	wil.local_latt[i] = lsize[i];
	
	if( machsize[i] == 1 ) {
	  wil.local_comm[i] = 1;

	}
	else {
	  wil.local_comm[i] = 0;
	}
      }

#ifdef CHROMA_USE_SLOPPY_BAGEL_DSLASH
      wil.BGL = true;
      wil.SloppyPrecision = true;
#endif 

      wfm_vec_init(&wil);

#ifdef CHROMA_WFM_NO_END
      // If we never call the end, 
      // the refcount is a lock variable
      // and we never increase it
      PABDslashEnv::refcount = 1;
#endif
      
    }

#ifndef CHROMA_WFM_NO_END
    // Correct refcounting 
    PABDslashEnv::refcount++;
#endif

    END_CODE();
  }

  PABWilsonDslashArray::~PABWilsonDslashArray(void) 
  {

#ifndef CHROMA_WFM_NO_END
    // Only deal with refcount defreases if the refcount is at least 1
    if( PABDslashEnv::refcount > 0 ) {

      // Always decrease the refcount
      PABDslashEnv::refcount--;

      // Only call free when the refcount == 0
      if( PABDslashEnv::refcount == 0 ) { 
	wfm_vec_end(&wil);
      }
    }
#endif 

    // Free the gauge field tho
    QDP::Allocator::theQDPAllocator::Instance().free(packed_gauge);
    
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
  PABWilsonDslashArray::apply (multi1d<LatticeFermion>& chi, 
			       const multi1d<LatticeFermion>& psi, 
			       enum PlusMinus isign, int cb) const
  {
    START_CODE();

    if( chi.size() != N5 ) chi.resize(N5);

    multi1d<Float*> chis(N5);
    multi1d<Float*> psis(N5);
    multi1d<int>    cbs(N5);

    for(int i=0; i < N5; i++) { 
      chis[i]=(Float*)&(chi[i].elem(wil_cbsize*cb).elem(0).elem(0).real());
      psis[i]=(Float*)&(psi[i].elem(wil_cbsize*(1-cb)).elem(0).elem(0).real());
      cbs[i]=1-cb;
    }

    wfm_dslash_vec(N5, 
		   (Float**)&(chis[0]), 
		   (Float *)&(packed_gauge[0]),
		   (Float**)&(psis[0]),
		   (int *)&(cbs[0]),
		   (isign == (enum PlusMinus)PLUS ? 0 : 1));
		     

    END_CODE();
  }

  //! Apply Wilson-Dirac dslash
  /*! \ingroup linop
   * Wilson dslash
   *
   * Arguments:
   *
   *  \param chi	      Result                           (Write)
   *  \param psi	      Pseudofermion field              (Read)
   *  \param isign      D'^dag or D' ( MINUS | PLUS ) resp.    (Read)
   *  \param cb	      Checkerboard of OUTPUT vector            (Read) 
   */
  void
  PABWilsonDslashArray::apply (LatticeFermion& chi, const LatticeFermion& psi, 
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
    /* sse_su3dslash_wilson((SSEREAL *)&(packed_gauge[0]),
       (SSEREAL *)&(psi.elem(0).elem(0).elem(0).real()),
       (SSEREAL *)&(chi.elem(0).elem(0).elem(0).real()),
       (int)isign, (1-cb));
    */

    wfm_dslash((Float *)&(chi.elem(wil_cbsize*(cb)).elem(0).elem(0).real()),
	       		  (Float *)&(packed_gauge[0]),
	       		  (Float *)&(psi.elem(wil_cbsize*(1-cb)).elem(0).elem(0).real()),
	       		  1-cb,
	       		  (isign == (enum PlusMinus)PLUS ? 0 : 1));

    getFermBC().modifyF(chi, QDP::rb[cb]);
	     
    END_CODE();
  }


} // End Namespace Chroma

