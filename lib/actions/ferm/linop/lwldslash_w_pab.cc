// $Id: lwldslash_w_pab.cc,v 3.0 2006-04-03 04:58:50 edwards Exp $
/*! \file
 *  \brief Wilson Dslash linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/lwldslash_w_pab.h"

#include <wfm.h>


namespace Chroma 
{
  namespace PABDslashEnv { 
    int refcount = 0;
  };
 
  //! Empty constructor
  PABWilsonDslash::PABWilsonDslash() {}
  
  //! Full constructor
  PABWilsonDslash::PABWilsonDslash(Handle< FermState<T,P,Q> > state)
  {
    create(state);
  }
  
  //! Full constructor with anisotropy
  PABWilsonDslash::PABWilsonDslash(Handle< FermState<T,P,Q> > state,
				   const AnisoParam_t& aniso_) 
  {
    create(state, aniso_);
  }

  //! Creation routine
  void PABWilsonDslash::create(Handle< FermState<T,P,Q> > state)
  {
    AnisoParam_t foo;
    create(state, foo);
  }

  //! Creation routine with anisotropy
  void PABWilsonDslash::create(Handle< FermState<T,P,Q> > state,
			       const AnisoParam_t& aniso_) 
  {
    START_CODE();

    anisoParam = aniso_;
    fbc = state->getFermBC();

    // Sanity check
    if (fbc.operator->() == 0)
    {
      QDPIO::cerr << "PABWilsonDslash: error: fbc is null" << endl;
      QDP_abort(1);
    }

    // Fold in anisotropy
    multi1d<LatticeColorMatrix> u = u_;
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

    // Rearrange the gauge
    // Allocate the packed gauge
    // Now use the QDP allocator system which calls qalloc on QCDOC
    bool got_fast=true;
    try { 
      packed_gauge=(PrimitiveSU3Matrix *)
	QDP::Allocator::theQDPAllocator::Instance().allocate(
			Nd*sizeof(PrimitiveSU3Matrix)*Layout::sitesOnNode(), 
			QDP::Allocator::FAST);
    }
    catch(std::bad_alloc ) { 
      got_fast = false;
    }

    if( got_fast == false ) { 
      try { 
	packed_gauge=(PrimitiveSU3Matrix *)
	  QDP::Allocator::theQDPAllocator::Instance().allocate(
			Nd*sizeof(PrimitiveSU3Matrix)*Layout::sitesOnNode(), 
			QDP::Allocator::DEFAULT);
      }
      catch(std::bad_alloc ) { 
	QDPIO::cerr << "Unable to allocate memory for packed array in lwldslash_w_pab.cc line 65" << endl << flush;
	QDP_abort(1);
      }
    }
    
    wil_cbsize=Layout::sitesOnNode()/2;

    // Rearrange gauge from  Dir,Site, Matrix 
    // to Site, Dir, Matrix
    // (and sites are ordered in cb order)
  
    int volume=Layout::sitesOnNode();
    for(int ix=0; ix < volume; ix++) { 
      for(int mu=0; mu < 4; mu++) { 
	packed_gauge[4*ix + mu] = transpose(u[mu].elem(ix).elem());
      }
    }

    if ( PABDslashEnv::refcount == 0 ) { 
         wfm_vec_init(&wil);
    }
    PABDslashEnv::refcount++;

    END_CODE();
  }


  PABWilsonDslash::~PABWilsonDslash() 
  {
    if( PABDslashEnv::refcount > 0 ) 
    {
      PABDslashEnv::refcount--;
    
      if( PABDslashEnv::refcount == 0 ) {
	wfm_vec_end(&wil);
      }
    }

    QDP::Allocator::theQDPAllocator::Instance().free(packed_gauge);
  }

  //! Apply Wilson-Dirac dslash
  /*! \ingroup linop
   * Wilson dslash
  
   * Arguments:
   *
   *  \param chi	      Result                           (Write)
   *  \param psi	      Pseudofermion field              (Read)
   *  \param isign      D'^dag or D' ( MINUS | PLUS ) resp.    (Read)
   *  \param cb	      Checkerboard of OUTPUT vector            (Read) 
   */
  void
  PABWilsonDslash::apply (LatticeFermion& chi, const LatticeFermion& psi, 
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
	     
    getFermBC().modifyF(chi, rb[cb]);

    END_CODE();
  }

}; // End Namespace Chroma

