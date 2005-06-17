// $Id: lwldslash_w_pab.cc,v 1.15 2005-06-17 15:17:53 bjoo Exp $
/*! \file
 *  \brief Wilson Dslash linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/lwldslash_w_pab.h"

#ifdef QDP_USE_QCDOC
#warning "Using QALLOC to allocate the packed gauge"
#include <qalloc.h>
#endif
#include <wfm.h>


namespace Chroma 
{
  namespace PABDslashEnv { 
    int refcount = 0;
  };
 
  //! Creation routine
  void PABWilsonDslash::create(const multi1d<LatticeColorMatrix>& u_)
  {

    // For now, keep an extra copy
    u = u_;

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
    // This is a dirty hack until we come up with a better way
#ifdef QDP_USE_QCDOC
    packed_gauge=(PrimitiveSU3Matrix *)qalloc(QFAST|QCOMMS, Nd*Layout::sitesOnNode()*sizeof(PrimitiveSU3Matrix));
    if( packed_gauge == 0x0 ) { 
      packed_gauge=(PrimitiveSU3Matrix *)qalloc(QCOMMS, Nd*Layout::sitesOnNode()*sizeof(PrimitiveSU3Matrix));
    }
#else 
    try {
      packed_gauge=new PrimitiveSU3Matrix[Nd*Layout::sitesOnNode()];
    }
    catch( bad_alloc ) { 
      packed_gauge = 0x0;
    }
#endif
    if( packed_gauge == 0x0 ) { 
      QDPIO::cout << "Unable to allocate packed gauge in PABWilsonDslash::create()" << endl;
      QDP_abort(1);
    }

    wil_cbsize=Layout::sitesOnNode()/2;

    // Rearrange gauge from  Dir,Site, Matrix 
    // to Site, Dir, Matrix
    // (and sites are ordered in cb order)
  
    int volume=Layout::sitesOnNode();
    for(int ix=0; ix < volume; ix++) { 
      for(int mu=0; mu < 4; mu++) { 
	packed_gauge[ 4*ix + mu ] = transpose(u_[mu].elem(ix).elem());
      }
    }

    if ( PABDslashEnv::refcount == 0 ) { 
         wfm_vec_init(&wil);
    }
    PABDslashEnv::refcount++;
  }


  PABWilsonDslash::~PABWilsonDslash(void) 
  {
    if( PABDslashEnv::refcount > 0 ) {
	PABDslashEnv::refcount--;
     

        if( PABDslashEnv::refcount == 0 ) {
	  wfm_vec_end(&wil);
        }
    }
#ifdef QDP_USE_QCDOC
    qfree(packed_gauge);
#else
    delete [] packed_gauge;
#endif
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
	     
  
  }

}; // End Namespace Chroma

