// $Id: lwldslash_w_sse.cc,v 1.22 2004-12-20 03:59:31 edwards Exp $
/*! \file
 *  \brief Wilson Dslash linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/lwldslash_w_sse.h"
#include <sse_config.h>

using namespace QDP;

// This is in C++ so it comes outside the extern "C" {}
extern void qdp_pack_gauge(const multi1d<LatticeColorMatrix>&_u, multi1d<PrimitiveSU3Matrix>& u_tmp);

namespace Chroma 
{ 
  extern "C" 
  {
    void init_sse_su3dslash(const int* latt_size);
    void free_sse_su3dslash(void);
    void sse_su3dslash_wilson(SSEREAL* u, SSEREAL *psi, SSEREAL *res, int isign, int cb);
  }

  //! Creation routine
  void SSEWilsonDslash::create(const multi1d<LatticeColorMatrix>& u_)
  {
    // Initialize internal structures for DSLASH
#if 0
    QDPIO::cout << "Calling init_sse_su3dslash()... " << endl;
#endif

    // Save a copy of the original fields
    u = u_;

    // Initialize using the total problem size
    init_sse_su3dslash(Layout::lattSize().slice());

    // Pack the gauge fields
    packed_gauge.resize( Nd * Layout::sitesOnNode() );

#if 0
    QDPIO::cout << "Done " << endl << flush;

    QDPIO::cout << "Calling pack_gauge_field..." << flush;
#endif

    qdp_pack_gauge(u_, packed_gauge);
  
#if 0
    QDPIO::cout << "Done" << endl << flush;
#endif
  }


  SSEWilsonDslash::~SSEWilsonDslash(void) 
  {
#if 0
    QDPIO::cout << "Calling free_sse_su3dslash()... " << endl;
#endif

    free_sse_su3dslash();
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
  SSEWilsonDslash::apply (LatticeFermion& chi, const LatticeFermion& psi, 
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
  

    END_CODE();
  }

}; // End Namespace Chroma

