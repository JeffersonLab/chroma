// $Id: lwldslash_w_sse.cc,v 1.15 2004-01-27 22:25:16 edwards Exp $
/*! \file
 *  \brief Wilson Dslash linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/lwldslash_w_sse.h"

using namespace QDP;

//! General Wilson-Dirac dslash
/*! \ingroup linop
 * DSLASH
 *
 * This routine is specific to Wilson fermions!
 *
 * Description:
 *
 * This routine applies the operator D' to Psi, putting the result in Chi.
 *
 *	       Nd-1
 *	       ---
 *	       \
 *   chi(x)  :=  >  U  (x) (1 - isign gamma  ) psi(x+mu)
 *	       /    mu			  mu
 *	       ---
 *	       mu=0
 *
 *	             Nd-1
 *	             ---
 *	             \    +
 *                +    >  U  (x-mu) (1 + isign gamma  ) psi(x-mu)
 *	             /    mu			   mu
 *	             ---
 *	             mu=0
 *
 */
extern "C" {
  
  void init_sse_su3dslash(const int* latt_size);
  void free_sse_su3dslash(void);
  void sse_su3dslash_wilson(SSEREAL* u, SSEREAL *psi, SSEREAL *res, int isign, int cb);

}

// This is in C++ so it comes outside the extern "C" {}
extern void qdp_pack_gauge(const multi1d<LatticeColorMatrix>&_u, multi1d<PrimitiveSU3Matrix>& u_tmp);

//! Creation routine
void SSEWilsonDslash::create(const multi1d<LatticeColorMatrix>& _u)
{
  // Initialize internal structures for DSLASH
#if 0
  QDPIO::cout << "Calling init_sse_su3dslash()... " << endl;
#endif

  // Make a copy of the lattice size
  int lat_size[Nd];
  for(int i=0; i < Nd; i++) 
    lat_size[i] = Layout::lattSize()[i];

  init_sse_su3dslash(lat_size);

  packed_gauge.resize( Nd * Layout::sitesOnNode() );

#if 0
  QDPIO::cout << "Done " << endl << flush;

  QDPIO::cout << "Calling pack_gauge_field..." << flush;
#endif

  qdp_pack_gauge(_u, packed_gauge);
  
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
  START_CODE("lWlDslash");

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
  
}

