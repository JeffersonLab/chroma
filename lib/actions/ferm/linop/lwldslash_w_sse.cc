// $Id: lwldslash_w_sse.cc,v 1.13 2003-12-15 21:46:48 edwards Exp $
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
  
  void init_sse_su3dslash(int *);
  void free_sse_su3dslash(void);
  void sse_su3dslash_wilson(SSEREAL* u, SSEREAL *psi, SSEREAL *res, int isign, int cb);
  void make_shift_tables(int *shift, int subgrid_vol_cb, int nrow[], int subgrid_cb_nrow[], int bound[], int Nd);

}

// This is in C++ so it comes outside the extern "C" {}
extern void qdp_pack_gauge(const multi1d<LatticeColorMatrix>&_u, multi1d<PrimitiveSU3Matrix>& u_tmp);

//! Creation routine
void SSEWilsonDslash::create(const multi1d<LatticeColorMatrix>& _u)
{
  // Initialize internal structures for DSLASH
#if 1
  QDPIO::cout << "Calling init_sse_su3dslash()... " << endl;
#endif

  const multi1d<int>& subgrid_size = Layout::subgridLattSize();

  // Subgrid size after checkerboarding
  int s_size[4];

  int i;
  for(i=0 ; i < Nd; i++) { 
    if ( subgrid_size[i] % 2 != 0 ) { 
      QDP_error_exit("This SSE Dslash does not work for odd sublattice. Here the sublattice is odd in dimension %d with length %d\n", i, subgrid_size[i]);
    }
    s_size[i]= subgrid_size[i];
  }

  s_size[0] /= 2;

  init_sse_su3dslash(s_size);

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
#if 1
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
   * SZIN usually passes 1 cb worth of fermions. Instead I have to 
   * Pass the right part of a full fermion.  CB is CB of output 
   * so pass cb*Layout::sitesOnNode()/2 for output and (1-cb)*Layout::sitesOnNode()/2 for input 
   *
   */
  sse_su3dslash_wilson((SSEREAL *)&(packed_gauge[0]),
		       (SSEREAL *)&(psi.elem((1-cb)*(Layout::sitesOnNode()/2)).elem(0).elem(0).real()),
		       (SSEREAL *)&(chi.elem(cb*(Layout::sitesOnNode()/2)).elem(0).elem(0).real()),
		       (int)isign, (1-cb));
  
}

