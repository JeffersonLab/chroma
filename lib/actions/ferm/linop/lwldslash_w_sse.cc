// $Id: lwldslash_w_sse.cc,v 1.2 2003-09-12 16:21:26 bjoo Exp $
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
  void pack_gauge_field(int, u_mat_array *, u_mat_array *);
  void init_sse_su3dslash(int);
  void free_sse_su3dslash(void);
  void sse_su3dslash_wilson(SSEREAL* u, SSEREAL *psi, SSEREAL *res, int isign, int cb);
  void make_shift_tables(int *shift, int subgrid_vol_cb, int nrow[], int subgrid_cb_nrow[], int bound[], int Nd);

}

typedef PColorMatrix < RComplex <REAL>, Nc > PrimitiveSU3Matrix;

void my_pack_gauge(const multi1d<LatticeColorMatrix>&_u, multi1d<PrimitiveSU3Matrix>& u_tmp)
{
  int ix, mu, cb, row, col;
  multi1d<PrimitiveSU3Matrix> v(8);;
  int volume = Layout::sitesOnNode()/2;
  
  
  for(cb = 0; cb < 2; cb++) { 
    for(ix = 0; ix < volume; ix+=2) {
      for(mu = 0; mu < 4; mu++) { 

	v[2*mu] = _u[mu].elem(ix + volume*cb).elem();
	v[2*mu+1] = _u[mu].elem(ix + 1 + volume*cb).elem();
      }

      for(mu = 0; mu < 4; mu++) {
	u_tmp[ mu + Nd*(ix + volume*cb) ] = v[mu];
	u_tmp[ mu + Nd*(ix + 1 + volume*cb) ] = v[Nd + mu]; 
      }
    }
  }
}

//! Creation routine
void SSEWilsonDslash::create(const multi1d<LatticeColorMatrix>& _u)
{
  // Initialize internal structures for DSLASH
  if( Layout::primaryNode() ) { 
    cout << "Calling init_sse_su3dslash()... " << flush;
  }

  init_sse_su3dslash(Layout::sitesOnNode()/2);

  packed_gauge.resize( Nd * Layout::sitesOnNode() );

  if( Layout::primaryNode() ) { 
    cout << "Done " << endl << flush;
  }

  if( Layout::primaryNode() ) { 
    cout << "Calling pack_gauge_field..." << flush;
  }
 
  
  my_pack_gauge(_u, packed_gauge);
  

  if( Layout::primaryNode() ) { 
    cout << "Done" << endl << flush;
  }
}


SSEWilsonDslash::~SSEWilsonDslash(void) {
  free_sse_su3dslash();
}

//! General Wilson-Dirac dslash
/*! \ingroup linop
 * Wilson dslash
 *
 * Arguments:
 *
 *  \param psi	      Pseudofermion field				(Read)
 *  \param isign      D'^dag or D' ( MINUS | PLUS ) resp.		(Read)
 *  \param cb	      Checkerboard of OUTPUT vector			(Read) 
 */
LatticeFermion SSEWilsonDslash::apply (const LatticeFermion& psi, enum LinOpSign isign, int cb) const
{
  START_CODE("lWlDslash");

  LatticeFermion chi;


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
		       (SSEREAL *)&(psi.elem((1-cb)*Layout::sitesOnNode()/2).elem(0).elem(0).real()),
		       (SSEREAL *)&(chi.elem((cb)*Layout::sitesOnNode()/2).elem(0).elem(0).real()),
		       (int)isign, (1-cb));
  
  return chi;
}

