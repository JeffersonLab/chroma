// $Id: lwldslash_w_sse.cc,v 1.1 2003-09-10 18:15:05 bjoo Exp $
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
  void init_sse_su3dslash(void);
  void free_sse_su3dslash(void);
  void sse_su3dslash_wilson(SSEREAL* u, SSEREAL *psi, SSEREAL *res, int isign, int cb);
  void make_shift_tables(int *shift, int subgrid_vol_cb, int nrow[], int subgrid_cb_nrow[], int bound[], int Nd);
}

//! Creation routine
void SSEWilsonDslash::create(const multi1d<LatticeColorMatrix>& _u)
{
  // Grab space for packed gauge

  multi3d<ColorMatrix> u_tmp(Nd,2,Layout::sitesOnNode()/2);
  for(int m=0; m < _u.size(); ++m)
  {
    multi1d<ColorMatrix> u_tt(Layout::sitesOnNode());
    QDP_extract(u_tt, _u[m], all);
    
    // Use this guy
    for(int i=0; i < u_tt.size(); ++i)
    {
      int cb = i / (Layout::sitesOnNode()/2);
      u_tmp[m][cb][i] = transpose(u_tt[i]);
    }
  }

  myu.resize(Nd, 2, Layout::sitesOnNode()/2);

  pack_gauge_field(Layout::sitesOnNode()/2,
		   (u_mat_array *)&(u_tmp[0][0][0].elem().elem().elem(0,0).real()),
		   (u_mat_array *)&(myu[0][0][0].elem().elem().elem(0,0).real()));

  init_sse_su3dslash();

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

  chi = zero;

  sse_su3dslash_wilson((SSEREAL *)&(myu[0][0][0].elem().elem().elem(0,0).real()),
		       (SSEREAL *)&(chi.elem((1-cb)*Layout::sitesOnNode()/2).elem(0).elem(0).real()),
		       (SSEREAL *)&(psi.elem(cb*Layout::sitesOnNode()/2).elem(0).elem(0).real()),
		       isign, cb);

  return chi;
}

