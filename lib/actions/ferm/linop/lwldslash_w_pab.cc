// $Id: lwldslash_w_pab.cc,v 1.1 2004-03-09 20:43:13 bjoo Exp $
/*! \file
 *  \brief Wilson Dslash linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/lwldslash_w_pab.h"

#include <wfm.h>

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
//! Creation routine
void PABWilsonDslash::create(const multi1d<LatticeColorMatrix>& _u)
{
  // I probably need the local size here...
  const multi1d<int>& lsize= Layout::subgridLattSize();
  const multi1d<int>& machsize = Layout::logicalSize();

  // Set up the wilson thingie
  wil.instruction_reg_num=1;

  for(int i = 0; i < Nd; i++) { 
    wil.local_latt[i] = lsize[i];
    
    if( machsize[i] == 1 ) {
      wil.local_comm[i] = 1;

    }
    else {
      wil.local_comm[i] = 0;
      wil.instruction_reg_num++;
    }
  }

  // Rearrange the gauge
  packed_gauge.resize(Nd * Layout::sitesOnNode());

  wil_cbsize=Layout::sitesOnNode()/2;

  // Rearrange gauge from  Dir,Site, Matrix 
  // to Site, Dir, Matrix
  // (and sites are ordered in cb order)
  
  int volume=Layout::sitesOnNode();
  for(int ix=0; ix < volume; ix++) { 
    for(int mu=0; mu < 4; mu++) { 
      packed_gauge[ 4*ix + mu ] = transpose(_u[mu].elem(ix).elem());
    }
  }


  wfm_init(&wil);
}


PABWilsonDslash::~PABWilsonDslash(void) 
{
  wfm_end(&wil);
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
PABWilsonDslash::apply (LatticeFermion& chi, const LatticeFermion& psi, 
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

