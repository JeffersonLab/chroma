// $Id: lwldslash_w_sse.cc,v 1.21 2004-12-17 17:45:04 bjoo Exp $
/*! \file
 *  \brief Wilson Dslash linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/lwldslash_w_sse.h"

using namespace QDP;

// This is in C++ so it comes outside the extern "C" {}
extern void qdp_pack_gauge(const multi1d<LatticeColorMatrix>&_u, multi1d<PrimitiveSU3Matrix>& u_tmp);

namespace Chroma 
{ 
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



  //! Take deriv of D
  /*! \return Computes   chi^dag * \dot(D} * psi  */
  void 
  SSEWilsonDslash::deriv(multi1d<LatticeColorMatrix>& ds_u,
			 const LatticeFermion& chi, const LatticeFermion& psi, 
			 enum PlusMinus isign, int cb) const
  {
    START_CODE();

    ds_u.resize(Nd);

    LatticeColorMatrix ds_tmp;
    for(int mu = 0; mu < Nd; ++mu)
    {
      LatticeFermion f_tmp;
      switch (isign)
      {
      case PLUS:
	// Undaggered:
        ds_u[mu][rb[cb]]    = u[mu]*traceSpin(outerProduct(shift(psi - Gamma(1 << mu)*psi, FORWARD, mu),chi));	
	ds_u[mu][rb[1-cb]]  = zero;

	// The piece that comes from the U^daggered term. 
	// This piece is just -dagger() of the piece from applying
	// this function on the opposite checkerboard. It therefore
	// only contributes a factor of 2 to the traceless antihermitian
	// part of the result. Which should be swept into the taproj
	// normalisation. Right now until then, I explicitly multiply
	// the result by 0.5 below.

	// ds_u[mu][rb[1-cb]]  = traceSpin(outerProduct(psi + Gamma(1 << mu)*psi,shift(chi, FORWARD, mu)))*adj(u[mu]);
       	// ds_u[mu][rb[1-cb]] *= -Real(1);
	//
	// From factor of 2 that comes from the U^daggered piece.
	// This maybe should be absorbed into the taproj normalisation
	//
	// ds_u[mu] *= Real(0.5);

	break;

      case MINUS:
	// Daggered:
	ds_u[mu][rb[cb]]    = u[mu]*traceSpin(outerProduct(shift(psi + Gamma(1 << mu)*psi, FORWARD, mu),chi));
	
	ds_u[mu][rb[1-cb]] = zero;
	
	// The piece that comes from the U^daggered term. 
	// This piece is just -dagger() of the piece from applying
	// this function on the opposite checkerboard. It therefore
	// only contributes a factor of 2 to the traceless antihermitian
	// part of the result. Which should be swept into the taproj
	// normalisation. Right now until then, I explicitly multiply
	// the result by 0.5 below.
	//
	//        ds_u[mu][rb[1-cb]]  = traceSpin(outerProduct(psi - Gamma(1 << mu)*psi,shift(chi, FORWARD, mu)))*adj(u[mu]);
	//        ds_u[mu][rb[1-cb]] *= -Real(1);
	//	 
	// From factor of 2 that comes from the U^daggered piece.
	// This maybe should be absorbed into the taproj normalisation
	//
	// ds_u[mu] *= Real(0.5);

	break;

      default:
	QDP_error_exit("unknown case");
      }
    }
    
    END_CODE();
  }


}; // End Namespace Chroma

