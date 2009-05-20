// $Id: lwldslash_w_cppd.cc,v 3.3 2009-05-20 19:28:11 bjoo Exp $
/*! \file
 *  \brief Wilson Dslash linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/lwldslash_w_cppd.h"
#include "cpp_dslash.h"
#include "cpp_dslash_qdp_packer.h"

using namespace CPlusPlusWilsonDslash;

namespace Chroma 
{ 
  //! Initialization routine
  void CPPWilsonDslashD::init()
  {
    // Initialize internal structures for DSLASH
#if 0
    QDPIO::cout << "Calling init_sse_su3dslash()... " << endl;
#endif


      // Initialize using the total problem size
      D=new Dslash<double>(Layout::lattSize().slice(),
			Layout::QDPXX_getSiteCoords,
                        Layout::QDPXX_getLinearSiteIndex,
                        Layout::QDPXX_nodeNumber);
      
  }


  //! Empty constructor
  CPPWilsonDslashD::CPPWilsonDslashD()
  {
    init();
  }
  
  //! Full constructor
  CPPWilsonDslashD::CPPWilsonDslashD(Handle< FermState<T,P,Q> > state)
  { 
    init();
    create(state);
  }
  
  //! Full constructor with anisotropy
  CPPWilsonDslashD::CPPWilsonDslashD(Handle< FermState<T,P,Q> > state,
				   const AnisoParam_t& aniso_) 
  {
    init();
    create(state, aniso_);
  }

  //! Full constructor with general coefficients
  CPPWilsonDslashD::CPPWilsonDslashD(Handle< FermState<T,P,Q> > state,
				   const multi1d<Real>& coeffs_)
  {
    init();
    create(state, coeffs_);
  }

  //! Creation routine
  void CPPWilsonDslashD::create(Handle< FermState<T,P,Q> > state)
  {
    multi1d<Real> cf(Nd);
    cf = 1.0;
    create(state, cf);
  }

  //! Creation routine with anisotropy
  void CPPWilsonDslashD::create(Handle< FermState<T,P,Q> > state,
			       const AnisoParam_t& anisoParam) 
  {
    START_CODE();

    create(state, makeFermCoeffs(anisoParam));

    END_CODE();
  }

  //! Full constructor with general coefficients
  void CPPWilsonDslashD::create(Handle< FermState<T,P,Q> > state,
			       const multi1d<Real>& coeffs_)
  {
    START_CODE();

    // Save a copy of the aniso params original fields and with aniso folded in
    coeffs = coeffs_;

    // Save a copy of the fermbc
    fbc = state->getFermBC();

    // Sanity check
    if (fbc.operator->() == 0)
    {
      QDPIO::cerr << "CPPWilsonDslashD: error: fbc is null" << endl;
      QDP_abort(1);
    }

    // Fold in anisotropy
    multi1d<LatticeColorMatrixD> u = state->getLinks();
  
    // Rescale the u fields by the anisotropy
    for(int mu=0; mu < u.size(); ++mu)
    {
      u[mu] *= coeffs[mu];
    }

    // Pack the gauge fields
    packed_gauge.resize( Nd * Layout::sitesOnNode() );

#if 0
    QDPIO::cout << "Done " << endl << flush;

    QDPIO::cout << "Calling pack_gauge_field..." << flush;
#endif

    qdp_pack_gauge(u, packed_gauge);
  
#if 0
    QDPIO::cout << "Done" << endl << flush;
#endif

    END_CODE();
  }


  CPPWilsonDslashD::~CPPWilsonDslashD() 
  {
    START_CODE();

#if 0
    QDPIO::cout << "Calling free_sse_su3dslash()... " << endl;
#endif

    // Never free
    // free_sse_su3dslash();

    END_CODE();
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
  CPPWilsonDslashD::apply (T& chi, const T& psi, 
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
    int source_cb = 1 - cb;
    int target_cb = cb;
    int cbsites = QDP::Layout::sitesOnNode()/2;


    (*D)((double *)&(chi.elem(all.start()).elem(0).elem(0).real()),	  
	 (double *)&(psi.elem(all.start()).elem(0).elem(0).real()),
	 (double *)&(packed_gauge[0]),
	 isign, 
	 source_cb);

    // sse_su3dslash_wilson((SSEREAL *)&(packed_gauge[0]),
    //			 (SSEREAL *)&(psi.elem(0).elem(0).elem(0).real()),
    //			 (SSEREAL *)&(chi.elem(0).elem(0).elem(0).real()),
    //			 (int)isign, source_cb);
  
    getFermBC().modifyF(chi, QDP::rb[cb]);

    END_CODE();
  }

} // End Namespace Chroma

