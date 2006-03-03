// $Id: unprec_clover_linop_w.cc,v 2.6 2006-03-03 02:37:39 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned clover linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_clover_linop_w.h"

using namespace QDP::Hints;

namespace Chroma 
{ 
  //! Creation routine with Anisotropy
  /*!
   * \param u_ 	    gauge field     	       (Read)
   * \param param_  parameters   	       (Read)
   */
  void UnprecCloverLinOp::create(const multi1d<LatticeColorMatrix>& u_, 
				 const CloverFermActParams& param_)
  {
    //   QDPIO::cout << __PRETTY_FUNCTION__ << ": enter" << endl;

    param = param_;

    A.create(u_, param);

    if( param.ext_fieldP ) {
      LatticeInteger my_coord = Layout::latticeCoordinate(0);
      LatticeReal re = Real(1);

      // field strength * x coordinate
      LatticeReal im = param.ext_field_strength*LatticeReal(my_coord);
    
      // The actual complex number factor for applying the field
      LatticeComplex factor = cmplx(re,im);

      // Now get the modified field
      multi1d<LatticeColorMatrix> u_prime(Nd);

      for(int mu=0; mu < Nd; mu++) {
	
	// X dependent factor on Y links (dim =1) for field in Z direction
	if ( mu == 1 ) {

	  u_prime[mu] = factor*u_[mu];
	}
	else {
	  // Otherwise leave field alone for other directiosn
	  u_prime[mu] = u_[mu];
	}
      }
      D.create(u_prime, param.anisoParam);
    }
    else { 
      D.create(u_, param.anisoParam);
    }

    // QDPIO::cout << __PRETTY_FUNCTION__ << ": exit" << endl;
  }



  //! Apply unpreconditioned Clover fermion linear operator
  /*!
   * The operator acts on the entire lattice
   *
   * \param chi 	  Pseudofermion field          (Write)
   * \param psi 	  Pseudofermion field          (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void UnprecCloverLinOp::operator()(LatticeFermion & chi, 
				     const LatticeFermion& psi, 
				     enum PlusMinus isign) const
  {
    LatticeFermion tmp; moveToFastMemoryHint(tmp);
    Real mhalf = -0.5;

    //  chi   =  A . psi - 0.5 * D' . psi  */
    A(chi, psi, isign); 
    D(tmp, psi, isign);
    chi += mhalf * tmp;
  }


  void 
  UnprecCloverLinOp::deriv(multi1d<LatticeColorMatrix>& ds_u,
			   const LatticeFermion& chi, const LatticeFermion& psi, 
			   enum PlusMinus isign) const
  {
    // A. deriv will resize
    
    A.deriv(ds_u, chi, psi, isign);

    multi1d<LatticeColorMatrix> ds_tmp(Nd);

    ds_tmp = zero;
    D.deriv(ds_tmp, chi, psi, isign);
    for(int mu=0; mu < Nd; mu++) { 
      ds_u[mu] -= Real(0.5)*ds_tmp[mu];
    }
    
  }


  //! Return flops performed by the operator()
  unsigned long UnprecCloverLinOp::nFlops() const
  {
    unsigned long site_flops = D.nFlops()+clov.nFlops()+4*Nc*Ns;
    return site_flops*Layout::sitesOnNode();
  }

} // End Namespace Chroma
