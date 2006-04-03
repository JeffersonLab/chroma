// $Id: unprec_clover_linop_w.cc,v 3.0 2006-04-03 04:58:51 edwards Exp $
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
   * \param fs 	    gauge field     	       (Read)
   * \param param_  parameters   	       (Read)
   */
  void UnprecCloverLinOp::create(Handle< FermState<T,P,Q> > fs,
				 const CloverFermActParams& param_)
  {
    //   QDPIO::cout << __PRETTY_FUNCTION__ << ": enter" << endl;

    param = param_;

    A.create(fs, param);
    D.create(fs, param.anisoParam);

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

    getFermBC().modifyF(chi);
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
    
    getFermBC().zero(ds_u);
  }


  //! Return flops performed by the operator()
  unsigned long UnprecCloverLinOp::nFlops() const
  {
    unsigned long site_flops = D.nFlops()+A.nFlops()+4*Nc*Ns;
    return site_flops*Layout::sitesOnNode();
  }

} // End Namespace Chroma
