// $Id: unprec_clover_linop_w.cc,v 2.2 2005-12-18 23:53:26 edwards Exp $
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
    param = param_;

    A.create(u_, param);
    D.create(u_, param.anisoParam);
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
    Real mquarter = -0.25;

    //  chi   =  A . psi - 0.25 * D' . psi  */
    A(chi, psi, isign);
    D(tmp, psi, isign);
    chi -= mquarter * tmp;
  }


  void 
  UnprecCloverLinOp::deriv(multi1d<LatticeColorMatrix>& ds_u,
			   const LatticeFermion& chi, const LatticeFermion& psi, 
			   enum PlusMinus isign) const
  {
    QDPIO::cerr << "Unprec clover deriv not implemented" << endl;
    QDP_abort(1);
  }


  //! Return flops performed by the operator()
  unsigned long UnprecCloverLinOp::nFlops() const
  {
    // DO NOT REMEMBER THE CORRECT VALUE HERE!!!
    return 0;    // NEED CORRECT VALUE
  }

} // End Namespace Chroma
