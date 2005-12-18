// $Id: prec_clover_linop_w.cc,v 2.2 2005-12-18 23:53:26 edwards Exp $
/*! \file
 *  \brief Even-odd preconditioned clover linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/prec_clover_linop_w.h"

using namespace QDP::Hints;

namespace Chroma 
{ 
  // Anonymous namespace
  namespace
  {
    //! Hack for now
    void 
    notImplemented(multi1d<LatticeColorMatrix>& ds_u, 
		   const LatticeFermion& chi, const LatticeFermion& psi, 
		   enum PlusMinus isign)
    {
      QDPIO::cerr << "Prec clover deriv not implemented" << endl;
      QDP_abort(1);
    }

  }


  //! Creation routine with Anisotropy
  /*!
   * \param u_ 	    gauge field     	       (Read)
   * \param param_  fermion kappa   	       (Read)
   */
  void EvenOddPrecCloverLinOp::create(const multi1d<LatticeColorMatrix>& u_, 
				      const CloverFermActParams& param_)
  {
    param = param_;

    clov.create(u_, param);
    invclov = clov;  // make a copy
    invclov.choles(0);  // invert the cb=0 part
    
    D.create(u_, param.anisoParam);
  }


  //! Apply the the even-even block onto a source vector
  void 
  EvenOddPrecCloverLinOp::evenEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
					enum PlusMinus isign) const
  {
    clov.apply(chi, psi, isign, 0);
  }

  //! Apply the inverse of the even-even block onto a source vector
  void 
  EvenOddPrecCloverLinOp::evenEvenInvLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
					   enum PlusMinus isign) const
  {
    invclov.apply(chi, psi, isign, 0);
  }
  
  //! Apply the the odd-odd block onto a source vector
  void 
  EvenOddPrecCloverLinOp::oddOddLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
				      enum PlusMinus isign) const
  {
    invclov.apply(chi, psi, isign, 1);
  }

  //! Apply even-odd linop component
  /*!
   * The operator acts on the entire even sublattice
   *
   * \param chi 	  Pseudofermion field     	       (Write)
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void 
  EvenOddPrecCloverLinOp::evenOddLinOp(LatticeFermion& chi, 
				       const LatticeFermion& psi, 
				       enum PlusMinus isign) const
  {
    START_CODE();

    Real mhalf = -0.5;

    D.apply(chi, psi, isign, 0);
    chi[rb[0]] *= mhalf;
  
    END_CODE();
  }

  //! Apply odd-even linop component
  /*!
   * The operator acts on the entire odd sublattice
   *
   * \param chi 	  Pseudofermion field     	       (Write)
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void 
  EvenOddPrecCloverLinOp::oddEvenLinOp(LatticeFermion& chi, 
				       const LatticeFermion& psi, 
				       enum PlusMinus isign) const
  {
    START_CODE();

    Real mhalf = -0.5;

    D.apply(chi, psi, isign, 1);
    chi[rb[1]] *= mhalf;
  
    END_CODE();
  }


  //! Apply even-odd preconditioned Clover fermion linear operator
  /*!
   * \param chi 	  Pseudofermion field     	       (Write)
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void EvenOddPrecCloverLinOp::operator()(LatticeFermion & chi, 
					  const LatticeFermion& psi, 
					  enum PlusMinus isign) const
  {
    LatticeFermion tmp1; moveToFastMemoryHint(tmp1);
    LatticeFermion tmp2; moveToFastMemoryHint(tmp2);
    Real mquarter = -0.25;
  
    /*  chi_o  =  A_oo  psi_o  -  tmp1_o  */
    invclov.apply(chi, psi, isign, 1);

    //  tmp1_o  =  D_oe   A^(-1)_ee  D_eo  psi_o
    D.apply(tmp1, psi, isign, 0);
    invclov.apply(tmp2, tmp1, isign, 0);
    invclov.apply(tmp1, tmp2, isign, 1);

    chi[rb[1]] += mquarter*tmp1;
  }


  //! Apply the even-even block onto a source vector
  void 
  EvenOddPrecCloverLinOp::derivEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
					     const LatticeFermion& chi, const LatticeFermion& psi, 
					     enum PlusMinus isign) const
  {
    notImplemented(ds_u,chi,psi,isign);
  }

  //! Apply the the even-odd block onto a source vector
  void 
  EvenOddPrecCloverLinOp::derivEvenOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
					    const LatticeFermion& chi, const LatticeFermion& psi, 
					    enum PlusMinus isign) const
  {
    notImplemented(ds_u,chi,psi,isign);
  }
 
  //! Apply the the odd-even block onto a source vector
  void 
  EvenOddPrecCloverLinOp::derivOddEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
					    const LatticeFermion& chi, const LatticeFermion& psi, 
					    enum PlusMinus isign) const
  {
    notImplemented(ds_u,chi,psi,isign);
  }

  //! Apply the the odd-odd block onto a source vector
  void 
  EvenOddPrecCloverLinOp::derivOddOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
					   const LatticeFermion& chi, const LatticeFermion& psi, 
					   enum PlusMinus isign) const
  {
    notImplemented(ds_u,chi,psi,isign);
  }

  //! Override virtual function for efficiency. 
  void 
  EvenOddPrecCloverLinOp::deriv(multi1d<LatticeColorMatrix>& ds_u,
				const LatticeFermion& chi, const LatticeFermion& psi, 
				enum PlusMinus isign) const
  {
    notImplemented(ds_u,chi,psi,isign);
  }


  //! Return flops performed by the operator()
  unsigned long EvenOddPrecCloverLinOp::nFlops() const
  {
    // DO NOT REMEMBER THE CORRECT VALUE HERE!!!
    return 0;    // NEED CORRECT VALUE
  }


} // End Namespace Chroma
