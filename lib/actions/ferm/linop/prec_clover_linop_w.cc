// $Id: prec_clover_linop_w.cc,v 2.5 2006-02-16 02:24:46 bjoo Exp $
/*! \file
 *  \brief Even-odd preconditioned clover linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/prec_clover_linop_w.h"

using namespace QDP::Hints;

namespace Chroma 
{ 


  //! Creation routine with Anisotropy
  /*!
   * \param u_ 	    gauge field     	       (Read)
   * \param param_  fermion kappa   	       (Read)
   */
  void EvenOddPrecCloverLinOp::create(const multi1d<LatticeColorMatrix>& u_, 
				      const CloverFermActParams& param_)
  {
    // QDPIO::cout << __PRETTY_FUNCTION__ << ": enter" << endl;

    param = param_;

    clov.create(u_, param);
 
    invclov = clov;  // make a copy
    invclov.choles(0);  // invert the cb=0 part

#if 0
    for(int mu=0; mu < Nd; mu++) { 

      for(int nu = mu+1; nu < Nd; nu++) { 

	int mu_nu_index = (1 << mu) + (1 << nu); // 2^{mu} 2^{nu}
	LatticeColorMatrix sigma_XY_dag=zero;
	
	// Need Sigma On Both Checkerboards
	invclov.triacntr(sigma_XY_dag, mu_nu_index, 0);
	
	std::ostringstream filename;
	filename << "triacntr" << mu_nu_index;
	XMLFileWriter outxml(filename.str());
	push(outxml, "root");
	write(outxml, "mu_nu_index", mu_nu_index);
	write(outxml, "sigma_TrAinv", sigma_XY_dag);
	pop(outxml);
	
	outxml.close();
      }
    }
#endif

    D.create(u_, param.anisoParam);
    // QDPIO::cout << __PRETTY_FUNCTION__ << ": exit" << endl;
  }

  //! Apply the the odd-odd block onto a source vector
  void 
  EvenOddPrecCloverLinOp::oddOddLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
				      enum PlusMinus isign) const
  {
    clov.apply(chi, psi, isign, 1);
  }


  //! Apply the the even-even block onto a source vector
  void 
  EvenOddPrecCloverLinOp::evenEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
					enum PlusMinus isign) const
  {
    // Nuke for testing
    clov.apply(chi, psi, isign, 0);
  }

  //! Apply the inverse of the even-even block onto a source vector
  void 
  EvenOddPrecCloverLinOp::evenEvenInvLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
					   enum PlusMinus isign) const
  {
    invclov.apply(chi, psi, isign, 0);
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
  
    //  tmp1_o  =  D_oe   A^(-1)_ee  D_eo  psi_o
    D.apply(tmp1, psi, isign, 0);
    invclov.apply(tmp2, tmp1, isign, 0);
    D.apply(tmp1, tmp2, isign, 1);

    //  chi_o  =  A_oo  psi_o  -  tmp1_o
    clov.apply(chi, psi, isign, 1);


    chi[rb[1]] += mquarter*tmp1;
  }


  //! Apply the even-even block onto a source vector
  void 
  EvenOddPrecCloverLinOp::derivEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
					     const LatticeFermion& chi, const LatticeFermion& psi, 
					     enum PlusMinus isign) const
  {
    clov.deriv(ds_u, chi, psi, isign, 0);
  }

  //! Apply the even-even block onto a source vector
  void 
  EvenOddPrecCloverLinOp::derivEvenEvenLogDet(multi1d<LatticeColorMatrix>& ds_u,
					      enum PlusMinus isign) const
  {
    
    // Testing Odd Odd Term - get nothing from even even term
    invclov.derivTrLn(ds_u, isign, 0);
  }

  //! Apply the the even-odd block onto a source vector
  void 
  EvenOddPrecCloverLinOp::derivEvenOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
					    const LatticeFermion& chi, const LatticeFermion& psi, 
					    enum PlusMinus isign) const
  {
    START_CODE();
    ds_u.resize(Nd);
    D.deriv(ds_u, chi, psi, isign, 0);
    for(int mu=0; mu < Nd; mu++) { 
      ds_u[mu]  *= Real(-0.5);
    }
    END_CODE();
  }
 
  //! Apply the the odd-even block onto a source vector
  void 
  EvenOddPrecCloverLinOp::derivOddEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
					    const LatticeFermion& chi, const LatticeFermion& psi, 
					    enum PlusMinus isign) const
  {

    START_CODE();
    ds_u.resize(Nd);

    D.deriv(ds_u, chi, psi, isign, 1);
    for(int mu=0; mu < Nd; mu++) { 
     ds_u[mu]  *= Real(-0.5);
    }
    END_CODE();
  }

  // Inherit this
  //! Apply the the odd-odd block onto a source vector
  void 
  EvenOddPrecCloverLinOp::derivOddOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
					   const LatticeFermion& chi, const LatticeFermion& psi, 
					   enum PlusMinus isign) const
  {   
    clov.deriv(ds_u, chi, psi, isign, 1);
  }


  //! Return flops performed by the operator()
  unsigned long EvenOddPrecCloverLinOp::nFlops() const
  {
    // DO NOT REMEMBER THE CORRECT VALUE HERE!!!
    return 0;    // NEED CORRECT VALUE
  }

  //! Get the log det of the even even part
  // BUt for now, return zero for testing.
  Double EvenOddPrecCloverLinOp::LogDetEvenEven(void) const  {
    return invclov.cholesDet(0);
  }
} // End Namespace Chroma
