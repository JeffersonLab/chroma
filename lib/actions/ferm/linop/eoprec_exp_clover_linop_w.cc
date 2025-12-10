/*! \file
 *  \brief Even-odd preconditioned exponentiated clover linear operator
 */

#include "actions/ferm/linop/eoprec_exp_clover_linop_w.h"



namespace Chroma 
{ 

   using namespace QDP::Hints;

 //! Creation routine with Anisotropy
  /*!
   * \param u_ 	    gauge field     	       (Read)
   * \param param_  fermion kappa   	       (Read)
   */
  void EvenOddPrecExpCloverLinOp::create(Handle< FermState<T,P,Q> > fs, 
				      const CloverFermActParams& param_)
  {
    START_CODE();
    // QDPIO::cout << __PRETTY_FUNCTION__ << ": enter" << std::endl;

    QDPIO::cout << "Using even-odd preconditioned exponentiated clover\n";

    param = param_;

    clov.create(fs, param);


    //clov. makeExpClov(PLUS,0);
    //makeExpClov(PLUS,1);

 
    invclov.create(fs,param,clov);  // make a copy

#if 0
    invclov.choles(0);  // invert the cb=0 part
#else
        invclov.makeExpClov(PLUS,0,0);
        invclov.makeExpClov(PLUS,1,0);
        
        invclov.makeExpClov(MINUS,0,1);
        invclov.makeExpClov(MINUS,1,1);
#endif


    D.create(fs, param.anisoParam);

    clov_deriv_time = 0;
    clov_apply_time = 0;

    moveToFastMemoryHint(tmp1);
    moveToFastMemoryHint(tmp2);
     
    // QDPIO::cout << __PRETTY_FUNCTION__ << ": exit" << std::endl;
    END_CODE();
  }

  //! Apply the the odd-odd block onto a source std::vector
  void 
  EvenOddPrecExpCloverLinOp::oddOddLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
				      enum PlusMinus isign) const
  {
    START_CODE();

    swatch.reset(); swatch.start();
    clov.apply(chi, psi, isign, 1);
    chi *= (Real(Nd) + param.Mass);

    swatch.stop();
    clov_apply_time += swatch.getTimeInSeconds();

    END_CODE();
  }


  //! Apply the the even-even block onto a source std::vector
  void 
  EvenOddPrecExpCloverLinOp::evenEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
					enum PlusMinus isign) const
  {
    START_CODE();

    // Nuke for testing
    swatch.reset(); swatch.start();
    clov.apply(chi, psi, isign, 0);
    chi *= (Real(Nd) + param.Mass);

    swatch.stop();
    clov_apply_time += swatch.getTimeInSeconds();
    
    END_CODE();
  }

  //! Apply the inverse of the even-even block onto a source std::vector
  void 
  EvenOddPrecExpCloverLinOp::evenEvenInvLinOp(LatticeFermion& chi, const LatticeFermion& psi, 
					   enum PlusMinus isign) const
  {
    START_CODE();

    swatch.reset(); swatch.start();

    clov.applyInv(chi, psi, isign, 0);
    chi /= (Real(Nd) + param.Mass);

    swatch.stop();
    clov_apply_time += swatch.getTimeInSeconds();
    
    END_CODE();
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
  EvenOddPrecExpCloverLinOp::evenOddLinOp(LatticeFermion& chi, 
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
  EvenOddPrecExpCloverLinOp::oddEvenLinOp(LatticeFermion& chi, 
				       const LatticeFermion& psi, 
				       enum PlusMinus isign) const
  {
    START_CODE();

    Real mhalf = -0.5;

    D.apply(chi, psi, isign, 1);
    chi[rb[1]] *= mhalf;
  
    END_CODE();
  }


  //! Apply even-odd preconditioned ExpClover fermion linear operator
  /*!
   * \param chi 	  Pseudofermion field     	       (Write)
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void EvenOddPrecExpCloverLinOp::operator()(LatticeFermion & chi, 
					  const LatticeFermion& psi, 
					  enum PlusMinus isign) const
  {
    START_CODE();

    Real mquarter = -0.25;
  
    //  tmp1_o  =  D_oe   A^(-1)_ee  D_eo  psi_o
    D.apply(tmp1, psi, isign, 0);

    swatch.reset(); swatch.start();
    clov.applyInv(tmp2, tmp1, isign, 0);   
    tmp2 /= (Real(Nd) + param.Mass);

    swatch.stop();
    clov_apply_time += swatch.getTimeInSeconds();

    D.apply(tmp1, tmp2, isign, 1);

    //  chi_o  =  A_oo  psi_o  -  tmp1_o
    swatch.reset(); swatch.start();
    clov.apply(chi, psi, isign, 1);
    chi *= (Real(Nd) + param.Mass);

    swatch.stop();
    clov_apply_time += swatch.getTimeInSeconds();

    chi[rb[1]] += mquarter*tmp1;

    // Twisted Term?
    if( param.twisted_m_usedP ){ 
      // tmp1 = i mu gamma_5 tmp1
      tmp1[rb[1]] = (Gamma(15) * timesI(psi));
      
      if( isign == PLUS ) {
	chi[rb[1]] += param.twisted_m * tmp1;
      }
      else {
	chi[rb[1]] -= param.twisted_m * tmp1;
      }
    }

    END_CODE();
  }


  //! Apply the even-even block onto a source std::vector
  void 
  EvenOddPrecExpCloverLinOp::derivEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
					     const LatticeFermion& chi, const LatticeFermion& psi, 
					     enum PlusMinus isign) const
  {
    START_CODE();
    
    swatch.reset(); swatch.start();
    clov.deriv(ds_u, chi, psi, isign, 0);
    for (int mu = 0; mu < Nd; mu++)
    {
      ds_u[mu] *= (Real(Nd) + param.Mass);
    }

    swatch.stop();
    clov_deriv_time  += swatch.getTimeInSeconds();

    END_CODE();
  }

  //! Apply the even-even block onto a source std::vector
  void 
  EvenOddPrecExpCloverLinOp::derivEvenEvenLinOpMP(multi1d<LatticeColorMatrix>& ds_u, 
					       const multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
					       enum PlusMinus isign) const
  {
    START_CODE();
    
    swatch.reset(); swatch.start();
    clov.derivMultipole(ds_u, chi, psi, isign, 0);

    for (int mu = 0; mu < Nd; mu++)
    {
      ds_u[mu] *= (Real(Nd) + param.Mass);
    }

    swatch.stop();
    clov_deriv_time  += swatch.getTimeInSeconds();

    END_CODE();
  }

  //! Apply the even-even block onto a source std::vector
  void 
  EvenOddPrecExpCloverLinOp::derivLogDetEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u,
						   enum PlusMinus isign) const
  {
    START_CODE();

    //invclov.derivTrLn(ds_u, isign, 0);
    // Testing Odd Odd Term - get nothing from even even term
    clov.derivTrLn(ds_u, isign, 0);
    for (int mu = 0; mu < Nd; mu++)
    {
      ds_u[mu] *= (Real(Nd) + param.Mass);
    }
        

    END_CODE();
  }

  //! Apply the the even-odd block onto a source std::vector
  void 
  EvenOddPrecExpCloverLinOp::derivEvenOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
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
 
  //! Apply the the odd-even block onto a source std::vector
  void 
  EvenOddPrecExpCloverLinOp::derivOddEvenLinOp(multi1d<LatticeColorMatrix>& ds_u, 
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
  //! Apply the the odd-odd block onto a source std::vector
  void 
  EvenOddPrecExpCloverLinOp::derivOddOddLinOp(multi1d<LatticeColorMatrix>& ds_u, 
					   const LatticeFermion& chi, const LatticeFermion& psi, 
					   enum PlusMinus isign) const
  {   
    START_CODE();

    swatch.reset(); swatch.start();
    clov.deriv(ds_u, chi, psi, isign, 1);
    for (int mu = 0; mu < Nd; mu++)
    {
      ds_u[mu] *= (Real(Nd) + param.Mass);
    }

    swatch.stop();
    clov_deriv_time += swatch.getTimeInSeconds();
    
    END_CODE();
  }

  void 
  EvenOddPrecExpCloverLinOp::derivOddOddLinOpMP(multi1d<LatticeColorMatrix>& ds_u, 
					       const multi1d<LatticeFermion>& chi, const multi1d<LatticeFermion>& psi, 
					       enum PlusMinus isign) const
  {
    START_CODE();
    
    swatch.reset(); swatch.start();
    clov.derivMultipole(ds_u, chi, psi, isign, 1);
    for (int mu = 0; mu < Nd; mu++)
    {
     ds_u[mu] *= (Real(Nd) + param.Mass);
    }

    swatch.stop();
    clov_deriv_time  += swatch.getTimeInSeconds();

    END_CODE();
  }

  //! Return flops performed by the operator()
  unsigned long EvenOddPrecExpCloverLinOp::nFlops() const
  {
    unsigned long cbsite_flops = 2*D.nFlops()+2*clov.nFlops()+4*Nc*Ns;
    if(  param.twisted_m_usedP ) { 
      cbsite_flops += 4*Nc*Ns; // a + mu*b : a = chi, b = g_5 I psi
    }
    return cbsite_flops*(Layout::sitesOnNode()/2);
  }

  //! Get the log det of the even even part
  // BUt for now, return zero for testing.
  Double EvenOddPrecExpCloverLinOp::logDetEvenEvenLinOp(void) const  {

    return invclov.cholesDet(0);
    //return clov.cholesDet(0);

  }
} // End Namespace Chroma
