/*! \file
 *  \brief Even-odd preconditioned clover linear operator
 */

#include "actions/ferm/linop/eoprec_clover_linop_w.h"



namespace Chroma 
{ 

   using namespace QDP::Hints;

 //! Creation routine with Anisotropy
  /*!
   * \param u_ 	    gauge field     	       (Read)
   * \param param_  fermion kappa   	       (Read)
   */
  template<typename T_, typename U_, typename DslashT, typename CloverT>
  void EvenOddPrecCloverLinOpT<T_,U_,DslashT,CloverT>::create(Handle< FermState<T,P,Q> > fs, 
				      const CloverFermActParams& param_)
  {
    START_CODE();
    // QDPIO::cout << __PRETTY_FUNCTION__ << ": enter" << std::endl;

    param = param_;

    clov.create(fs, param);
 
    invclov.create(fs,param,clov);  // make a copy
    invclov.choles(0);  // invert the cb=0 part

    D.create(fs, param.anisoParam);

    clov_deriv_time = 0;
    clov_apply_time = 0;

    moveToFastMemoryHint(tmp1);
    moveToFastMemoryHint(tmp2);
     
    // QDPIO::cout << __PRETTY_FUNCTION__ << ": exit" << std::endl;
    END_CODE();
  }

  //! Apply the the odd-odd block onto a source std::vector
  template<typename T_, typename U_, typename DslashT, typename CloverT>
  void 
  EvenOddPrecCloverLinOpT<T_,U_,DslashT,CloverT>::oddOddLinOp(T& chi, const T& psi, 
				      enum PlusMinus isign) const
  {
    START_CODE();

    swatch.reset(); swatch.start();
    clov.apply(chi, psi, isign, 1);
    swatch.stop();
    clov_apply_time += swatch.getTimeInSeconds();

    END_CODE();
  }


  //! Apply the the even-even block onto a source std::vector
  template<typename T_, typename U_, typename DslashT, typename CloverT>
  void 
  EvenOddPrecCloverLinOpT<T_,U_,DslashT,CloverT>::evenEvenLinOp(T& chi, const T& psi, 
					enum PlusMinus isign) const
  {
    START_CODE();

    // Nuke for testing
    swatch.reset(); swatch.start();
    clov.apply(chi, psi, isign, 0);
    swatch.stop();
    clov_apply_time += swatch.getTimeInSeconds();
    
    END_CODE();
  }

  //! Apply the inverse of the even-even block onto a source std::vector
  template<typename T_, typename U_, typename DslashT, typename CloverT>
  void 
  EvenOddPrecCloverLinOpT<T_,U_,DslashT,CloverT>::evenEvenInvLinOp(T& chi, const T& psi, 
					   enum PlusMinus isign) const
  {
    START_CODE();

    swatch.reset(); swatch.start();
    invclov.apply(chi, psi, isign, 0);
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
  template<typename T_, typename U_, typename DslashT, typename CloverT>
  void 
  EvenOddPrecCloverLinOpT<T_,U_,DslashT,CloverT>::evenOddLinOp(T& chi, 
				       const T& psi, 
				       enum PlusMinus isign) const
  {
    START_CODE();

    RealT mhalf = -0.5;

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
  template<typename T_, typename U_, typename DslashT, typename CloverT>
  void 
  EvenOddPrecCloverLinOpT<T_,U_,DslashT,CloverT>::oddEvenLinOp(T& chi, 
				       const T& psi, 
				       enum PlusMinus isign) const
  {
    START_CODE();

    RealT mhalf = -0.5;

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
  template<typename T_, typename U_, typename DslashT, typename CloverT>
  void EvenOddPrecCloverLinOpT<T_,U_,DslashT,CloverT>::operator()(T & chi, 
					  const T& psi, 
					  enum PlusMinus isign) const
  {
    START_CODE();

    RealT mquarter = -0.25;
  
    //  tmp1_o  =  D_oe   A^(-1)_ee  D_eo  psi_o
    D.apply(tmp1, psi, isign, 0);

    swatch.reset(); swatch.start();
    invclov.apply(tmp2, tmp1, isign, 0);
    swatch.stop();
    clov_apply_time += swatch.getTimeInSeconds();

    D.apply(tmp1, tmp2, isign, 1);

    //  chi_o  =  A_oo  psi_o  -  tmp1_o
    swatch.reset(); swatch.start();
    clov.apply(chi, psi, isign, 1);
    swatch.stop();
    clov_apply_time += swatch.getTimeInSeconds();

    chi[rb[1]] += mquarter*tmp1;

    // Twisted Term?
    if( param.twisted_m_usedP ){ 
      RealT twisted_m = toDouble(param.twisted_m);
      // tmp1 = i mu gamma_5 tmp1
      tmp1[rb[1]] = (Gamma(15) * timesI(psi));
      
      if( isign == PLUS ) {
	chi[rb[1]] += twisted_m * tmp1;
      }
      else {
	chi[rb[1]] -= twisted_m * tmp1;
      }
    }

    END_CODE();
  }


  //! Apply the even-even block onto a source std::vector
  template<typename T_, typename U_, typename DslashT, typename CloverT>
  void 
  EvenOddPrecCloverLinOpT<T_,U_,DslashT,CloverT>::derivEvenEvenLinOp(P& ds_u, 
					     const T& chi, const T& psi, 
					     enum PlusMinus isign) const
  {
    START_CODE();
    
    swatch.reset(); swatch.start();
    clov.deriv(ds_u, chi, psi, isign, 0);
    swatch.stop();
    clov_deriv_time  += swatch.getTimeInSeconds();

    END_CODE();
  }

  //! Apply the even-even block onto a source std::vector
  template<typename T_, typename U_, typename DslashT, typename CloverT>
  void 
  EvenOddPrecCloverLinOpT<T_,U_,DslashT,CloverT>::derivEvenEvenLinOpMP(P& ds_u, 
					       const multi1d<T>& chi, const multi1d<T>& psi, 
					       enum PlusMinus isign) const
  {
    START_CODE();
    
    swatch.reset(); swatch.start();
    clov.derivMultipole(ds_u, chi, psi, isign, 0);
    swatch.stop();
    clov_deriv_time  += swatch.getTimeInSeconds();

    END_CODE();
  }

  //! Apply the even-even block onto a source std::vector
  template<typename T_, typename U_, typename DslashT, typename CloverT>
  void 
  EvenOddPrecCloverLinOpT<T_,U_,DslashT,CloverT>::derivLogDetEvenEvenLinOp(P& ds_u,
						   enum PlusMinus isign) const
  {
    START_CODE();

    // Testing Odd Odd Term - get nothing from even even term
    invclov.derivTrLn(ds_u, isign, 0);
    
    END_CODE();
  }

  //! Apply the the even-odd block onto a source std::vector
  template<typename T_, typename U_, typename DslashT, typename CloverT>
  void 
  EvenOddPrecCloverLinOpT<T_,U_,DslashT,CloverT>::derivEvenOddLinOp(P& ds_u, 
					    const T& chi, const T& psi, 
					    enum PlusMinus isign) const
  {
    START_CODE();
    ds_u.resize(Nd);
    D.deriv(ds_u, chi, psi, isign, 0);
    for(int mu=0; mu < Nd; mu++) { 
      ds_u[mu]  *= RealT(-0.5);
    }
    END_CODE();
  }
 
  //! Apply the the odd-even block onto a source std::vector
  template<typename T_, typename U_, typename DslashT, typename CloverT>
  void 
  EvenOddPrecCloverLinOpT<T_,U_,DslashT,CloverT>::derivOddEvenLinOp(P& ds_u, 
					    const T& chi, const T& psi, 
					    enum PlusMinus isign) const
  {
    START_CODE();
    ds_u.resize(Nd);

    D.deriv(ds_u, chi, psi, isign, 1);
    for(int mu=0; mu < Nd; mu++) { 
     ds_u[mu]  *= RealT(-0.5);
    }
    END_CODE();
  }

  // Inherit this
  //! Apply the the odd-odd block onto a source std::vector
  template<typename T_, typename U_, typename DslashT, typename CloverT>
  void 
  EvenOddPrecCloverLinOpT<T_,U_,DslashT,CloverT>::derivOddOddLinOp(P& ds_u, 
					   const T& chi, const T& psi, 
					   enum PlusMinus isign) const
  {   
    START_CODE();

    swatch.reset(); swatch.start();
    clov.deriv(ds_u, chi, psi, isign, 1);
    swatch.stop();
    clov_deriv_time += swatch.getTimeInSeconds();
    
    END_CODE();
  }

  template<typename T_, typename U_, typename DslashT, typename CloverT>
  void 
  EvenOddPrecCloverLinOpT<T_,U_,DslashT,CloverT>::derivOddOddLinOpMP(P& ds_u, 
					       const multi1d<T>& chi, const multi1d<T>& psi, 
					       enum PlusMinus isign) const
  {
    START_CODE();
    
    swatch.reset(); swatch.start();
    clov.derivMultipole(ds_u, chi, psi, isign, 1);
    swatch.stop();
    clov_deriv_time  += swatch.getTimeInSeconds();

    END_CODE();
  }

  //! Return flops performed by the operator()
  template<typename T_, typename U_, typename DslashT, typename CloverT>
  unsigned long EvenOddPrecCloverLinOpT<T_,U_,DslashT,CloverT>::nFlops() const
  {
    unsigned long cbsite_flops = 2*D.nFlops()+2*clov.nFlops()+4*Nc*Ns;
    if(  param.twisted_m_usedP ) { 
      cbsite_flops += 4*Nc*Ns; // a + mu*b : a = chi, b = g_5 I psi
    }
    return cbsite_flops*(Layout::sitesOnNode()/2);
  }

  //! Get the log det of the even even part
  // BUt for now, return zero for testing.
  template<typename T_, typename U_, typename DslashT, typename CloverT>
  Double EvenOddPrecCloverLinOpT<T_,U_,DslashT,CloverT>::logDetEvenEvenLinOp(void) const  {
    return invclov.cholesDet(0);
  }


  template class EvenOddPrecCloverLinOpT<LatticeFermion, LatticeColorMatrix, WilsonDslash, CloverTerm>;
#if BASE_PRECISION == 64
  template class EvenOddPrecCloverLinOpT<LatticeFermionF, LatticeColorMatrixF, WilsonDslashF, CloverTermF>;
#endif

} // End Namespace Chroma
