// $Id: eoprec_clover_dumb_linop_w.cc,v 3.1 2009-04-17 02:05:33 bjoo Exp $
/*! \file
 *  \brief Even-odd preconditioned clover linear operator
 */

#include "actions/ferm/linop/eoprec_clover_dumb_linop_w.h"



namespace Chroma 
{
#if 0
  //using namespace QDP::Hints;

  //! Creation routine with Anisotropy
  /*!
   * \param u_ 	    gauge field     	       (Read)
   * \param param_  fermion kappa   	       (Read)
   */
  void EvenOddPrecDumbCloverFLinOp::create(Handle< FermState<T,P,Q> > fs, 
					   const CloverFermActParams& param_)
  {
    START_CODE();

    // QDPIO::cout << __PRETTY_FUNCTION__ << ": enter" << endl;

    param = param_;
    clov.create(fs, param);
    invclov.create(fs,param,clov);  // make a copy
    invclov.choles(0);  // invert the cb=0 part
    D.create(fs, param.anisoParam);

    // QDPIO::cout << __PRETTY_FUNCTION__ << ": exit" << endl;
    END_CODE();
  }



  //! Apply even-odd preconditioned Clover fermion linear operator
  /*!
   * \param chi 	  Pseudofermion field     	       (Write)
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void EvenOddPrecDumbCloverFLinOp::operator()(T& chi, 
					  const T& psi, 
					  enum PlusMinus isign) const
  {
    START_CODE();

    T tmp1; 
    T tmp2; 
    Real mquarter = -0.25;
  
    //  tmp1_o  =  D_oe   A^(-1)_ee  D_eo  psi_o
    D.apply(tmp1, psi, isign, 0);

    invclov.apply(tmp2, tmp1, isign, 0);

    D.apply(tmp1, tmp2, isign, 1);

    //  chi_o  =  A_oo  psi_o  -  tmp1_o
    clov.apply(chi, psi, isign, 1);

    chi[rb[1]] += mquarter*tmp1;

    // Twisted Term?
    if( param.twisted_m_usedP ){ 
      // tmp1 = i mu gamma_5 tmp1
      tmp1[rb[1]] = (GammaConst<Ns,Ns*Ns-1>() * timesI(psi));
      
      if( isign == PLUS ) {
	chi[rb[1]] += param.twisted_m * tmp1;
      }
      else {
	chi[rb[1]] -= param.twisted_m * tmp1;
      }
    }

    END_CODE();
  }



  //! Return flops performed by the operator()
  unsigned long EvenOddPrecDumbCloverFLinOp::nFlops() const
  {
    unsigned long cbsite_flops = 2*D.nFlops()+2*clov.nFlops()+4*Nc*Ns;
    if(  param.twisted_m_usedP ) { 
      cbsite_flops += 4*Nc*Ns; // a + mu*b : a = chi, b = g_5 I psi
    }
    return cbsite_flops*(Layout::sitesOnNode()/2);
  }

  const Subset& EvenOddPrecDumbCloverFLinOp::subset() const {
    return rb[1];
  }


  //! Creation routine with Anisotropy
  /*!
   * \param u_ 	    gauge field     	       (Read)
   * \param param_  fermion kappa   	       (Read)
   */
  void EvenOddPrecDumbCloverDLinOp::create(Handle< FermState<T,P,Q> > fs, 
					   const CloverFermActParams& param_)
  {
    START_CODE();

    // QDPIO::cout << __PRETTY_FUNCTION__ << ": enter" << endl;

    param = param_;
    clov.create(fs, param);
    invclov.create(fs,param,clov);  // make a copy
    invclov.choles(0);  // invert the cb=0 part
    D.create(fs, param.anisoParam);

    // QDPIO::cout << __PRETTY_FUNCTION__ << ": exit" << endl;
    END_CODE();
  }



  //! Apply even-odd preconditioned Clover fermion linear operator
  /*!
   * \param chi 	  Pseudofermion field     	       (Write)
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void EvenOddPrecDumbCloverDLinOp::operator()(T& chi, 
					  const T& psi, 
					  enum PlusMinus isign) const
  {
    START_CODE();

    T tmp1; moveToFastMemoryHint(tmp1);
    T tmp2; moveToFastMemoryHint(tmp2);
    Real mquarter = -0.25;
  
    //  tmp1_o  =  D_oe   A^(-1)_ee  D_eo  psi_o
    D.apply(tmp1, psi, isign, 0);

    invclov.apply(tmp2, tmp1, isign, 0);

    D.apply(tmp1, tmp2, isign, 1);

    //  chi_o  =  A_oo  psi_o  -  tmp1_o
    clov.apply(chi, psi, isign, 1);

    chi[rb[1]] += mquarter*tmp1;

    // Twisted Term?
    if( param.twisted_m_usedP ){ 
      // tmp1 = i mu gamma_5 tmp1
      tmp1[rb[1]] = (GammaConst<Ns,Ns*Ns-1>() * timesI(psi));
      
      if( isign == PLUS ) {
	chi[rb[1]] += param.twisted_m * tmp1;
      }
      else {
	chi[rb[1]] -= param.twisted_m * tmp1;
      }
    }

    END_CODE();
  }



  //! Return flops performed by the operator()
  unsigned long EvenOddPrecDumbCloverDLinOp::nFlops() const
  {
    unsigned long cbsite_flops = 2*D.nFlops()+2*clov.nFlops()+4*Nc*Ns;
    if(  param.twisted_m_usedP ) { 
      cbsite_flops += 4*Nc*Ns; // a + mu*b : a = chi, b = g_5 I psi
    }
    return cbsite_flops*(Layout::sitesOnNode()/2);
  }

  const Subset& EvenOddPrecDumbCloverDLinOp::subset() const {
    return rb[1];
  }

#endif

} // End Namespace Chroma
