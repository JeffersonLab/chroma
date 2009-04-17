// -*- C++ -*-
// $Id: eoprec_clover_dumb_linop_w.h,v 3.1 2009-04-17 02:05:33 bjoo Exp $
/*! \file
 *  \brief Even-odd preconditioned Clover fermion linear operator
 */

#ifndef __prec_clover_dumb_linop_w_h__
#define __prec_clover_dumb_linop_w_h__

#include "state.h"
#include "fermbc.h"
#include "eoprec_logdet_linop.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/clover_term_w.h"


namespace Chroma 
{ 

#if 1 
  //! Even-odd preconditioned Clover-Dirac operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   *
   * The kernel for Clover fermions is
   *
   *      M  =  A + (d+M) - (1/2) D'

   * This is a dumb version with only a constructor and an 
   * apply method. It is also fixed to single precision
   */
  class EvenOddPrecDumbCloverFLinOp : public LinearOperator<LatticeFermionF>
  {
  public:
    typedef LatticeFermionF T;
    typedef LatticeColorMatrixF U;
    typedef multi1d<U> P;
    typedef multi1d<U> Q;

    //! Full constructor
    EvenOddPrecDumbCloverFLinOp(Handle< FermState<LatticeFermionF,P,Q> > fs,
				const CloverFermActParams& param_)
      {create(fs,param_);}

    //! Destructor is automatic
    ~EvenOddPrecDumbCloverFLinOp() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<LatticeFermionF,P,Q>& getFermBC() const {return D.getFermBC();}

    //! Creation routine
    void create(Handle< FermState<LatticeFermionF,P,Q> > fs,
		const CloverFermActParams& param_) {

      param = param_;
      clov.create(fs, param);
      invclov.create(fs,param,clov);  // make a copy
      invclov.choles(0);  // invert the cb=0 part
      D.create(fs, param.anisoParam);
      
    }


    // Override inherited one with a few more funkies
    void operator()(LatticeFermionF& chi, const LatticeFermionF& psi, 
		    enum PlusMinus isign) const 
    {
      START_CODE();
      
      LatticeFermionF tmp1; 
      LatticeFermionF tmp2; 
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
    unsigned long nFlops() const 
    {
      unsigned long cbsite_flops = 2*D.nFlops()+2*clov.nFlops()+4*Nc*Ns;
      if(  param.twisted_m_usedP ) { 
	cbsite_flops += 4*Nc*Ns; // a + mu*b : a = chi, b = g_5 I psi
      }
      return cbsite_flops*(Layout::sitesOnNode()/2);
    }

    const Subset& subset() const 
    {
      return rb[1];
    }

  private:
    CloverFermActParams param;
    WilsonDslashF D;
    CloverTermF   clov;
    CloverTermF   invclov;  // uggh, only needed for evenEvenLinOp
  };

#endif
#if 1
  //! Even-odd preconditioned Clover-Dirac operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   *
   * The kernel for Clover fermions is
   *
   *      M  =  A + (d+M) - (1/2) D'

   * This is a dumb version with only a constructor and an 
   * apply method. It is also fixed to double precision
   */
  class EvenOddPrecDumbCloverDLinOp : public LinearOperator<LatticeFermionD>
  {
  public:
    typedef LatticeFermionD T;
    typedef LatticeColorMatrixD U;
    typedef multi1d<U> P;
    typedef multi1d<U> Q;

    //! Full constructor
    EvenOddPrecDumbCloverDLinOp(Handle< FermState<T,P,Q> > fs,
				const CloverFermActParams& param_)
      {create(fs,param_);}

    //! Destructor is automatic
    ~EvenOddPrecDumbCloverDLinOp() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return D.getFermBC();}

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > fs,
		const CloverFermActParams& param_) {
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

 
    // Override inherited one with a few more funkies
    void operator()(T& chi, const T& psi, 
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
    unsigned long nFlops() const
    {
      unsigned long cbsite_flops = 2*D.nFlops()+2*clov.nFlops()+4*Nc*Ns;
      if(  param.twisted_m_usedP ) { 
	cbsite_flops += 4*Nc*Ns; // a + mu*b : a = chi, b = g_5 I psi
      }
      return cbsite_flops*(Layout::sitesOnNode()/2);
    }

    const Subset& subset() const 
    {
      return rb[1];
    }

  private:
    CloverFermActParams param;
    WilsonDslashD D;
    CloverTermD   clov;
    CloverTermD   invclov;  // uggh, only needed for evenEvenLinOp
  };

#endif

} // End Namespace Chroma


#endif
