/*! \file
 *  \brief Symmetric even-odd preconditioned clover linear operator
 */

#include "actions/ferm/linop/seoprec_clover_linop_w.h"

namespace Chroma 
{ 
  using namespace QDP::Hints;

  //! Creation routine with Anisotropy
  /*!
   * \param u_ 	    gauge field     	       (Read)
   * \param param_  fermion kappa   	       (Read)
   */
  void SymEvenOddPrecCloverLinOp::create(Handle< FermState<T,P,Q> > fs, 
					 const CloverFermActParams& param_)
  {
    START_CODE();
    // QDPIO::cout << __PRETTY_FUNCTION__ << ": enter" << std::endl;

    param = param_;
    QDPIO::cout << "Using Twisted Mass: " << param.twisted_m_usedP << std::endl;
    if( param.twisted_m_usedP) {
    	QDPIO::cout << "Twisted Mass is " << param.twisted_m << std::endl;
    }

    clov.create(fs, param);
 
    invclov.create(fs,param,clov);  // make a copy
    invclov.choles(0);  // invert the cb=0 part
    invclov.choles(1);  // invert the cb=1 part

    D.create(fs, param.anisoParam);

    clov_deriv_time = 0;
    clov_apply_time = 0;

    // QDPIO::cout << __PRETTY_FUNCTION__ << ": exit" << std::endl;
    END_CODE();
  }

  //! Apply the the odd-odd block onto a source std::vector
  void 
  SymEvenOddPrecCloverLinOp::unprecOddOddLinOp(LatticeFermion& chi, const LatticeFermion& psi,
					 enum PlusMinus isign) const
  {
    START_CODE();


    clov.apply(chi, psi, isign, 1);

    END_CODE();
  }


  //! Apply the inverse of the odd-odd block onto a source std::vector
  void 
  SymEvenOddPrecCloverLinOp::unprecOddOddInvLinOp(LatticeFermion& chi, const LatticeFermion& psi,
					    enum PlusMinus isign) const
  {
    START_CODE();

    invclov.apply(chi, psi, isign, 1);
    END_CODE();
  }
  

  //! Apply the the even-even block onto a source std::vector
  void 
  SymEvenOddPrecCloverLinOp::unprecEvenEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi,
					   enum PlusMinus isign) const
  {
    START_CODE();
    clov.apply(chi, psi, isign, 0);
    END_CODE();
  }

  //! Apply the inverse of the even-even block onto a source std::vector
  void 
  SymEvenOddPrecCloverLinOp::unprecEvenEvenInvLinOp(LatticeFermion& chi, const LatticeFermion& psi,
					      enum PlusMinus isign) const
  {
    START_CODE();

    invclov.apply(chi, psi, isign, 0);
    
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
  SymEvenOddPrecCloverLinOp::unprecEvenOddLinOp(LatticeFermion& chi,
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
  SymEvenOddPrecCloverLinOp::unprecOddEvenLinOp(LatticeFermion& chi,
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
  void SymEvenOddPrecCloverLinOp::operator()(T& chi,
					     const T& psi,
					     enum PlusMinus isign) const
  {
    START_CODE();

    T tmp1; moveToFastMemoryHint(tmp1);
    T tmp2; moveToFastMemoryHint(tmp2);
    Real mquarter = -0.25;
 

    if( isign == PLUS) {  
       //  tmp2_o  =  A^(-1)_oo  D_oe  A^(-1)_ee  D_eo  psi_o
       D.apply(tmp1, psi, isign, 0);
       invclov.apply(tmp2, tmp1, isign, 0);

       D.apply(tmp1, tmp2, isign, 1);
       invclov.apply(tmp2, tmp1, isign, 1);
     }
     else { 
       invclov.apply(tmp1, psi, isign, 1);
       D.apply(tmp2, tmp1, isign, 0);
       invclov.apply(tmp1, tmp2, isign, 0);
       D.apply(tmp2, tmp1, isign, 1);
     }

    //  chi_o  =  psi_o  -  tmp2_o
    chi[rb[1]] = psi  +  mquarter*tmp2;
#if 1
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
#endif
    getFermBC().modifyF(chi);
    END_CODE();
  }


  //! Apply the even-even block onto a source std::vector
  void 
  SymEvenOddPrecCloverLinOp::derivUnprecEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u,
						const LatticeFermion& chi, const LatticeFermion& psi, 
						enum PlusMinus isign) const
  {
    START_CODE();
    

    clov.deriv(ds_u, chi, psi, isign, 0);

    END_CODE();
  }

  // Inherit this
  //! Apply the the odd-odd block onto a source std::vector
  void
  SymEvenOddPrecCloverLinOp::derivUnprecOddOddLinOp(multi1d<LatticeColorMatrix>& ds_u,
					      const LatticeFermion& chi, const LatticeFermion& psi,
					      enum PlusMinus isign) const
  {
    START_CODE();

    clov.deriv(ds_u, chi, psi, isign, 1);

    END_CODE();
  }

  //! Apply the even-even block onto a source std::vector
  void 
  SymEvenOddPrecCloverLinOp::derivLogDetEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u,
						      enum PlusMinus isign) const
  {
    START_CODE();

    invclov.derivTrLn(ds_u, isign, 0);
    
    END_CODE();
  }

  //! Apply the even-even block onto a source std::vector
  void 
  SymEvenOddPrecCloverLinOp::derivLogDetOddOddLinOp(multi1d<LatticeColorMatrix>& ds_u,
						    enum PlusMinus isign) const
  {
    START_CODE();

    invclov.derivTrLn(ds_u, isign, 1);
    
    END_CODE();
  }

  //! Apply the the even-odd block onto a source std::vector
  void 
  SymEvenOddPrecCloverLinOp::derivUnprecEvenOddLinOp(multi1d<LatticeColorMatrix>& ds_u,
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
  SymEvenOddPrecCloverLinOp::derivUnprecOddEvenLinOp(multi1d<LatticeColorMatrix>& ds_u,
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

	//! Deriv
	void SymEvenOddPrecCloverLinOp::deriv(P& ds_u, const T& chi, const T& psi,
			enum PlusMinus isign) const
	{
		T M_eo_psi; moveToFastMemoryHint(M_eo_psi);
		T M_oe_dag_chi; moveToFastMemoryHint(M_oe_dag_chi);
		T tmp;
		T Ainv_left;
		T Ainv_right;

		enum PlusMinus msign = ( isign == PLUS ) ? MINUS : PLUS;

		if( isign == PLUS ) {
			D.apply(tmp, psi, PLUS, 0);
			invclov.apply(M_eo_psi, tmp, PLUS, 0);

			invclov.apply(tmp, chi, MINUS,1);
			D.apply(M_oe_dag_chi,tmp,MINUS,0);
		}
		else {
			invclov.apply(tmp,psi,MINUS,1);
			D.apply(M_eo_psi,tmp,MINUS,0);

			D.apply(tmp,chi,PLUS,0);
			invclov.apply(M_oe_dag_chi,tmp,PLUS,0);
		}
		M_eo_psi[rb[0]] *= Real(-0.5);
		M_oe_dag_chi[rb[0]] *= Real(-0.5);

		ds_u.resize(Nd);
		ds_u = zero;
		P ds_tmp;
		ds_tmp.resize(Nd);

		//derivOddEvenLinOp(ds_tmp,chi,M_eo_psi,isign);
		if( isign == PLUS ) {

			//
		   D.apply(tmp,M_eo_psi, PLUS,1);
		   invclov.apply(Ainv_right,tmp, PLUS,1);
		   invclov.apply(Ainv_left,chi,MINUS,1);
		   clov.deriv(ds_u,Ainv_left,Ainv_right, PLUS, 1);


		   invclov.apply(tmp, chi, MINUS,1);
		   D.deriv(ds_tmp, tmp,M_eo_psi, PLUS,1);
		   ds_u -= ds_tmp;


		   D.apply(tmp,psi, PLUS,0);
		   invclov.apply(Ainv_right,tmp,PLUS,0);
		   invclov.apply(Ainv_left, M_oe_dag_chi,PLUS,0);
		   clov.deriv(ds_tmp,Ainv_left, Ainv_right,PLUS,0);
		   ds_u += ds_tmp;


		   invclov.apply(tmp, M_oe_dag_chi, MINUS,0);
		   D.deriv(ds_tmp, tmp,psi, PLUS,0);
		   ds_u -= ds_tmp;
		}
		else {
			D.apply(tmp,chi,PLUS,0);
			invclov.apply(Ainv_left,tmp,PLUS,0);
			invclov.apply(Ainv_right,M_eo_psi,MINUS,0);
			clov.deriv(ds_u,Ainv_left,Ainv_right,MINUS,0);

			invclov.apply(tmp, M_eo_psi, MINUS,0);
			D.deriv(ds_tmp, chi, tmp, MINUS,1);
			ds_u -= ds_tmp;

	    	invclov.apply(tmp, psi, MINUS,1);
	    	D.deriv(ds_tmp, M_oe_dag_chi, tmp, MINUS,0);
	    	ds_u -= ds_tmp;

	    	D.apply(tmp,M_oe_dag_chi,PLUS,1);
	    	invclov.apply(Ainv_left,tmp,PLUS,1);
	    	invclov.apply(Ainv_right,psi,MINUS,1);
	    	clov.deriv(ds_tmp, Ainv_left,Ainv_right, MINUS,1);
	    	ds_u += ds_tmp;
		}

		for(int mu=0; mu < Nd; ++mu) {
			ds_u[mu] *= Real(-0.5);
		}
		getFermBC().zero(ds_u);
	}


  //! Return flops performed by the operator()
  unsigned long SymEvenOddPrecCloverLinOp::nFlops() const
  {
    unsigned long cbsite_flops = 2*D.nFlops()+2*clov.nFlops()+4*Nc*Ns;
    if(  param.twisted_m_usedP ) { 
      cbsite_flops += 4*Nc*Ns; // a + mu*b : a = chi, b = g_5 I psi
    }
    return cbsite_flops*(Layout::sitesOnNode()/2);
  }

  //! Get the log det of the even even part
  // BUt for now, return zero for testing.
  Double SymEvenOddPrecCloverLinOp::logDetEvenEvenLinOp(void) const  {
    return invclov.cholesDet(0);
  }

  //! Get the log det of the odd odd part
  // BUt for now, return zero for testing.
  Double SymEvenOddPrecCloverLinOp::logDetOddOddLinOp(void) const  {
    return invclov.cholesDet(1);
  }

} // End Namespace Chroma
