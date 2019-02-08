/*! \file
 *  \brief Even-odd preconditioned clover linear operator
 */

#include "chromabase.h"
#include "actions/ferm/linop/unprec_clover_plus_igmuA2_linop_w.h"

using namespace QDP::Hints;

namespace Chroma 
{ 
  //! Creation routine with Anisotropy
  /*!
   * \param fs 	    gauge field     	       (Read)
   * \param param_  parameters   	       (Read)
   */
  void UnprecCloverPlusIG5MuA2Linop::create(Handle< FermState<T,P,Q> > fs,
				 const CloverFermActParams& param_)
  {
    //   QDPIO::cout << __PRETTY_FUNCTION__ << ": enter" << std::endl;

    param = param_;

    A.create(fs, param);
    D.create(fs, param.anisoParam);

    // QDPIO::cout << __PRETTY_FUNCTION__ << ": exit" << std::endl;
  }



  //! Apply unpreconditioned Clover fermion linear operator
  /*!
   * The operator acts on the entire lattice
   *
   * \param chi 	  Pseudofermion field          (Write)
   * \param psi 	  Pseudofermion field          (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void UnprecCloverPlusIG5MuA2Linop::operator()(LatticeFermion & chi,
				     const LatticeFermion& psi, 
				     enum PlusMinus isign) const
  {
    LatticeFermion tmp,tmp2; moveToFastMemoryHint(tmp);
    Real mhalf = -0.5;

    //  chi   =  A . psi - 0.5 * D' . psi  */
    A(chi, psi, isign); 
    D(tmp, psi, isign);
    chi += mhalf * tmp;

    if( param.twisted_m_usedP ) {
    	// tmp2 =  A^2 psi
    	A.apply(tmp,psi,isign,1);
    	A.apply(tmp2,tmp,isign,1);

    	// tmp_even = zero
    	tmp2[rb[0]] = zero;

    	// tmp_odd = gamma_5 tmp2 = gamma_5 A^2 psi
    	tmp[rb[1]] = Gamma(15)*tmp2;

    	if ( isign == PLUS ) {
    		chi[rb[1]] -= param.twisted_m * timesI(tmp);
    	}
    	else {
    		chi[rb[1]] += param.twisted_m * timesI(tmp);
    	}
    }

    getFermBC().modifyF(chi);
  }


  void 
  UnprecCloverPlusIG5MuA2Linop::deriv(multi1d<LatticeColorMatrix>& ds_u,
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
    
    if( param.twisted_m_usedP ) {
    	//  i X^\dag d/dt [ g_5 A^2 ] Y
    	// = X^dag g_5 [ \dot{A^\dag}A +  A^\dag\dot{A} ] Y
    	// = X^dag g_5 \dot{A^dag} AY + X^dag g_5 A^dag \dot{A} Y
    	// = W^dag \dot{A^dag} Z  + W' ^dag \dot{A} Y
    	//
    	// write code here...

    	// first term: X^dag g_5 \dot{A^dag} AY = W^dag \dot{A^dag} Z
    	//  g_5 X = W   A Y = Z

    	LatticeFermion tmp_W = zero;
    	LatticeFermion tmp_Z = zero;
    	tmp_W[rb[1]] = Gamma(15)*chi;
    	A.apply(tmp_Z,psi,isign,1);


    	ds_tmp = zero;
    	A.deriv(ds_tmp,tmp_W, tmp_Z, isign, 1);
    	for(int mu=0; mu < Nd; ++mu) {
    		if( isign == PLUS ) {
    			ds_u[mu] -= param.twisted_m * timesI(ds_tmp[mu]);
    		}
    		else {
    			ds_u[mu] += param.twisted_m * timesI(ds_tmp[mu]);
    		}
    	}

    	// Second term: X^dag g_5 A^dag \dot{A} Y = Z^\dag \dot{A} Y
    	// Z^\dag = X^dag g_5 A^\dag => A g_5 X = A W = Z

    	A.apply(tmp_Z, tmp_W, isign,1);
    	ds_tmp = zero;
    	A.deriv(ds_tmp, tmp_Z, chi, isign, 1);
    	for(int mu=0; mu < Nd; ++mu) {
    		if( isign == PLUS ) {
    			ds_u[mu] -= param.twisted_m * timesI(ds_tmp[mu]);
    		}
    		else {
    			ds_u[mu] += param.twisted_m * timesI(ds_tmp[mu]);
    		}
    	}
    } // if twisted_m_usedP

    getFermBC().zero(ds_u);
  }


  //! Return flops performed by the operator()
  unsigned long UnprecCloverPlusIG5MuA2Linop::nFlops() const
  {
    unsigned long site_flops = D.nFlops()+A.nFlops()+4*Nc*Ns;
    // A^2 = 2A_flops  -/+ imu_g5 = 2 real flops
    //                              x complex component x spin x color
    unsigned long twist_site_flops = 2*A.nFlops() + 4*Nc*Ns;
    return site_flops*twist_site_flops*Layout::sitesOnNode();
  }

} // End Namespace Chroma
