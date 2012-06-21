// -*- C++ -*-
// $Id: clover_term_base_w.h,v 3.8 2009-04-17 02:05:32 bjoo Exp $
/*! \file
 *  \brief Clover term linear operator
 */

#ifndef __clover_term_base_w_h__
#define __clover_term_base_w_h__

#include "linearop.h"


namespace Chroma 
{ 
  //! Clover term
  /*!
   * \ingroup linop
   *
   */

  template<typename T, typename U>
	   class CloverTermBase : public DslashLinearOperator<T,
							      multi1d<U>,
							      multi1d<U> >
  {
  public:
    //! No real need for cleanup here
    virtual ~CloverTermBase() {}

    //! Subset is all here
    const Subset& subset() const {return all;}


    virtual void applySite(T& chi, const T& psi, enum PlusMinus isign, int site) const = 0;

    //! Invert
    /*!
     * Computes the inverse of the term on cb using Cholesky
     */
    virtual void choles(int cb) = 0;

    //! Invert
    /*!
     * Computes the determinant of the term
     *
     * \return logarithm of the determinant  
     */
    virtual Double cholesDet(int cb) const = 0;

    //! Take deriv of D
    /*!
     * \param chi     left vector                                 (Read)
     * \param psi     right vector                                (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     *
     * \return Computes   \f$chi^\dag * \dot(D} * psi\f$
     */
    void deriv(multi1d<U>& ds_u, 
	       const T& chi, const T& psi, 
	       enum PlusMinus isign) const;

    //! Take deriv of D
    /*!
     * \param chi     left vector on cb                           (Read)
     * \param psi     right vector on 1-cb                        (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     * \param cb      Checkerboard of chi vector                  (Read)
     *
     * \return Computes   \f$chi^\dag * \dot(D} * psi\f$
     */
    void deriv(multi1d<U>& ds_u, 
	       const T& chi, const T& psi, 
	       enum PlusMinus isign, int cb) const;

    //! Take deriv of D
    /*!
     * \param chi     left vectors                           (Read)
     * \param psi     right vectors                         (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     * \param cb      Checkerboard of chi vector                  (Read)
     *
     * \return Computes   \f$chi^\dag * \dot(D} * psi\f$
     */
    void derivMultipole(multi1d<U>& ds_u,
			const multi1d<T>& chi, const multi1d<T>& psi,
			enum PlusMinus isign) const;

    //! Take deriv of D
    /*!
     * \param chi     left vectors on cb                           (Read)
     * \param psi     right vectors on cb                        (Read)
     * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
     * \param cb      Checkerboard of chi vector                  (Read)
     *
     * \return Computes   \f$chi^\dag * \dot(D} * psi\f$
     */

    void derivMultipole(multi1d<U>& ds_u,
			const multi1d<T>& chi, const multi1d<T>& psi,
			enum PlusMinus isign, int cb) const;


    //! Take derivative of TrLn D
    void derivTrLn(multi1d<U>& ds_u, 
		   enum PlusMinus isign, int cb) const;


    void deriv_loops(const int u, const int mu, const int cb,
		     U& ds_u_mu,
		     U& ds_u_nu,
		     const U& Lambda) const;

    //! Return flops performed by the operator()
    unsigned long nFlops() const;

    //! Calculates Tr_D ( Gamma_mat L )
    virtual void triacntr(U& B, int mat, int cb) const = 0;

  protected:

    //! Get the u field
    virtual const multi1d<U>& getU() const = 0;

    //! get the clover coefficient 
    virtual Real getCloverCoeff(int mu, int nu) const = 0;

  };

  //! Return flops performed by the operator()
  template<typename T, typename U>
  unsigned long 
  CloverTermBase<T,U>::nFlops() const {return 552;}


  //! Take deriv of D
  /*!
   * \param chi     left vector                                 (Read)
   * \param psi     right vector                                (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   *
   * \return Computes   \f$\chi^\dag * \dot(D} * \psi\f$
   */
  template<typename T, typename U>
  void CloverTermBase<T,U>::deriv(multi1d<U>& ds_u, 
			     const T& chi, const T& psi, 
			     enum PlusMinus isign) const
  {
    START_CODE();

    // base deriv resizes.
    // Even even checkerboard
    deriv(ds_u, chi, psi, isign,0);
    
    // Odd Odd checkerboard
    multi1d<U> ds_tmp;
    deriv(ds_tmp, chi, psi, isign,1);
    
    ds_u += ds_tmp;
    
    END_CODE();
  }

  template<typename T, typename U>
  void CloverTermBase<T,U>::derivMultipole(multi1d<U>& ds_u, 
			     const multi1d<T>& chi, const multi1d<T>& psi, 
			     enum PlusMinus isign) const
  {
    START_CODE();

    // base deriv resizes.
    // Even even checkerboard
    derivMultipole(ds_u, chi, psi, isign,0);
    
    // Odd Odd checkerboard
    multi1d<U> ds_tmp;
    derivMultipole(ds_tmp, chi, psi, isign,1);
    
    ds_u += ds_tmp;
    
    END_CODE();
  }



  template<typename T, typename U>
  void CloverTermBase<T,U>::deriv_loops(const int mu, const int nu, const int cb,
				   U& ds_u_mu,
				   U& ds_u_nu,
				   const U& Lambda) const
  {
    START_CODE();

    const multi1d<U>& u = getU();

    // New thingie - now assume Lambda lives only on sites with checkerboard 
    // CB
    //            Lambda
    //   0           X           0           x = cb, O = 1-cb
    //
    //
    // Lambda                 Lambda 
    //   X           0           X
    //
    //
    //            Lambda 
    //   0           X           0
    //
    // So I can only construct 4 out of the 8 staples on the sites
    // that have CB and the OTHER 4 of the 8 staples on sites with 
    // 1-cb
    //
    
    // Sites with CB first:
    //

    U staple_for;
    U staple_back;
    U staple_left;
    U staple_right;

    U u_nu_for_mu = shift(u[nu],FORWARD, mu); // Can reuse these later
    U u_mu_for_nu = shift(u[mu],FORWARD, nu);
    U Lambda_xplus_mu = shift(Lambda, FORWARD, mu);
    U Lambda_xplus_nu = shift(Lambda, FORWARD, nu);
    U Lambda_xplus_muplusnu = shift(Lambda_xplus_mu, FORWARD, nu);

    U u_tmp3;

    U ds_tmp_mu;
    U ds_tmp_nu;
    {
      U up_left_corner;
      U up_right_corner;
      U low_right_corner;
      U low_left_corner;

      //   u_tmp1 =   <-------
      //              |
      //              |                       ON ALL CHECKERBOARDS
      //              |                       (Because it's used in staples)  
      //              V
      up_left_corner = adj(u_mu_for_nu)*adj(u[nu]);
      
      
      //
      //              <------^
      //                     |
      //                     |
      //                     |
      
      up_right_corner = u_nu_for_mu*adj(u_mu_for_nu);
      
      //                       |
      //                       |
      //                       |
      //                       V
      //                <------
      low_right_corner = adj(u_nu_for_mu)*adj(u[mu]);
      
      //
      //                    ^
      //  low left corner=  |                         ON ALL CHECKBERBOARDS
      //                    |                         (Because it's used in the staples)
      //                    |  
      //                     <-------
      low_left_corner = adj(u[mu])*u[nu];
      
      
      // Now compute the terms of the force:
      // 
      // Altogether 8 terms. 4 Upwards with + sign, and 4 Downwards with - sign
      //                     4 terms use staples and 4 don't
      
      // NON STAPLE TERMS FIRST:
      
      // 1) mu links
      //
      //    <-------  X (CB)                      <--------
      //    |         ^                           |
      //    |         |        re  use  u_tmp1 =  | 
      //    V         |                           V
      //   CB       1-CB       
      u_tmp3[rb[cb]] = u_nu_for_mu*Lambda_xplus_muplusnu;
      ds_u_mu[rb[cb]] = u_tmp3*up_left_corner;
      
      //    nu links
      //    X
      //     <------
      //     | 
      //     |
      //     |
      //     V-----> CB
      //   (1-CB)   
      //
      u_tmp3[rb[1-cb]] = adj(u_mu_for_nu)*Lambda_xplus_nu;
      
      // accumulate into ds_tmp_nu and shift everything together at the end
      ds_tmp_nu[rb[1-cb]] = u_tmp3*adj(low_left_corner);
      
      
      
      // 2)  mu links
      //
      //  CB
      //    X <------ 
      //    |        ^       re use u[nu](x+mu) = u_nu_for_mu
      //    |        |       re use u[mu](x+nu) = u_mu_for_nu
      //    V        |
      //    1-CB    CB        
      u_tmp3[rb[1-cb]] = Lambda_xplus_nu*adj(u[nu]);
      ds_u_mu[rb[1-cb]] = up_right_corner * u_tmp3;
      
      //      nu_links
      //    
      //     <------
      //     | 
      //     |
      //     |
      //   X V----->1-CB
      //   (CB)   
      //
      u_tmp3[rb[cb]] = up_left_corner*Lambda;
      //
      // accumulate into ds_tmp_nu and shift everything together at the end
      ds_tmp_nu[rb[cb]] = u_tmp3*u[mu];
      
      
      
      
      //
      // Terms 3) and 4)
      //
      // These last two can be done on the other checkerboard and then shifted together. at the very end...
      //
      //  CB      1-CB          
      //    ^       |   
      //    |       |     
      //    |       V
      //    <-------X CB
      
      
      // 3) Mu links
      //
      //  Compunte with low_left_corner:     ^           |
      //                                     |           | 
      //               low_left_corner    =  |           | 
      //                                     |           V
      //                              (1-CB) <--------   X CB
      u_tmp3[rb[1-cb]] = adj(u_nu_for_mu)*Lambda_xplus_mu;
      //
      // accumulate into ds_tmp_mu and shift at the end.
      ds_tmp_mu[rb[1-cb]] = u_tmp3*low_left_corner;
      
      // Nu links
      // 
      //  CB    ------>                                     ------->
      //                |                                          |
      //                |      reuse adj(up_right_corner):         |
      //                |                                          |
      //                V                                          V
      //  1-CB   <------X
      u_tmp3[rb[1-cb]] = adj(up_right_corner)*Lambda_xplus_mu;
      ds_u_nu[rb[1-cb]] = u_tmp3*adj(u[mu]);
      
      
      
      // 4) Mu links
      //
      //  1-CB      CB
      //   ^        |
      //   |        |        reuse = u[nu](x+mu) = u_nu_for_mu
      //   |        |
      //   X <----- V 1-CB
      //   CB
      u_tmp3[rb[cb]] = low_right_corner*Lambda;
      //
      // accumulate into ds_tmp_mu and shift at the end.
      ds_tmp_mu[rb[cb]] = u_tmp3*u[nu];
      
      
      
      // Nu links
      // 
      //  1-CB   ------> X                               
      //               |                                          |
      //               |       reuse low_right_corner:            |
      //               |                                          |
      //               V                                          V
      //   CB   <------                                    <------
      u_tmp3[rb[cb]] =    u_mu_for_nu*Lambda_xplus_muplusnu;
      ds_u_nu[rb[cb]] =   u_tmp3*low_right_corner;
      
      
      //  ds_tmp_mu now holds the last 2 terms, one on each of its checkerboards, but Now I need
      //  to shift them both together onto ds_u_mu
      //  I'll keep them in ds_tmp_mu right, bearing in mind I'll need to bring
      //  them in with a -ve contribution...
      

      // STAPLE TERMS:   
      
      // Construct the staples

      //  Staple_for =  <--------
      //                |       ^
      //                |       |             ON ALL CHECKERBOARDS
      //                |       |             
      //                V       |
      staple_for = u_nu_for_mu*up_left_corner;


      // Staple_right =   <-----             ON ALL CHECKERBOARDS
      //                 |
      //                 |
      //                 V
      //                 ----->
      staple_right = up_left_corner*u[mu];
      
      
      
      //                 ----->
      //                       |
      //                       |
      //                       |
      //                <----- V
      staple_left  = u_mu_for_nu*low_right_corner;
      
      
      
      
      //  Staple_back =  ^       |
      //                 |       |            ON ALL CHECKERBOARDS
      //                 |       |
      //                 <------ V
      //                      
      staple_back = adj(u_nu_for_mu)*low_left_corner;
      
    }  // Corner pieces go away here
    
    // 5) Mu links
    //
    //    <------- 
    //    |        ^
    //    |        |     use computed staple
    //    V        |
    //    x        
    //   CB       1-CB  
    ds_u_mu[rb[cb]] += staple_for*Lambda;


    //  Nu links
    //
    //     CB   <---- 1-CB
    //        |
    //        |                use staple_right
    //        V
    //    1-CB  -----> X CB
    //
    //  Accumulate into ds_tmp_nu and shift at the end.

    ds_tmp_nu[rb[1-cb]] += staple_right*Lambda_xplus_mu;


    // 6)  Mu links
    //
    //    <------- 
    //    |        ^
    //    |        |    re use computed staple 
    //    V        |
    //   1-CB      X CB	  

    ds_u_mu[rb[1-cb]] += Lambda_xplus_mu*staple_for;



    //  Nu links
    // 
    //      <----  X CB
    //     |
    //     |                     use adj(staple_right)
    //     |
    // CB  V ----> (1-CB)
    ds_tmp_nu[rb[cb]] += Lambda_xplus_muplusnu * staple_right;


    // 7) Mu links
    //
    //   CB      1-CB
    //  X         
    //    ^       |
    //    |       |   re use computed staple 
    //    |       |
    //    <-------V 
    //
    //  Accumulate into ds_tmp_mu and shift at the end.
    ds_tmp_mu[rb[1-cb]] += staple_back*Lambda_xplus_nu;

    //   Now for nu
    //
    //   (1-CB)  -----> CB         use adj(staple_left)
    //                |
    //                |
    //                V
    //      CB X <----
    //
    ds_u_nu[rb[cb]] += staple_left*Lambda;

    // 8) Mu links
    //
    //  1-CB      X CB
    //    ^       |
    //    |       |  reuse computed staple 
    //    |       |
    //    <-------V
    //
    // Accumulate into ds_tmp_mu and shift at the end
    ds_tmp_mu[rb[cb]] += Lambda_xplus_muplusnu * staple_back;

    // Now for Nu
    // 
    //    CB X ------> (1-CB)
    //               |
    //               |
    //               |
    //               V
    // 1-CB  <------- CB
    ds_u_nu[rb[1-cb]] += Lambda_xplus_nu * staple_left;

    // Now shift the accumulated pieces to mu and nu
    // 
    // Hope that this is not too slow as an expression
    ds_u_mu -= shift(ds_tmp_mu, BACKWARD, nu);
    ds_u_nu -= shift(ds_tmp_nu, BACKWARD, mu); 

    END_CODE();
  }


  //! Take deriv of D
  /*!
   * \param chi     left vector on cb                           (Read)
   * \param psi     right vector on 1-cb                        (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   * \param cb      Checkerboard of chi vector                  (Read)
   *
   * \return Computes   \f$\chi^\dag * \dot(D} * \psi\f$
   */
  template<typename T, typename U>
  void CloverTermBase<T,U>::deriv(multi1d<U>& ds_u, 
			     const T& chi, const T& psi, 
			     enum PlusMinus isign, int cb) const
  {
    START_CODE();


    // Do I still need to do this?
    if( ds_u.size() != Nd ) { 
      ds_u.resize(Nd);
    }

    ds_u = zero;

    // Get the links
    const multi1d<U>& u = getU();


    // Now compute the insertions
    for(int mu=0; mu < Nd; mu++) {
      for(int nu = mu+1; nu < Nd; nu++) {
	
	// These will be appropriately overwritten - no need to zero them.
	// Contributions to mu links from mu-nu clover piece
	U ds_tmp_mu; 

	// -ve contribs  to the nu_links from the mu-nu clover piece 
	// -ve because of the exchange of gamma_mu gamma_nu <-> gamma_nu gamma_mu
	U ds_tmp_nu;

	// The weight for the terms
	Real factor = (Real(-1)/Real(8))*getCloverCoeff(mu,nu);

	// Get gamma_mu gamma_nu psi -- no saving here, from storing shifts because
	// I now only do every mu, nu pair only once.

	int mu_nu_index = (1 << mu) + (1 << nu); // 2^{mu} 2^{nu}
	T ferm_tmp = Gamma(mu_nu_index)*psi;
	U s_xy_dag = traceSpin( outerProduct(ferm_tmp,chi));
	s_xy_dag *= Real(factor);

	// Compute contributions
	deriv_loops(mu, nu, cb, ds_tmp_mu, ds_tmp_nu, s_xy_dag);

	// Accumulate them
	ds_u[mu] += ds_tmp_mu;
	ds_u[nu] -= ds_tmp_nu;


      }
    }


    // Clear out the deriv on any fixed links
    (*this).getFermBC().zero(ds_u);
    END_CODE();
  }

  template<typename T, typename U>
  void CloverTermBase<T,U>::derivMultipole(multi1d<U>& ds_u, 
					   const multi1d<T>& chi, const multi1d<T>& psi, 
					   enum PlusMinus isign, int cb) const
  {
    // Multipole deriv
    START_CODE();
    
    // Do I still need to do this?
    if( ds_u.size() != Nd ) { 
      ds_u.resize(Nd);
    }

    ds_u = zero;

    // Get the links
    const multi1d<U>& u = getU();


    // Now compute the insertions
    for(int mu=0; mu < Nd; mu++) {
      for(int nu = mu+1; nu < Nd; nu++) {
	
	// These will be appropriately overwritten - no need to zero them.
	// Contributions to mu links from mu-nu clover piece
	U ds_tmp_mu; 

	// -ve contribs  to the nu_links from the mu-nu clover piece 
	// -ve because of the exchange of gamma_mu gamma_nu <-> gamma_nu gamma_mu
	U ds_tmp_nu;

	// The weight for the terms
	Real factor = (Real(-1)/Real(8))*getCloverCoeff(mu,nu);

	// Get gamma_mu gamma_nu psi -- no saving here, from storing shifts because
	// I now only do every mu, nu pair only once.

	int mu_nu_index = (1 << mu) + (1 << nu); // 2^{mu} 2^{nu}

	// Accumulate all the trace spin outer products 
	U s_xy_dag = zero;
	for(int i=0; i < chi.size(); i++) { 
	  T ferm_tmp = Gamma(mu_nu_index)*psi[i];
	  s_xy_dag += traceSpin( outerProduct(ferm_tmp,chi[i]));
	}

	s_xy_dag *= Real(factor);

	// Compute contributions
	deriv_loops(mu, nu, cb, ds_tmp_mu, ds_tmp_nu, s_xy_dag);

	// Accumulate them
	ds_u[mu] += ds_tmp_mu;
	ds_u[nu] -= ds_tmp_nu;


      }
    }


    // Clear out the deriv on any fixed links
    (*this).getFermBC().zero(ds_u);
    END_CODE();
  }


  //! Take deriv of D using Trace Log
  /*!
   * \param chi     left vector on cb                           (Read)
   * \param psi     right vector on 1-cb                        (Read)
   * \param isign   D'^dag or D'  ( MINUS | PLUS ) resp.        (Read)
   * \param cb      Checkerboard of chi vector                  (Read)
   *
   * \return Computes   \f$\chi^\dag * \dot(D} * \psi\f$  
   */
  template<typename T, typename U>
  void CloverTermBase<T,U>::derivTrLn(multi1d<U>& ds_u, 
				 enum PlusMinus isign, int cb) const
  {
    START_CODE();
    
    // Do I still need to do this?
    if( ds_u.size() != Nd ) { 
      ds_u.resize(Nd);
    }
    
    ds_u = zero;

    for(int mu=0; mu < Nd; mu++) {
      for(int nu = mu+1; nu < Nd; nu++) { 

	  // Index 
	  int mu_nu_index = (1 << mu) + (1 << nu); // 2^{mu} 2^{nu}

	  // The actual coefficient factor
	  Real factor = Real(-1)*getCloverCoeff(mu,nu)/Real(8);
	  
	  U sigma_XY_dag=zero;

	  // Get  weight*Tr_spin gamma_mu gamma_nu A^{-1} piece
	  triacntr(sigma_XY_dag, mu_nu_index, cb);
	  sigma_XY_dag[rb[cb]] *= factor;

	  // These will be overwritten so no need to initialize to zero
	  U ds_tmp_mu;
	  U ds_tmp_nu;

	  // Get contributions from the loops and insersions
	  deriv_loops(mu, nu, cb, ds_tmp_mu, ds_tmp_nu, sigma_XY_dag);

	  // Accumulate
	  ds_u[mu] += ds_tmp_mu;
	  // -ve weight for nu from gamma_mu gamma_nu -> gamma_nu gamma_mu
	  // commutation.
	  ds_u[nu] -= ds_tmp_nu;

      } // End loop over nu

    } // end of loop over mu
    

    // Not sure this is needed here, but will be sure
    (*this).getFermBC().zero(ds_u);
    
    END_CODE();
  }
      


} // End Namespace Chroma


#endif
