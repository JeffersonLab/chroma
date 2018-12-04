// -*- C++ -*-
/*! \file
 *  \brief Symmetric even-odd preconditioned Clover fermion linear operator
 *  \with shifted mass term i*\mu*\gamma5*A_oo
 */

#ifndef __shiftedseo_clover_linop_w_h__
#define __shiftedseo_clover_linop_w_h__

#include "state.h"
#include "fermbc.h"
#include "seoprec_logdet_linop.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/clover_term_w.h"

namespace Chroma
{
  //! Symmetric even-odd preconditioned Clover-Dirac operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   *
   * The kernel for Clover fermions is
   *
   *      M  =  A + (d+M) - (1/2) D'
   */
	class ShiftSymEvenOddPrecCloverLinOp: 
		public SymEvenOddPrecLogDetLinearOperator<LatticeFermion, 
		multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
	{
		public:
			// Typedefs to save typing
			typedef LatticeFermion               T;
			typedef multi1d<LatticeColorMatrix>  P;
			typedef multi1d<LatticeColorMatrix>  Q;

			//! Partial constructor
			ShiftSymEvenOddPrecCloverLinOp(){}

			//! Full constructor with shifted mass mu
		    ShiftSymEvenOddPrecCloverLinOp(Handle< FermState<T,P,Q> > fs,
					      const CloverFermActParams& param_, const Real& mu_):mu(mu_)
		    {
		      create(fs,param_);
		    }
		
		    //! Destructor is automatic
		    ~ShiftSymEvenOddPrecCloverLinOp() {
		    }
		
		    //! Return the fermion BC object for this linear operator
		    const FermBC<T,P,Q>& getFermBC() const {return D.getFermBC();}
		
		    //! Creation routine
		    void create(Handle< FermState<T,P,Q> > fs,
				const CloverFermActParams& param_);
		
		    //! Apply the the even-even block onto a source std::vector
		    void unprecEvenEvenLinOp(T& chi, const T& psi,
				       enum PlusMinus isign) const;
		
		    //! Apply the inverse of the even-even block onto a source std::vector
		    void unprecEvenEvenInvLinOp(T& chi, const T& psi,
					  enum PlusMinus isign) const;
		  
		    //! Apply the the odd-odd block onto a source std::vector
		     void unprecOddOddLinOp(T& chi, const T& psi,
		 		     enum PlusMinus isign) const;
		
		     //! Apply the inverse of the odd-odd block onto a source std::vector
		     void unprecOddOddInvLinOp(T& chi, const T& psi,
		 			enum PlusMinus isign) const;
		
		    //! Apply the the even-odd block onto a source std::vector
		    void unprecEvenOddLinOp(T& chi, const T& psi,
				      enum PlusMinus isign) const ;
		
		    //! Apply the the odd-even block onto a source std::vector
		    void unprecOddEvenLinOp(T& chi, const T& psi,
				      enum PlusMinus isign) const ;

			// Deriv of A_ee
			virtual void derivUnprecEvenEvenLinOp(P& ds_u, const T& chi, const T& psi, enum PlusMinus isign) const;

			// Deriv of A_oo
			virtual void derivUnprecOddOddLinOp(P& ds_u, const T& chi, const T& psi,
					enum PlusMinus isign) const;

			// Deriv of D_eo
			virtual void derivUnprecEvenOddLinOp(P& ds_u, const T& chi, const T& psi,
					enum PlusMinus isign) const;

			// Deriv of D_oe
			virtual void derivUnprecOddEvenLinOp(P& ds_u, const T& chi, const T& psi,
					enum PlusMinus isign) const;

			//! Apply the even-even block onto a source std::vector
			void derivLogDetEvenEvenLinOp(P& ds_u,
					enum PlusMinus isign) const;

			//! Apply the odd-odd block onto a source std::vector
			void derivLogDetOddOddLinOp(P& ds_u,
					enum PlusMinus isign) const;
			
			void operator()(T& chi, const T& psi,
					enum PlusMinus isign) const override;

			void deriv(P& ds_u, const T& chi, const T& psi,
					enum PlusMinus isign) const override;

			void derivMultipole(P& ds_u, const multi1d<T>& chi, 
					const multi1d<T>& psi,
					enum PlusMinus isign) const override;

			//! Return flops performed by the operator()
			unsigned long nFlops() const;

			//! Get the log det of the even even part
			Double logDetEvenEvenLinOp(void) const;

			//! Get the log det of the odd odd part
			Double logDetOddOddLinOp(void) const;

		private:
			CloverFermActParams param;
			WilsonDslash D;
			CloverTerm clov;
			CloverTerm invclov;
			Real mu; // shifted mass 
			mutable double clov_apply_time;
			mutable double clov_deriv_time;
			mutable StopWatch swatch;


	};

} // End NameSpace Chroma

#endif
