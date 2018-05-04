// -*- C++ -*-
/*! \file
 *  \brief Solve a M*psi=chi linear system by FGMRES_DR
 */

#ifndef __syssolver_linop_fgmres_dr_h__
#define __syssolver_linop_fgmres_dr_h__

#include "chroma_config.h"
#include "handle.h"
#include "state.h"
#include "syssolver.h"
#include "linearop.h"

#include "util/gauge/reunit.h"
#include "actions/ferm/invert/syssolver_linop.h"
#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_fgmres_dr_params.h"

namespace Chroma
{

  //! FGMRESDR system solver namespace
  namespace LinOpSysSolverFGMRESDREnv
  {
    //! Register the syssolver
    bool registerAll();
  }


  class Givens {
  public:

    // Givens rotation.
    //   There are a variety of ways to choose the rotations 
    //   which can do the job. I employ the method given by Saad
    //   in Iterative Methods (Sec. 6.5.9 (eq. 6.80)
    //
    //  [  conj(c) conj(s)  ] [ h_jj    ] = [ r ]
    //  [    -s       c     ] [ h_j+1,j ]   [ 0 ]
    // 
    //  We know that h_j+1,j is a vector norm so imag(h_{j+1,j} = 0)
    //
    //  we have: s = h_{j+1,j} / t
    //           c = h_jj / t
    //
    //   t=sqrt( norm2(h_jj) + h_{j+1,j}^2 ) is real and nonnegative
    //
    //  so in this case s, is REAL and nonnegative (since t is and h_{j+1,j} is
    //  but  c is in general complex valued.
    //
    // 
    //  using this we find r = conj(c) h_jj + conj(s) h_{j+1,j}
    //                        = (1/t)* [  conj(h_jj)*h_jj + h_{j+1,j}*h_{j+1,h} ]
    //                        = (1/t)* [  norm2(h_jj) + h_{j+1,j}^2 ]
    //                        = (1/t)* [ t^2 ] = t 
    //
    //  Applying this to a general 2 vector 
    //    
    //   [ conj(c) conj(s) ] [ a ] = [ r_1 ]
    //   [   -s      c     ] [ b ]   [ r_2 ]
    //
    //   we have r_1 = conj(c)*a + conj(s)*b  = conj(c)*a + s*b  since s is real, nonnegative
    //      and  r_2 = -s*a + c*b 
    //
    //  NB: In this setup we choose the sine 's' to be real in the rotation.
    //      This is in contradistinction from LAPACK which typically chooses the cosine 'c' to be real
    //
    //
    // There are some special cases: 
    //   if  h_jj and h_{j+1,j} are both zero we can choose any s and c as long as |c|^2 + |s|^2 =1
    //   Keeping with the notation that s is real and nonnegative we choose
    //   if    h_jj != 0 and h_{j+1,h) == 0 => c = sgn(conj(h_jj)), s = 0, r = | h_jj |
    //   if    h_jj == 0 and h_{j+1,j} == 0 => c = 0, s = 1,  r = 0 = h_{j+1,j} 
    //   if    h_jj == 0 and h_{j+1,j} != 0 => c = 0, s = 1,  r = h_{j+1,j} = h_{j+1,j}
    //   else the formulae we computed.

    /*! Given  a marix H, construct the rotator so that H(row,col) = r and H(row+1,col) = 0
     *
     *  \param col  the column Input
     *  \param  H   the Matrix Input
     */

    Givens(int col, const multi2d<DComplex>& H) : col_(col) 
    {
      DComplex f = H(col_,col_);
      DComplex g = H(col_,col_+1);

      if( toBool( real(f) == Double(0) && imag(f) == Double(0) ) ) {
       
	// h_jj is 0
	c_ = DComplex(0);
	s_ = DComplex(1);
	r_ = g;  // Handles the case when g is also zero
      }
      else {
	if( toBool( real(g) == Double(0) && imag(g) == Double(0) ) ) {
	  s_ = DComplex(0);
	  c_ = conj(f)/sqrt( norm2(f) ); //   sgn( conj(f) ) = conj(f) / | conj(f) |  = conj(f) / | f | 
	  r_ = DComplex( sqrt(norm2(f)) );
	}
	else {
	  // Revisit this with 
	  Double t = sqrt( norm2(f) + norm2(g) );
	  r_ = t;
	  c_  = f/t;
	  s_  = g/t;
	}
      }
    }
    
    /*! Apply the rotation to column col of the matrix H. The 
     *  routine affects col and col+1.
     *
     *  \param col  the columm
     *  \param  H   the matrix
     */

    void operator()(int col,  multi2d<DComplex>& H) {
      if ( col == col_ ) {
	// We've already done this column and know the answer
	H(col_,col_) = r_;
	H(col_,col_+1) = 0;
      }
      else {
	int row = col_; // The row on which the rotation was defined
	DComplex a = H(col,row);
	DComplex b = H(col,row+1);
	H(col,row) = conj(c_)*a + conj(s_)*b;
	H(col,row+1) = -s_*a + c_*b;
      }
    }

    /*! Apply rotation to Column Vector v */
    void operator()(multi1d<DComplex>& v) {
      	DComplex a =  v(col_);
	DComplex b =  v(col_+1);
	v(col_) = conj(c_)*a + conj(s_)*b;
	v(col_+1) =  -s_*a + c_*b;

    }

  private:
    int col_;
    DComplex s_;
    DComplex c_;
    DComplex r_; 
  };


  //! Solve a M*psi=chi linear system by FGMRESDR
  /*! \ingroup invert
   */

  class LinOpSysSolverFGMRESDR : public LinOpSystemSolver<LatticeFermion>
  {
  public:
    using T = LatticeFermion;
    using U = LatticeColorMatrix;
    using Q = multi1d<U>;
 
    //! Constructor
    /*!
     * \param A_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    LinOpSysSolverFGMRESDR(Handle< LinearOperator<T> > A,
			   Handle< FermState<T,Q,Q> > state,
			   const SysSolverFGMRESDRParams& invParam) : 
      A_(A), state_(state), invParam_(invParam)
      {
	  
	// Use factory to construct preconditioner
	// The preconditioner has to look like a system solver for now.
	try { 
	  std::istringstream precond_xml(invParam_.PrecondParams.xml );
	  XMLReader precond_xml_reader(precond_xml);

	  preconditioner_ = TheLinOpFermSystemSolverFactory::Instance().createObject( invParam_.PrecondParams.id, precond_xml_reader, invParam_.PrecondParams.path, state_, A_);
	}
	catch(...) { 
	  QDPIO::cout << "Unable to create Preconditioner" << std::endl;
	  QDP_abort(1);
	}

#ifndef BUILD_LAPACK
	QDPIO::cout << "WARNING: LAPACK Is not built!!!" << std::endl;
	QDPIO::cout << "WARNING: SUBSPACE Augmentation requires it!!!!" << std::endl;
	QDPIO::cout << "WARNING: Using the solver with NDefl > 0 will fail!!! " << std::endl;
	QDPIO::cout << "WARNING: Run this version with NDefl = 0 only!!!!" << std::endl;
	QDPIO::cout << "WARNING: OR reconfigure with --enable-lapack=lapack" << std::endl;

	if( invParam_.NDefl > 0 ) { 
	  QDPIO::cout << " invParam.NDefl > 0: This mode is not supported without LAPACK" <<std::endl;
	  QDP_abort(1);
	}
#endif
	  

	// Initialize stuff
	InitMatrices();

	
      }


    //! Initialize the internal matrices
    void InitMatrices();
    

    //! Destructor is automatic
    ~LinOpSysSolverFGMRESDR() {}
    
    //! Return the subset on which the operator acts
    const Subset& subset() const {return A_->subset();}

    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     */
    SystemSolverResults_t operator() (T& psi, const T& chi) const;


    void FlexibleArnoldi(int n_krylov,
			 int n_deflate,
			 const Real& rsd_target,
			 multi1d<T>& V,
			 multi1d<T>& Z,
			 multi2d<DComplex>& H, 
			 multi2d<DComplex>& R,
			 multi1d< Handle<Givens> >& givens_rots,
			 multi1d<DComplex>& g,
			 multi2d<DComplex>& Qk, 
			 multi1d<DComplex>& Qk_tau,
			 int&  ndim_cycle) const;

    void LeastSquaresSolve(const multi2d<DComplex>& R, 
			   const multi1d<DComplex>& rhs,
			   multi1d<DComplex>& eta,
			   int dim) const;

    void GetEigenvectors(int total_dim, 
			 const multi2d<DComplex>& H,
			 multi1d<DComplex>& f_m,
			 multi2d<DComplex>& evecs,
			 multi1d<DComplex>& evals,
			 multi1d<int>& order_array) const;

  private:
    // Hide default constructor
    LinOpSysSolverFGMRESDR() {}
    Handle< LinearOperator<T> > A_;
    Handle< FermState<T,Q,Q> > state_;
    SysSolverFGMRESDRParams invParam_;
    Handle< LinOpSystemSolver<T> > preconditioner_;

    // These can become state variables, as they will need to be
    // handed around
    mutable multi2d<DComplex> H_; // The H matrix
    mutable multi2d<DComplex> R_; // R = H diagonalized with Givens rotations
    mutable multi1d<T> V_;  // K(A)
    mutable multi1d<T> Z_;  // K(MA)
    mutable multi1d< Handle<Givens> > givens_rots_;

    // This is the c = V^H_{k+1} r vector (c is frommers Notation)
    // For regular FGMRES I need to keep only the basis transformed
    // version of it, for rotating by the Q's (I call this 'g')
    // However, for FGMRES-DR I need this because I need
    //  || c - H eta || to extend the G_k matrix to form V_{k+1}
    // it is made of the previous v_{k+1} by explicitly evaluating c
    // except at the end of the first cycle when the V_{k+1} are not
    // yet available

    // Once G_k+1 is available, I can form the QR decomposition
    // and set  my little g = Q^H and then as I do the Arnoldi 
    // I can then work it with the Givens rotations...
    
    mutable multi1d<DComplex> c_;                                 
    mutable multi1d<DComplex> eta_;
    mutable multi1d<DComplex> g_;
    mutable multi2d<DComplex> Hk_QR_;
    mutable multi1d<DComplex> Hk_QR_taus_;


  };

} // End namespace

#endif 

