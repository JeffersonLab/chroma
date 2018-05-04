/*! \file
 *  \brief Solve a M*psi=chi linear system by MR
 */
#include <algorithm>
#include <functional>
#include <vector>
#include "chromabase.h"
#include "qdp-lapack.h"
#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"

#include "actions/ferm/invert/syssolver_linop_fgmres_dr.h"

namespace Chroma
{

  //! FGMRESDR system solver namespace
  namespace LinOpSysSolverFGMRESDREnv
  {
    //! Callback function
    LinOpSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new LinOpSysSolverFGMRESDR(A,state, SysSolverFGMRESDRParams(xml_in, path));
    }


    //! Name to be used
    const std::string name("FGMRESDR_INVERTER");

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheLinOpFermSystemSolverFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }
  }


  /*! Flexible Arnodli Iteration.
   *
   * Works on the system A (dx) = r
   *
   * resulting from recasting  A ( x_0 + dx ) = b
   *
   * This routine, produces the Krylov Subspace V
   * and the Flexible subspace Z
   *
   * It also constructs the relevant H matrix
   * and currently diagonalizes it with Givens rotations to yield R
   *
   * as it is proceeding it builds up the vector g, and tracks the accumulated
   * residuum. 
   *  
   *  The workspaces: V, Z, H, R, givens_rots and 'g' must be initialized
   *  outside this routine. This routine does no size checking at the moment.
   *
   *
   *
   * \param  n_krylov is the length of the 'cycle'
   * \param  n_deflate is the size of the augmented subspace ( not yet used )
   * \param  rsd_target is the desired target relative residuum
   * \param  A is the linear operator for the solve
   * \param  M is the preconditioner (a SystemSolver)

   */
  template< typename T >
  void FlexibleArnoldiT(int n_krylov, 
			int n_deflate,                  // Size of augmented space. Not yet used.
			const Real& rsd_target,
			const LinearOperator<T>& A,     // Operator
			const LinOpSystemSolver<T>& M,  // Preconditioner
			multi1d<T>& V,
			multi1d<T>& Z,
			multi2d<DComplex>& H,
			multi2d<DComplex>& R,
			multi1d< Handle<Givens> >& givens_rots,
			multi1d<DComplex>& g,
			multi2d<DComplex>& Qk,
			multi1d<DComplex>& Qk_tau,
			int& ndim_cycle)

  {
    const Subset& s = M.subset();      // Linear Operator Subset
    ndim_cycle = 0;     

    const int total_dim=n_krylov+n_deflate;
    

    // Work by columns:
    for(int j=n_deflate; j < total_dim; ++j) {
     
      M( Z[j], V[j] );  // z_j = M v_j
      T w;
      A( w, Z[j], PLUS);  // w  = A z_

      // Fill out column j
      for(int i=0; i <= j ;  ++i ) {
	H(j,i) = innerProduct(V[i], w, s);
	w[s] -= H(j,i)* V[i];
      }
      
      Double wnorm=sqrt(norm2(w,s));
      H(j,j+1) = DComplex(wnorm);

      // In principle I should check w_norm to be 0, and if it is I should
      // terminate: Code Smell -- get rid of 1.0e-14.
      if ( toBool( fabs( wnorm ) < Double(1.0e-14) ) ) {

	// If wnorm = 0 exactly, then we have converged exactly
	// Replay Givens rots here, how to test?

	QDPIO::cout << "Converged at iter = " << j-n_deflate+1 << std::endl;
	ndim_cycle = j;
	return;
      }

      Double invwnorm = Double(1)/wnorm;
      V[j+1] = invwnorm*w;

      // Done construcing H.
      // Now to incrementally update the R matrix
      // and compute the residuum
      if( n_deflate > 0 ) { 

	// If we have deflation, Qk and Qk_tau
	// should hold the QR decomposition of the 
	// Deflated part of H, called H_k that is constructed
	// outside the Arnoldi.
	// To incrementally update the 'R' matrix we need to 
	// apply the Q out of this, to the current column, before
	// adding the Givens rotations
	
	// Copy column of H into R_col
	multi1d<DComplex> R_col(j+2);
	for(int i=0; i <= j+1; i++) { 
	  R_col[i] = H(j,i);
	}

	// R_col <- Qk^H R_col
	char side = 'L';
        char trans = 'C';
	QDPLapack::zunmqrv(side,trans, n_deflate+1, n_deflate+1, Qk, Qk_tau, R_col );
	// Copy into R
	for(int i=0; i <= j+1; i++) { 
	  R(j,i)=R_col[i];
	}
      }
      else {
	// If no deflation, just copy the column into R
	for(int i=0; i <= j+1; ++i) { 
	  R(j,i) = H(j,i);
	}
      }

      // Apply Existing Givens Rotations to this column of R
      for(int i=0;i < j; ++i) {
	(*givens_rots[i])(j,R);
      }

      // Compute next Givens Rot for this column
      givens_rots[j] = new Givens(j,R); 

      (*givens_rots[j])(j,R); // Apply it to R
      (*givens_rots[j])(g);   // Apply it to the g vector

      Double accum_resid = fabs(real(g[j+1]));

#if 0
      // Debug
      for(int rows=0; rows <= j+1; ++rows) {
	QDPIO::cout << " g["<<rows<<"]="<<g[rows]<< std::endl;
      }
#endif 
      // j-ndeflate is the 0 based iteration count
      // j-ndeflate+1 is the 1 based human readable iteration count
      QDPIO::cout << "Iter " << j-n_deflate+1 << " || r || = " << accum_resid << " Target=" <<rsd_target << std::endl;

      ndim_cycle = j+1;
      if ( toBool( accum_resid <= rsd_target ) ) { 
      	QDPIO::cout << "Flexible Arnoldi Cycle Converged at iter = " << j-n_deflate << std::endl;
	return;
      }
    }
  }

  /*! Initialize the key matrices (resize and zero)
   */
  void LinOpSysSolverFGMRESDR::InitMatrices()
  {
    int total_dim = invParam_.NKrylov + invParam_.NDefl;
    H_.resize(total_dim, total_dim+1); // This is odd. Shouldn't it be
    R_.resize(total_dim, total_dim+1); // resize (n_cols, n_rows)?
                                      // Doing that give weird double free 
                                      //errors

    V_.resize(total_dim+1);
    Z_.resize(total_dim+1);
    givens_rots_.resize(total_dim+1);
    g_.resize(total_dim+1);
    c_.resize(total_dim+1);
    eta_.resize(total_dim);

    
    for(int col =0; col < total_dim; col++) { 
      for(int row = 0; row < total_dim+1; row++) { 
	H_(col,row) = zero;
	R_(col,row) = zero;
      }
    }

    for(int row = 0; row < total_dim+1; row++) { 
      V_[row] = zero;
      Z_[row] = zero;
      g_[row] = zero;
      c_[row] = zero;
    }

    for(int row = 0; row < total_dim; row++) { 
      eta_[row] = zero;
    }
    

  }


  /*! Solve least squares system to get the coefficients for 
   *  Updating the solution at the end of an GGMRES-DR Cycle.
   *  In this instantiation the system solved is
   *
   *    R \eta = rhs
   *
   * where R is the matrix H from GMRES, reduced to upper triangular
   * form by Givens rotations, and rhs is the result of the Givens rotations
   * applied to the orignal g = [ beta, 0, 0 ... ]^T vector.
   *
   * Hence this operation is a straightforward back substitution
   * NB: the matrix R has dim+1 rows and dim columns, but the dim+1-th row has 
   * All zero entries and should be deleted. d
   * 
   * dim is passed in instead of using matrix dimensions in case
   * the cycle terminated early.
   *
   */

  void
  LinOpSysSolverFGMRESDR::LeastSquaresSolve(const multi2d<DComplex>& R, 
					    const multi1d<DComplex>& Q_rhs,
					    multi1d<DComplex>& eta,
					    int n_cols) const
  {
    /* Assume here we have a square matrix with an extra row.
       Hence the loop counters are the columns not the rows.
       NB: For an augmented system this will change */
    eta[n_cols-1] = Q_rhs[n_cols-1]/R(n_cols-1,n_cols-1);
    for(int row = n_cols-2; row >= 0; --row) {
      eta[row] = Q_rhs[row];
      for(int col=row+1; col <  n_cols; ++col) { 
	eta[row] -= R(col,row)*eta[col];
      }
      eta[row] /= R(row,row);
    }
  }

#if 0 
  // This version for no incremental update of R. It takes H and c
  // and does an explicit QR on H. May be more stable?
  void
  LinOpSysSolverFGMRESDR::LeastSquaresSolve(const multi2d<DComplex>& H, 
					    const multi1d<DComplex>& rhs,
					    multi1d<DComplex>& eta,
					    int n_cols) const
  {
    // QR decompose H
    multi2d<DComplex> R(n_cols,n_cols+1);
    for(int col=0; col < n_cols; ++col) { 
      for(int row=0; row < n_cols+1; ++row) { 
	R(col,row) = H(col,row);
      }
    }

    multi1d<DComplex> Q_rhs(n_cols+1);
    for(int row=0; row < n_cols+1; ++row) { 
      Q_rhs[row] = rhs[row];
    }

    multi1d<DComplex> R_taus;
    QDPLapack::zgeqrf(n_cols+1, n_cols, R, R_taus);
    char side='L';
    char trans='C';
    QDPLapack::zunmqrv(side,
		       trans,
		       n_cols+1,
		       n_cols+1,
		       R, 
		       R_taus,
		       Q_rhs
		       );

    eta[n_cols-1] = Q_rhs[n_cols-1]/R(n_cols-1,n_cols-1);
    for(int row = n_cols-2; row >= 0; --row) {
      eta[row] = Q_rhs[row];
      for(int col=row+1; col <  n_cols; ++col) { 
	eta[row] -= R(col,row)*eta[col];
      }
      eta[row] /= R(row,row);
    }
  }
#endif

  /*! Get the eigenvectors of the matrix 
   *    H_tilde = H_m + f_m h^{H}_m 
   *
   *  of which we will keep the first NDefl corresponding
   *  to the NDefl eigenvalues of smallest modulus as our 
   *  matrix G. In turn the Harmonic Ritz vectors are 
   *  V_m G_m. 
   *
   *  We need a two stage process here. 
   *  One: first we need to solve H_m^H f_m = h^H_m
   *  where h_m is the last row of the \hat{H}_m
   *
   *  second, once in posession of f_m we need to constuct
   *  H_tilde and get its eigenvalues and eigenvectors.
   *
   *  NB: We use lapack here, specifically the routines:
   *   zgetrf/zgetrs to solve for f_m 
   *
   *   zgeev to get the eigenvalues/eigenvectors
   *  
   */
  void 
  LinOpSysSolverFGMRESDR::GetEigenvectors(int total_dim,
					  const multi2d<DComplex>& H,
					  multi1d<DComplex>& f_m,
					  multi2d<DComplex>& evecs,
					  multi1d<DComplex>& evals,
					  multi1d<int>& order_array) const
  {
    // Lapack routines overwrite the matrix 
    // with factorizations so I want a copy here.
    evals.resize(total_dim);
    evecs.resize(total_dim,total_dim);
    order_array.resize(total_dim);
    f_m.resize(total_dim);

    multi2d<DComplex> H_tilde(total_dim,total_dim);
    multi1d<DComplex> h_m(total_dim);
    for(int col=0; col < total_dim; ++col)  {
      for(int row=0; row < total_dim; ++row) { 
	H_tilde(col,row)=H(col,row);
      }
    }
    for(int col=0; col < total_dim; ++col) { 
      // h_m is a column vector and is the hermitian conjugate
      // of the last row of H (ie H_{row=total_dim}, col=:})
      // Hence I conjugate here.
      h_m[col] = conj(H(col,total_dim)); 
      f_m[col] = h_m[col];   // f_m will be overwritten
    }

    multi1d<int> ipiv(total_dim);
    int info;
    QDPLapack::zgetrf(total_dim,total_dim,H_tilde,total_dim,ipiv,info);
    if (info != 0) { 
      QDPIO::cout << "ZGETRF reported failure: info=" << info << std::endl;
      QDP_abort(1);
    }
    char trans='C'; // We will solve with H^H. C=solve with Herm. Conj.
    QDPLapack::zgetrs(trans,total_dim,1,H_tilde,total_dim,ipiv,f_m,total_dim,info);
    if (info != 0) { 
      QDPIO::cout << "ZGETRS reported failure: info=" << info << std::endl;
      QDP_abort(1);
    }

    // OK now we can construct the matrix whose eigenvalues we seek
    // We can reuse H_tilde for this
    for(int col=0; col < total_dim; ++col) { 
      for(int row=0; row < total_dim; ++row) { 
	H_tilde(col,row) = H(col,row) + f_m[row]*conj(h_m[col]);
      }
    }
     
    QDPLapack::zgeev(total_dim, H_tilde, evals, evecs);

    
    // Now I need to sort the eigenvalues in terms of smallest modulus.
    // I will use the C++ std library here.
    // It needs iterators tho.
    std::vector<int> std_order_array(total_dim);
    for(int i=0; i < total_dim; ++i) {
      std_order_array[i]=i;
    }
    
    // This calls the standard library sort function
    // with a C++ lambda for the comparison function.
    // 
    sort(std_order_array.begin(), std_order_array.end(),
	 // This is the lambda below: [&evals] brings evals into the closure
         // (so called capture region
	 [&evals](int i, int j)->bool {                
	   return toBool( norm2(evals[i]) < norm2(evals[j])); 
	 });

    for(int i=0; i < total_dim; ++i) {
      order_array[i]=std_order_array[i];
    }
    
  }

 

  /*! Flexible Arnolid Process. 
   *  Currently not using augmentation (n_deflate ignored)
   *
   *  NB: This currentl just forwards to a templated free function
   */
  void
  LinOpSysSolverFGMRESDR::FlexibleArnoldi(int n_krylov,
					  int n_deflate,
					  const Real& rsd_target,
					  multi1d<T>& V,
					  multi1d<T>& Z, 
					  multi2d<DComplex>&H,
					  multi2d<DComplex>&R,
					  multi1d<Handle<Givens> >& givens_rots,
					  multi1d<DComplex> &g,
					  multi2d<DComplex> &Qk,
					  multi1d<DComplex> &Qk_tau,
					  int& ndim_cycle) const
  {
    FlexibleArnoldiT<>(n_krylov, 
		       n_deflate,
		       rsd_target,
		       (*A_),
		       (*preconditioner_),
		       V,Z,H,R,givens_rots,g, Qk, Qk_tau, ndim_cycle);
  }


  /*! Solve the linear system  A psi = chi  via FGMRES-DR
   *  Right now the DR part is not implemented. 
   * 
   *
   *
   */

  SystemSolverResults_t 
  LinOpSysSolverFGMRESDR::operator() (T& psi, const T& chi) const
  {
    START_CODE();
    SystemSolverResults_t res; // Value to return
    
    const Subset& s = A_->subset();
    Double norm_rhs = sqrt(norm2(chi,s));   //  || b ||
    Double target = norm_rhs * invParam_.RsdTarget; // Target  || r || < || b || RsdTarget
 

    // Compute ||r|| 
    T r = zero; T tmp = zero;
    r[s] = chi;
    (*A_)(tmp, psi, PLUS);
    r[s] -=tmp;

    // The current residuum
    Double r_norm = sqrt(norm2(r,s));

    // Initialize iterations
    int iters_total = 0;
    int n_cycles = 0;
    int prev_dim = invParam_.NKrylov; // This is the default previous dimension

    // We are done if norm is sufficiently accurate, 
    bool finished = toBool( r_norm <= target ) ;

    // We keep executing cycles until we are finished
    while( !finished ) { 

      // If not finished, we should do another cycle with RHS='r' to find dx to add to psi.
      ++n_cycles;
      
  
      int dim; // dim at the end of cycle (in case we terminate in-cycle
      int n_deflate = invParam_.NDefl;
      int n_krylov  = invParam_.NKrylov;

      // Set up the input g vector (c-H eta?)
      // 
      if (n_cycles == 1 || n_deflate == 0 ) { 
	// We are either first cycle, or
	// We have no deflation subspace ie we are regular FGMRES
	// and we are just restarting
	//
	// Set up initial vector c = [ beta, 0 ... 0 ]^T
	// In this case beta should be the r_norm, ie || r || (=|| b || for the first cycle)
	//
	// NB: We will have a copy of this called 'g' onto which we will
	// apply Givens rotations to get an inline estimate of the residuum
	for(int j=0; j < g_.size(); ++j) { 
	  c_[j] = DComplex(0); 
	  g_[j] = DComplex(0);
	}
	c_[0] = r_norm;
	g_[0] = r_norm;
	
	// Set up initial V[0] = rhs / || r^2 ||
	// and since we are solving for A delta x = r
	// the rhs is 'r'
	//
	Double beta_inv = Double(1)/r_norm;
	V_[0][s] = beta_inv * r;
	
	n_deflate = 0; // Override n_deflate for first cycle
	               // So Arnoldi process starts from right place
      }
      else { 
	QDPIO::cout << "AUGMENTING SUBSPACE: n_deflate=" << n_deflate << " prev_dim=" << prev_dim << std::endl;

	// This is a cycle where we need to augment the space 
	multi1d<DComplex> f_m(prev_dim);
	multi2d<DComplex> evecs(prev_dim,prev_dim);
	multi1d<DComplex> evals(prev_dim);
	multi1d<int> order_array(prev_dim);
	
	(*this).GetEigenvectors(prev_dim, 
				H_,
				f_m,
				evecs,
				evals,
				order_array);

	
	// This is where we will store G_{k = [ g_1 | .. | g_k ]
	// 
	// G_{k+1} = [ G_k c - H \eta ]
        //           [  0             ]
	//
	// Then we will perform the QR decomposition of Gkplus1
	// to get Qplus1. 
	// 
	// NB: The LAPACK QR decomposition will overwrite 
	// the original G_{k+1} matrix so I will call it 
	// Qkplus1 (for Q_{k+1}) right away.
	multi2d<DComplex> Qkplus1(n_deflate+1,prev_dim+1);
	
	// First copy in the eigenvectors
	for(int col=0; col < n_deflate; ++col) { 
	  for(int row=0; row < prev_dim; ++row) { 
	    Qkplus1(col,row) = evecs( order_array[col], row);
	  }
	  Qkplus1(col,prev_dim) = zero;
	}

	// Q_{k+1} + G_{k+1} =[ G_k |  c - H \eta ]
	//                    [  0  |             ]

	for(int row = 0; row < prev_dim+1; ++row) { 
	  Qkplus1(n_deflate,row) = c_[row];
	  for(int col=0; col < prev_dim; ++col) { 
	    Qkplus1(n_deflate,row) -= H_(col,row)*eta_(col);
	  }
	}

	// QR Decompose Qkplus1
	multi1d<DComplex> tau_kplus1; // QR Decomposition Factors 
	QDPIO::cout << "QR Decomposing Gk+1" << std::endl;
	QDPLapack::zgeqrf(prev_dim+1, n_deflate+1, Qkplus1, tau_kplus1);

	/// Copy H into H_copy
	multi2d<DComplex> H_copy(prev_dim, prev_dim+1);
	for(int col=0; col < prev_dim; ++col) {
	  for(int row=0; row < prev_dim + 1; ++row) { 
	    H_copy(col,row) = H_(col,row);
	  }
	}

	char side='R';
	char trans='N';
	QDPIO::cout << "Post multiplying with Q_k using ZUNMQR 1" << std::endl;
	
	// !!!! NB: The dimensions here are still the original dimensions of H_copy,
	// !!!! even tho at the end of this result only the (rows=prev_dim+1)x(cols=n_deflate) portion is 
	// !!!! valid
	QDPLapack::zunmqr2(side,trans, prev_dim+1, prev_dim, n_deflate, Qkplus1, tau_kplus1, H_copy);

	QDPIO::cout << "Pre multiply with Q^{H}_{k+1} usign ZUNMQR2" << std::endl;
	trans='C'; // Multiply with Herm Conjugate
	side='L';  // Multiply from Left
	// !!!! NB: The dimensions here are still the original dimensions of H_copy2 
	// !!!! Even tho when we are done here, really only the (rows=k+1)x(cols=k) portion of it 
	// !!!! is what is relevant
	QDPLapack::zunmqr2(side,trans, prev_dim+1, prev_dim, n_deflate+1, Qkplus1, tau_kplus1, H_copy);

	// H_copy is the (rows=n_deflate+1, cols=n_deflate) part of 'H' with which we will start the cycle.
	
	// Now I need to form Q_{k+1} explicitly to form the new bases  V_{k+1} and Z_{k}
	QDPLapack::zungqr(prev_dim+1,n_deflate+1,n_deflate+1, Qkplus1, tau_kplus1);


	multi1d<LatticeFermion> new_V(n_deflate+1);
	for(int i=0; i < n_deflate+1; ++i) {
	  new_V[i][s] = zero;
	  for(int j=0; j < prev_dim+1; ++j) { 
	    new_V[i][s] += V_[j]*Qkplus1(i,j);
	  }
	}

	multi1d<LatticeFermion> new_Z(n_deflate);
	for(int i=0; i < n_deflate; ++i) {
	  new_Z[i][s] = zero;
	  for(int j=0; j < prev_dim; ++j) { 
	    new_Z[i][s] += Z_[j]*Qkplus1(i,j);
	  }
	}

	// Reinit things:
	int total_dim = n_krylov+n_deflate;

	for(int col=0; col < total_dim; ++col) { 
	  for(int row=0; row < total_dim; ++row) { 
	    // Nuke H
	    H_(col,row) = zero;
	    R_(col,row) = zero;
	  }
	}

	for(int row=0; row < total_dim; ++row) { 
	  c_[row] = zero;
	  g_[row] = zero;
	}

	// Reinit space -- copy in the first new Z's and zero the rest
	for(int cols=0; cols < n_deflate; ++cols) { 
	  Z_[cols][s] = new_Z[cols];
	}
	for(int cols=n_deflate; cols < total_dim; ++cols) {
	  Z_[cols][s] = zero;
	}

	// Reinit space -- copy in the new V's and zero the rest
	for(int cols=0; cols < n_deflate+1; ++cols) { 
	  V_[cols][s] = new_V[cols];
	}

	for(int cols=n_deflate+1; cols < total_dim+1; ++cols) {
	  V_[cols][s] = zero;
	}

	// c = [ V^H_{k+1} r ] 
	//     [      0      ]
	for(int row=0; row < n_deflate+1; ++row) { 
	  c_[row] = innerProduct( V_[row],r,s );  
	  g_[row] = c_[row];
	}
	for(int row=n_deflate+1; row < total_dim+1; ++row) {
	  c_[row] = zero;
	  g_[row] = zero;
	}
	
	Hk_QR_.resize(n_deflate,n_deflate+1);
	Hk_QR_taus_.resize(n_deflate+1);

	for(int col=0; col < n_deflate; ++col) { 
	  for(int row=0; row < n_deflate+1; ++row) { 
	    H_(col,row) = H_copy(col,row);
	    Hk_QR_(col,row) = H_(col,row);
	  }
	}

	// Now I want to do a QR decomposition of Hk_QR_
	QDPLapack::zgeqrf(n_deflate+1, n_deflate, Hk_QR_, Hk_QR_taus_);
	for(int col=0; col < n_deflate; col++) {
	  for(int row=col; row < n_deflate; row++) { 
	    R_(col,row) = Hk_QR_(col,row);
	  }
	}

	// g = Q_H c  for incremental residuum 
        side = 'L';
	trans = 'C';
	QDPLapack::zunmqrv(side,trans, n_deflate+1, n_deflate+1, Hk_QR_, Hk_QR_taus_, g_);

#if 0
	for(int row=0; row < total_dim+1; ++row) { 
	  QDPIO::cout << " g[" << row << "]=" << g_[row] << std::endl;
	}
#endif
      }

      // Carry out Flexible Arnoldi process for the cycle

      // We are solving for the defect:   A dx = r
      // so the RHS in the Arnoldi process is 'r'
      // NB: We recompute a true 'r' after every cycle
      // So in the cycle we could in principle 
      // use reduced precision... TBInvestigated.
      FlexibleArnoldi(n_krylov,
		      n_deflate,
		      target,
		      V_,
		      Z_, 
		      H_, 
		      R_,
		      givens_rots_,
		      g_,
		      Hk_QR_,
		      Hk_QR_taus_,
		      dim);
      
      int iters_this_cycle = dim - n_deflate;

      QDPIO::cout << "Arnoldi cycle finished: dim=" << dim << std::endl;
      // Solve the least squares system for the lsq_coeffs coefficients
#if 0
      LeastSquaresSolve(H_,c_,eta_, dim); // Solve Least Squares System
#endif
      LeastSquaresSolve(R_,g_,eta_, dim); // Solve Least Squares System
      // Compute the correction dx = sum_j  eta_j Z_j
      LatticeFermion dx = zero;
      for(int j=0; j < dim; ++j) { 
	dx[s] += eta_[j]*Z_[j];
      }

      // Update psi
      psi[s] += dx;


      QDPIO::cout << "FGMRESDR: Cycle finished with " << iters_this_cycle << " iterations" << std::endl;

      // Recompute r
      r[s] = chi;
      (*A_)(tmp, psi, PLUS);
      r[s] -= tmp;  // This 'r' will be used in next cycle as the || rhs ||
      
      // Recompute true norm
      r_norm = sqrt(norm2(r,s));
      QDPIO::cout << "FGMRESDR: || r || = " << r_norm <<  " target = " << target << std::endl;

      // Update total iters
      iters_total += iters_this_cycle;

      // Check if we are done either via convergence, or runnign out of iterations
      finished = toBool( r_norm <= target ) || (iters_total >= invParam_.MaxIter);
      prev_dim = dim;

    }

    // Either we've exceeded max iters, or we have converged in either case set res:
    res.n_count = iters_total;
    res.resid = r_norm;
    QDPIO::cout << "FGMRESDR: Done. Cycles=" << n_cycles << ", Iters=" << iters_total << " || r ||/|| b ||=" << r_norm / norm_rhs << " Target=" << invParam_.RsdTarget << std::endl;
    END_CODE();
    return res;

  }




    

}


