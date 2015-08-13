/*! \file
 *  \brief Solve a M*psi=chi linear system by MR
 */

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
   *
   */
  template< typename T >
  void FlexibleArnoldiT(int n_krylov, 
			int n_deflate,                  // Size of augmented space. Not yet used.
			const Real& rsd_target,
			const Double& r_norm,           // Caller computes this, so save me doing it.
			const LinearOperator<T>& A,     // Operator
			const LinOpSystemSolver<T>& M,  // Preconditioner
			const T& rhs,
			multi1d<T>& V,
			multi1d<T>& Z,
			multi2d<DComplex>& H,
			multi2d<DComplex>& R,
			multi1d< Handle<Givens> >& givens_rots,
			multi1d<DComplex>& g,
			int& ndim_cycle)

  {
    const Subset& s = M.subset();      // Linear Operator Subset
    ndim_cycle = 0;     


    // Set up initial vector g = [ beta, 0 ... 0 ]^T
    // In this case beta should be the r_norm, ie || r || (=|| b || for the first cycle)
    for(int j=0; j < g.size(); ++j) { 
      g[j] = DComplex(0); 
    }
    g[0] = r_norm;
    
    // Set up initial V[0] = rhs / || r^2 || -- 
    Double beta_inv = Double(1)/r_norm;
    V[0][s] = beta_inv * rhs;
    
    // NB: Because the internal loop goes < j
    // I need 'j' to start from 1 - a la Fortran
    // However when used as an index I need it ala C so j-1
    // In this game, j is a column
    // and i is a row.

    // Work by columns:
    for(int j=0; j < n_krylov; ++j) {
     
      M( Z[j], V[j] );  // z_j = M v_j
      T w;
      A( w, Z[j], PLUS);  // w  = A z_

      // Fill out column j
      for(int i=0; i <= j ;  ++i ) {
	H(i,j) = innerProduct(V[i], w, s);
	w[s] -= H(i,j)* V[i];
      }
      
      Double wnorm=sqrt(norm2(w,s));
      H(j+1,j) = DComplex(wnorm);

      // In principle I should check w_norm to be 0, and if it is I should
      // terminate: Code Smell -- get rid of 1.0e-14.
      if ( toBool( fabs( wnorm ) < Double(1.0e-14) ) ) {

	// If wnorm = 0 exactly, then we have converged exactly
	// Replay Givens rots here, how to test?

	QDPIO::cout << "Converged at iter = " << j << std::endl;
	ndim_cycle = j;
	return;
      }

      Double invwnorm = Double(1)/wnorm;
      V[j+1] = invwnorm*w;


      // Done construcing H. Construct relevant column of R
      for(int i=0; i <= j+1; ++i) { 
	R(i,j) = H(i,j);
      }

      // Apply Existing Givens Rotations to this column of R
      for(int i=0;i < j; ++i) {
	(*givens_rots[i])(j,R);
      }
      givens_rots[j] = new Givens(j,R); // Compute next Givens Rot
      (*givens_rots[j])(j,R); // Apply it to R
      (*givens_rots[j])(g);   // Apply it to the g vector

      Double accum_resid = fabs(real(g[j+1]));
      QDPIO::cout << "Iter " << j << " || r || = " << accum_resid << " Target=" <<rsd_target << std::endl;

      ndim_cycle = j+1;
      if ( toBool( accum_resid <= rsd_target ) ) { 
	QDPIO::cout << "Converged at iter = " << j << std::endl;
	return;
      }
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
						 const multi1d<DComplex>& rhs,
						 multi1d<DComplex>& eta,
						 int n_cols) const
  {

    /* Assume here we have a square matrix with an extra row.
       Hence the loop counters are the columns not the rows.
       NB: For an augmented system this will change */
    eta[n_cols-1] = rhs[n_cols-1]/R(n_cols-1,n_cols-1);
    for(int row = n_cols-2; row >= 0; --row) {
      eta[row] = rhs[row];
      for(int col=row+1; col <  n_cols; ++col) { 
	eta[row] -= R(row,col)*eta[col];
      }
      eta[row] /= R(row,row);
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
					  const Double& r_norm,
					  const T& rhs, 
					  multi1d<T>& V,
					  multi1d<T>& Z, 
					  multi2d<DComplex>&H,
					  multi2d<DComplex>&R,
					  multi1d<Handle<Givens> >& givens_rots,
					  multi1d<DComplex> &g,
					  int& ndim_cycle) const
  {
    FlexibleArnoldiT<>(n_krylov, 
		       n_deflate,
		       rsd_target,
		       r_norm,
		       (*A_),
		       (*preconditioner_),
		       rhs, 
		       V,Z,H,R,givens_rots,g, ndim_cycle);
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

    // We are done if norm is sufficiently accurate, 
    bool finished = toBool( r_norm <= target ) ;

    // We keep executing cycles until we are finished
    while( !finished ) { 

      // If not finished, we should do another cycle with RHS='r' to find dx to add to psi.
      ++n_cycles;
      
      // Declare workspace
      // Work spacesn_
      multi2d<DComplex> H(invParam_.NKrylov+1, invParam_.NKrylov); // The H matrix
      multi2d<DComplex> R(invParam_.NKrylov+1, invParam_.NKrylov); // R = H diagonalized with Givens rotations
      multi1d<T> V(invParam_.NKrylov+1);  // K(A)
      multi1d<T> Z(invParam_.NKrylov+1);  // K(MA)
      multi1d< Handle<Givens> > givens_rots(invParam_.NKrylov+1);
      multi1d<DComplex> g(invParam_.NKrylov+1);

      // Init the matrices
      for(int row = 0; row < invParam_.NKrylov+1; ++row) {
	for(int col = 0; col < invParam_.NKrylov; ++col) { 
	  H(row,col) = DComplex(0);
	  R(row,col) = DComplex(0);
	}
      }
  
      int dim; // dim at the end of cycle (in case we terminate in-cycle


      // Carry out Flexible Arnoldi process for the cycle

      // We are solving for the defect:   A dx = r
      // so the RHS in the Arnoldi process is 'r'
      // NB: We recompute a true 'r' after every cycle
      // So in the cycle we could in principle 
      // use reduced precision... TBInvestigated.
      FlexibleArnoldi(invParam_.NKrylov, 
		      invParam_.NDefl,
		      target,
		      r_norm,
		      r, 
		      V,
		      Z, 
		      H, 
		      R,
		      givens_rots,
		      g,
		      dim);
      

      // Solve the least squares system for the lsq_coeffs coefficients
      multi1d<DComplex> lsq_coeffs(dim);  
      LeastSquaresSolve(R,g,lsq_coeffs, dim); // Solve Least Squares System

      // Compute the correction dx = sum_j  lsq_coeffs_j Z_j
      LatticeFermion dx = zero;
      for(int j=0; j < dim; ++j) { 
	dx[s] += lsq_coeffs[j]*Z[j];
      }

      // Update psi
      psi[s] += dx;


      QDPIO::cout << "FGMRESDR: Cycle finished with " << dim << " iterations" << std::endl;

      // Recompute r
      r[s] = chi;
      (*A_)(tmp, psi, PLUS);
      r[s] -= tmp;  // This 'r' will be used in next cycle as the || rhs ||
      
      // Recompute true norm
      r_norm = sqrt(norm2(r,s));
      QDPIO::cout << "FGMRESDR: || r || = " << r_norm <<  " target = " << target << std::endl;

      // Update total iters
      iters_total += dim;

      // Check if we are done either via convergence, or runnign out of iterations
      finished = toBool( r_norm <= target ) || (iters_total >= invParam_.MaxIter);
    }

    // Either we've exceeded max iters, or we have converged in either case set res:
    res.n_count = iters_total;
    res.resid = r_norm;
    QDPIO::cout << "FGMRESDR: Done. Cycles=" << n_cycles << ", Iters=" << iters_total << " || r ||/|| b ||=" << r_norm / norm_rhs << " Target=" << invParam_.RsdTarget << std::endl;
    END_CODE();
    return res;

  }




    

}


