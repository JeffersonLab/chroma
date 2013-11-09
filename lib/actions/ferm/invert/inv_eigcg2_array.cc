// -*- C++ -*-
// $Id: inv_eigcg2_array.cc,v 1.3 2008-01-16 19:03:57 edwards Exp $
/*! \file
 *  \brief Conjugate-Gradient algorithm with eigenvector acceleration
 */

#include <qdp-lapack.h>
//#include "octave_debug.h"
//#include "octave_debug.cc"

#include "actions/ferm/invert/inv_eigcg2_array.h"

#include "actions/ferm/invert/containers.h"
#include "actions/ferm/invert/norm_gram_schm.h"

//#define DEBUG
#define DEBUG_FINAL


// NEEDS A LOT OF CLEAN UP
namespace Chroma 
{
  using namespace LinAlg;
  //using namespace Octave;

  namespace InvEigCG2ArrayEnv
  {
    // needs a little clean up... for possible exceptions if the size of H
    // is smaller than hind 
    // LatticeFermionF
    template<typename T>
    void SubSpaceMatrix_T(LinAlg::Matrix<DComplex>& H,
			  const LinearOperatorArray<T>& A,
			  const multi2d<T>& evec,
			  int Nvecs)
    {
      H.N = Nvecs;
      if(Nvecs>evec.size2()){
	H.N = evec.size2();
      }
      if(H.N>H.size()){
	QDPIO::cerr<<"OOPS! your matrix can't take this many matrix elements\n";
	exit(1);	
      }
      // fill H  with zeros
      H.mat = 0.0;
    
      multi1d<T> Ap;
      for(int i(0);i<H.N;i++)
      {
	A(Ap,evec[i],PLUS);
	for(int j(i);j<H.N;j++)
	{
	  H(j,i) = innerProduct(evec[j], Ap, A.subset());
	  //enforce hermiticity
	  H(i,j) = conj(H(j,i));
	  if(i==j) H(i,j) = real(H(i,j));
	}
      } 
    }


    template<typename T>
    SystemSolverResults_t InvEigCG2_T(const LinearOperatorArray<T>& A,
				      multi1d<T>& x, 
				      const multi1d<T>& b,
				      multi1d<Double>& eval, 
				      multi2d<T>& evec,
				      int Neig,
				      int Nmax,
				      const Real& RsdCG, int MaxCG)
    {
      START_CODE();

      FlopCounter flopcount;
      flopcount.reset();
      StopWatch swatch;
      swatch.reset();
      swatch.start();
    
      const int Ls = A.size();

      SystemSolverResults_t  res;
    
      multi1d<T> p; 
      multi1d<T> Ap; 
      multi1d<T> r;
      //multi1d<T> z;

      Double rsd_sq = (RsdCG * RsdCG) * Real(norm2(b,A.subset()));
      Double alphaprev, alpha,pAp;     
      Real beta; 
      Double betaprev ;
      Double r_dot_z, r_dot_z_old;

      int k = 0;
      A(Ap,x,PLUS);
      for(int l=0; l < Ls; ++l)
	r[l][A.subset()] = b[l] - Ap[l];


#if 1
      QDPIO::cout << "InvEigCG2: Nevecs(input) = " << evec.size2() << endl;
      QDPIO::cout << "InvEigCG2: k = " << k << "  res^2 = " << r_dot_z << endl;
#endif
      Matrix<DComplex> H(Nmax); // square matrix containing the tridiagonal
      VectorArrays<T> vec(Nmax,Ls); // contains the vectors we use...
      

      beta=0.0;
      alpha = 1.0;
      multi1d<T> Ap_prev;
      multi1d<T> tt;

      // Algorithm from page 529 of Golub and Van Loan
      // Modified to match the m-code from A. Stathopoulos
      while(toBool(r_dot_z>rsd_sq))
      {
	/** preconditioning algorithm **/
	//z[A.subset()]=r; //preconditioning can be added here
	//r_dot_z = innerProductReal(r,z,A.subset());
	//****//
	r_dot_z_old = r_dot_z;
	r_dot_z = norm2(r,A.subset());
	Double inv_sqrt_r_dot_z = 1.0/sqrt(r_dot_z);
	k++;
	if(k==1){
	  //p[A.subset()] = z;	
	  for(int l=0; l < Ls; ++l)
	    p[l][A.subset()] = r[l];	
	}
	else{
	  betaprev = beta;
	  beta = r_dot_z/r_dot_z_old;
	  //p[A.subset()] = z + beta*p; 
	  for(int l=0; l < Ls; ++l)
	    p[l][A.subset()] = r[l] + beta*p[l];
	}
	//-------- Eigenvalue eigenvector finding code ------
	if((Neig>0)&& (H.N == Nmax)) 
	  for(int l=0; l < Ls; ++l)
	    Ap_prev[l][A.subset()] = Ap[l];
	//---------------------------------------------------
	A(Ap,p,PLUS);
	
	//-------- Eigenvalue eigenvector finding code ------
	if(Neig>0)
	{
	  if(k>1)
	    H(H.N-1,H.N-1) = 1/alpha + betaprev/alphaprev;
	  if(vec.N==Nmax)
	  {
	    QDPIO::cout<<"MAGIC BEGINS: H.N ="<<H.N<<endl;
#ifdef DEBUG
	    {
	      stringstream tag;
	      tag<<"H"<<k;
	      OctavePrintOut(H.mat,Nmax,tag.str(),"Hmatrix.m");
	    }
	    {
	      Matrix<DComplex> tmp(Nmax); 
	      SubSpaceMatrix(tmp,A,vec.vec,vec.N);
	      stringstream tag;
	      tag<<"H"<<k<<"ex";
	      OctavePrintOut(tmp.mat,Nmax,tag.str(),"Hmatrix.m");
	    }
#endif
	    multi2d<DComplex> Hevecs(H.mat);
	    multi1d<Double> Heval;
	    char V = 'V'; char U = 'U';
	    QDPLapack::zheev(V,U,Nmax,Hevecs,Heval);
#ifdef DEBUG
	    {
	      stringstream tag;
	      tag<<"Hevecs"<<k;
	      OctavePrintOut(Hevecs,Nmax,tag.str(),"Hmatrix.m");
	    }
	    for(int i(0);i<Nmax;i++)
	      QDPIO::cout<<" eignvalue: "<<Heval[i]<<endl;
#endif
	    multi2d<DComplex> Hevecs_old(H.mat);
	    multi1d<Double> Heval_old;
	    QDPLapack::zheev(V,U,Nmax-1,Hevecs_old,Heval_old);
	    for(int i(0);i<Nmax;i++)	    
	      Hevecs_old(i,Nmax-1) = Hevecs_old(Nmax-1,i) = 0.0;
	    
	    //Reduce to 2*Neig vectors
	    H.N = Neig + Neig; // Thickness of restart 2*Neig

	    for(int i(Neig);i<2*Neig;i++) 
	      for(int j(0);j<Nmax;j++)	    
		Hevecs(i,j) = Hevecs_old(i-Neig,j);

	    // Orthogonalize the last Neig vecs (keep the first Neig fixed)
	    // zgeqrf(Nmax, 2*Neig, Hevecs, Nmax,
	    //  TAU_CMPLX_ofsize_2Neig, Workarr, 2*Neig*Lapackblocksize, info);
	    multi1d<DComplex> TAU;
	    QDPLapack::zgeqrf(Nmax,2*Neig,Hevecs,TAU);
	    char R = 'R'; char L = 'L'; char N ='N'; char C = 'C';
	    multi2d<DComplex> Htmp = H.mat;
	    QDPLapack::zunmqr(R,N,Nmax,Nmax,Hevecs,TAU,Htmp);
	    QDPLapack::zunmqr(L,C,Nmax,2*Neig,Hevecs,TAU,Htmp);

	    QDPLapack::zheev(V,U,2*Neig,Htmp,Heval);
#ifdef DEBUG
	    {
	      stringstream tag;
	      tag<<"Htmp"<<k;
	      OctavePrintOut(Htmp,Nmax,tag.str(),"Hmatrix.m");
	    }
#endif
	    for(int i(Neig); i< 2*Neig;i++ ) // mhpws prepei na einai 0..2*Neig
	      for(int j(2*Neig); j<Nmax; j++) 
		Htmp(i,j) =0.0;

#ifdef DEBUG
	    {
	      stringstream tag;
	      tag<<"HtmpBeforeZUM"<<k;
	      OctavePrintOut(Htmp,Nmax,tag.str(),"Hmatrix.m");
	    }
#endif
	    QDPLapack::zunmqr(L,N,Nmax,2*Neig,Hevecs,TAU,Htmp);
#ifdef DEBUG
	    {
	      stringstream tag;
	      tag<<"HtmpAfeterZUM"<<k;
	      OctavePrintOut(Htmp,Nmax,tag.str(),"Hmatrix.m");
	    }
#endif
	    
#ifndef USE_BLAS_FOR_LATTICEFERMIONS
	    multi2d<T> tt_vec(2*Neig,Ls);
	    for(int i(0);i<2*Neig;i++){
	      for(int l=0; l < Ls; ++l)
		tt_vec[i][l][A.subset()] = zero;
	      for(int j(0);j<Nmax;j++)
		for(int l=0; l < Ls; ++l)
		  tt_vec[i][l][A.subset()] +=Htmp(i,j)*vec[j][l]; 
	    }
	    for(int i(0);i<2*Neig;i++)
	      for(int l=0; l < Ls; ++l)
		vec[i][l][A.subset()] = tt_vec[i][l];
#else
	    // Here I need the restart_X bit
	    // zgemm("N", "N", Ns*Nc*Vol/2, 2*Neig, Nmax, 1.0,
	    //      vec.vec, Ns*Nc*Vol, Htmp, Nmax+1, 0.0, tt_vec, Ns*Nc*Vol);
	    // copy apo tt_vec se vec.vec

#endif
	    vec.N = 2*Neig; // restart the vectors to keep

	    H.mat = 0.0; // zero out H 
	    for (int i=0;i<2*Neig;i++) H(i,i) = Heval[i];

	    //A(tt,r,PLUS);
	    for(int l=0; l < Ls; ++l)
	      tt[l] = Ap[l] - beta*Ap_prev[l]; //avoid the matvec

#ifndef USE_BLAS_FOR_LATTICEFERMIONS
	    for (int i=0;i<2*Neig;i++){
	      H(2*Neig,i)=innerProduct(vec[i],tt,A.subset())*inv_sqrt_r_dot_z;
	      H(i,2*Neig)=conj(H(2*Neig,i));
	      //H(i,2*Neig)=innerProduct(vec[i],tt,A.subset())*inv_sqrt_r_dot_z;
	      //H(2*Neig,i)=conj(H(i,2*Neig));
	    }
#else  
	    // Optimized code for creating the H matrix row and column
	    // asume even-odd layout: Can I check if this is the case at 
	    // compile time? or even at run time? This is a good test to 
	    // have to avoid wrong results.
	    
	    
#endif
	  }//H.N==Nmax
	  else{
	    if(k>1)
	      H(H.N-1,H.N) = H(H.N,H.N-1) = -sqrt(beta)/alpha;
	  }
	  H.N++;
	  //vec.NormalizeAndAddVector(z,inv_sqrt_r_dot_z,A.subset());
	  vec.NormalizeAndAddVector(r,inv_sqrt_r_dot_z,A.subset());
	}

	pAp = innerProductReal(p,Ap,A.subset());
      	alphaprev = alpha;// additional line for Eigenvalue eigenvector code
	alpha = r_dot_z/pAp;
	for(int l=0; l < Ls; ++l)
	  x[l][A.subset()] += alpha*p[l];
	for(int l=0; l < Ls; ++l)
	  r[l][A.subset()] -= alpha*Ap[l];
	      
	
	//---------------------------------------------------
	if(k>MaxCG){
	  res.n_count = k;
	  res.resid   = sqrt(r_dot_z);
	  QDP_error_exit("too many CG iterations: count = %d", res.n_count);
	  END_CODE();
	  return res;
	}
#if 1
	QDPIO::cout << "InvEigCG2: k = " << k ;
	QDPIO::cout << "  r_dot_z = " << r_dot_z << endl;
#endif
      }//end CG loop

      if(Neig>0)
      {
	evec.resize(Neig,Ls);
	eval.resize(Neig);
	
#define  USE_LAST_VECTORS
#ifdef USE_LAST_VECTORS
	
	multi2d<DComplex> Hevecs(H.mat);
	multi1d<Double> Heval;
	char V = 'V'; char U = 'U';
	QDPLapack::zheev(V,U,H.N-1,Hevecs,Heval);
	
	for(int i(0);i<Neig;i++){
	  for(int l=0; l < Ls; ++l)
	    evec[i][l][A.subset()] = zero;
	  eval[i] = Heval[i];
	  for(int j(0);j<H.N-1;j++)
	    for(int l=0; l < Ls; ++l)
	      evec[i][l][A.subset()] +=Hevecs(i,j)*vec[j][l];
	}
#else
	for(int i(0);i<Neig;i++){
	  evec[i][A.subset()] = vec[i];
	  eval[i] = real(H(i,i));
	}
#endif
      }
      res.n_count = k;
      res.resid = sqrt(r_dot_z);
      swatch.stop();
      QDPIO::cout << "InvEigCG2: k = " << k << endl;
      flopcount.report("InvEigCG2", swatch.getTimeInSeconds());
      END_CODE();
      return res;
    }

   
    template<typename T>
    SystemSolverResults_t old_InvEigCG2_T(const LinearOperatorArray<T>& A,
					  multi1d<T>& x, 
					  const multi1d<T>& b,
					  multi1d<Double>& eval, 
					  multi2d<T>& evec,
					  int Neig,
					  int Nmax,
					  const Real& RsdCG, int MaxCG)
    {
      START_CODE();

      FlopCounter flopcount;
      flopcount.reset();
      StopWatch swatch;
      swatch.reset();
      swatch.start();

      const int Ls = A.size();
    
      SystemSolverResults_t  res;
    
      multi1d<T> p; 
      multi1d<T> Ap; 
      multi1d<T> r,z;

      Double rsd_sq = (RsdCG * RsdCG) * Real(norm2(b,A.subset()));
      Double alphaprev, alpha,pAp;
      Real beta;
      Double r_dot_z, r_dot_z_old;
      //Complex r_dot_z, r_dot_z_old,beta;
      //Complex alpha,pAp;

      int k = 0;
      A(Ap,x,PLUS);
      r[A.subset()] = b - Ap;
      Double r_norm2 = norm2(r,A.subset());

#if 1
      QDPIO::cout << "InvEigCG2: Nevecs(input) = " << evec.size2() << endl;
      QDPIO::cout << "InvEigCG2: k = " << k << "  res^2 = " << r_norm2 << endl;
#endif
      Real inorm(Real(1.0/sqrt(r_norm2)));
      bool FindEvals = (Neig>0);
      int tr; // Don't know what this does...
      bool from_restart;
      Matrix<DComplex> H(Nmax); // square matrix containing the tridiagonal
      VectorArrays<T> vec(Nmax,Ls); // contains the vectors we use...
      //-------- Eigenvalue eigenvector finding code ------
      vec.AddVectors(evec,A.subset());
      //QDPIO::cout<<"vec.N="<<vec.N<<endl;
      p[A.subset()] = inorm*r; // this is not needed GramSchmidt will take care of it
      vec.AddOrReplaceVector(p,A.subset());
      //QDPIO::cout<<"vec.N="<<vec.N<<endl;
      normGramSchmidt(vec.vec,vec.N-1,vec.N,A.subset());
      normGramSchmidt(vec.vec,vec.N-1,vec.N,A.subset()); //call twice: need to improve...
      SubSpaceMatrix(H,A,vec.vec,vec.N);
      from_restart = true;
      tr = H.N - 1; // this is just a flag as far as I can tell  
      //-------

      Double betaprev;
      beta=0.0;
      alpha = 1.0;
      // Algorithm from page 529 of Golub and Van Loan
      // Modified to match the m-code from A. Stathopoulos
      while(toBool(r_norm2>rsd_sq))
      {
	/** preconditioning algorithm **/
	z[A.subset()]=r; //preconditioning can be added here
	/**/
	r_dot_z = innerProductReal(r,z,A.subset());
	k++;
	betaprev = beta; // additional line for Eigenvalue eigenvector code
	if(k==1){
	  p[A.subset()] = z;	
	  H.N++;
	}
	else{
	  beta = r_dot_z/r_dot_z_old;
	  p[A.subset()] = z + beta*p; 
	  //-------- Eigenvalue eigenvector finding code ------
	  // fist block
	  if(FindEvals){
	    if(!((from_restart)&&(H.N == tr+1))){
	      if(!from_restart){
		H(H.N-2,H.N-2) = 1/alpha + betaprev/alphaprev;
	      } 
	      from_restart = false; 
	      H(H.N-1,H.N-2) = -sqrt(beta)/alpha;
	      H(H.N-2,H.N-1) = -sqrt(beta)/alpha;
	    }
	    H.N++;
	  }
	  //---------------------------------------------------
	}
	A(Ap,p,PLUS);
	pAp = innerProductReal(p,Ap,A.subset());
      
	alphaprev = alpha;// additional line for Eigenvalue eigenvector code
	alpha = r_dot_z/pAp;
	x[A.subset()] += alpha*p;
	r[A.subset()] -= alpha*Ap;
	r_norm2 =  norm2(r,A.subset());
	r_dot_z_old = r_dot_z;

      
	//-------- Eigenvalue eigenvector finding code ------
	// second block
	if(FindEvals){
	  if (vec.N==Nmax){//we already have stored the maximum number of vectors
	    // The magic begins here....
	    QDPIO::cout<<"MAGIC BEGINS: H.N ="<<H.N<<endl;
	    H(Nmax-1,Nmax-1) = 1/alpha + beta/alphaprev;

#ifdef DEBUG
	    {
	      stringstream tag;
	      tag<<"H"<<k;
	      OctavePrintOut(H.mat,Nmax,tag.str(),"Hmatrix.m");
	    }
	    {
	      Matrix<DComplex> tmp(Nmax); 
	      SubSpaceMatrix(tmp,A,vec.vec,vec.N);
	      stringstream tag;
	      tag<<"H"<<k<<"ex";
	      OctavePrintOut(tmp.mat,Nmax,tag.str(),"Hmatrix.m");
	    }
#endif
	    //exit(1);
	    multi2d<DComplex> Hevecs(H.mat);
	    multi1d<Double> Heval;
	    char V = 'V'; char U = 'U';
	    QDPLapack::zheev(V,U,Hevecs,Heval);
	    multi2d<DComplex> Hevecs_old(H.mat);

#ifdef DEBUG
	    {
	      multi1d<T> tt_vec(vec.size());
	      for(int i(0);i<Nmax;i++){
		tt_vec[i][A.subset()] = zero;
		for(int j(0);j<Nmax;j++)
		  tt_vec[i][A.subset()] +=Hevecs(i,j)*vec[j]; // NEED TO CHECK THE INDEXINT
	      }
	      // CHECK IF vec are eigenvectors...
	  
	      multi1d<T> Av;
	      for(int i(0);i<Nmax;i++){
		A(Av,tt_vec[i],PLUS);
		DComplex rq = innerProduct(tt_vec[i],Av,A.subset());
		Av[A.subset()] -= Heval[i]*tt_vec[i];
		Double tt = sqrt(norm2(Av,A.subset()));
		QDPIO::cout<<"1 error eigenvector["<<i<<"] = "<<tt<<" ";
		tt =  sqrt(norm2(tt_vec[i],A.subset()));
		QDPIO::cout<<" --- eval = "<<Heval[i]<<" ";
		QDPIO::cout<<" --- rq = "<<real(rq)<<" ";
		QDPIO::cout<<"--- norm = "<<tt<<endl ;
	      }
	    }
#endif
	    multi1d<Double> Heval_old;
	    QDPLapack::zheev(V,U,Nmax-1,Hevecs_old,Heval_old);
	    for(int i(0);i<Nmax;i++)	    
	      Hevecs_old(i,Nmax-1) = Hevecs_old(Nmax-1,i) = 0.0;
#ifdef DEBUG
	    {
	      stringstream tag;
	      tag<<"Hevecs_old"<<k;
	      OctavePrintOut(Hevecs_old,Nmax,tag.str(),"Hmatrix.m");
	    }
	  
#endif
	    tr = Neig + Neig; // v_old = Neig optimal choice
	   
	    for(int i(Neig);i<tr;i++)
	      for(int j(0);j<Nmax;j++)	    
		Hevecs(i,j) = Hevecs_old(i-Neig,j);
#ifdef DEBUG
	    {
	      stringstream tag;
	      tag<<"Hevecs"<<k;
	      OctavePrintOut(Hevecs,Nmax,tag.str(),"Hmatrix.m");
	    }
#endif

	    // Orthogonalize the last Neig vecs (keep the first Neig fixed)
	    // zgeqrf(Nmax, 2*Neig, Hevecs, Nmax,
	    //    TAU_CMPLX_ofsize_2Neig, Workarr, 2*Neig*Lapackblocksize, info);
	    multi1d<DComplex> TAU;
	    QDPLapack::zgeqrf(Nmax,2*Neig,Hevecs,TAU);
	    char R = 'R'; char L = 'L'; char N ='N'; char C = 'C';
	    multi2d<DComplex> Htmp = H.mat;
	    QDPLapack::zunmqr(R,N,Nmax,Nmax,Hevecs,TAU,Htmp);
	    QDPLapack::zunmqr(L,C,Nmax,2*Neig,Hevecs,TAU,Htmp);
	    // Notice that now H is of size 2Neig x 2Neig, 
	    // but still with LDA = Nmax 
#ifdef DEBUG
	    {
	      stringstream tag;
	      tag<<"Htmp"<<k;
	      OctavePrintOut(Htmp,Nmax,tag.str(),"Hmatrix.m");
	    }
#endif
	    QDPLapack::zheev(V,U,2*Neig,Htmp,Heval);
	    for(int i(Neig); i< 2*Neig;i++ ) 
	      for(int j(2*Neig); j<Nmax; j++)
		Htmp(i,j) =0.0;
#ifdef DEBUG
	    {
	      stringstream tag;
	      tag<<"evecstmp"<<k;
	      OctavePrintOut(Htmp,Nmax,tag.str(),"Hmatrix.m");
	    }
#endif
	    QDPLapack::zunmqr(L,N,Nmax,2*Neig,Hevecs,TAU,Htmp);
#ifdef DEBUG
	    {
	      stringstream tag;
	      tag<<"Yevecstmp"<<k;
	      OctavePrintOut(Htmp,Nmax,tag.str(),"Hmatrix.m");
	    }
#endif
	    // Here I need the restart_X bit
	    multi1d<T> tt_vec = vec.vec;
	    for(int i(0);i<2*Neig;i++){
	      vec[i][A.subset()] = zero;
	      for(int j(0);j<Nmax;j++)
		vec[i][A.subset()] +=Htmp(i,j)*tt_vec[j]; // NEED TO CHECK THE INDEXING
	    }

#ifdef DEBUG
	    // CHECK IF vec are eigenvectors...
	    {
	      multi1d<T> Av;
	      for(int i(0);i<2*Neig;i++){
		A(Av,vec[i],PLUS);
		DComplex rq = innerProduct(vec[i],Av,A.subset());
		Av[A.subset()] -= Heval[i]*vec[i];
		Double tt = sqrt(norm2(Av,A.subset()));
		QDPIO::cout<<"error eigenvector["<<i<<"] = "<<tt<<" ";
		tt =  sqrt(norm2(vec[i],A.subset()));
		QDPIO::cout<<"--- rq ="<<real(rq)<<" ";
		QDPIO::cout<<"--- norm = "<<tt<<endl ;
	      }

	    }
#endif
	    H.mat = 0.0; // zero out H 
	    for (int i=0;i<2*Neig;i++) H(i,i) = Heval[i];
	  
	    multi1d<T> Ar;
	    // Need to reorganize so that we aboid this extra matvec
	    A(Ar,r,PLUS);
	    Ar /= sqrt(r_norm2); 
	    // this has the oposite convension than the subspace matrix
	    // this is the reason the vector reconstruction in here does not
	    // need the conj(H) while in the Rayleigh Ritz refinement 
	    // it does need it. 
	    // In here the only complex matrix elements are these computined
	    // in the next few lines. This is why the inconsistency 
	    //  does not matter anywhere else. 
	    // In the refinement step though all matrix elements are complex
	    // hence things break down.
	    // Need to fix the convension so that we do not
	    // have these inconsistencies.
	    for (int i=0;i<2*Neig;i++){
	      H(2*Neig,i) = innerProduct(vec[i], Ar, A.subset());
	      H(i,2*Neig) = conj(H(2*Neig,i));
	    }
	    H(2*Neig,2*Neig) = innerProduct(r, Ar, A.subset())/sqrt(r_norm2);
	  
	    H.N = 2*Neig + 1 ; // why this ?
	    from_restart = true;
	    vec.N = 2*Neig; // only keep the lowest Neig. Is this correct???
#ifdef DEBUG
	    {
	      stringstream tag;
	      tag<<"finalH"<<k;
	      OctavePrintOut(H.mat,Nmax,tag.str(),"Hmatrix.m");
	    }
#endif
	  }
	  // Here we add a vector in the list
	  Double inorm = 1.0/sqrt(r_norm2);
	  vec.NormalizeAndAddVector(r,inorm,A.subset());
	  // Shouldn't be the z vector when a preconditioner is used?
	}
	//---------------------------------------------------

	if(k>MaxCG){
	  res.n_count = k;
	  res.resid   = sqrt(r_norm2);
	  QDP_error_exit("too many CG iterations: count = %d", res.n_count);
	  END_CODE();
	  return res;
	}
#if 1
	QDPIO::cout << "InvEigCG2: k = " << k << "  res^2 = " << r_norm2;
	QDPIO::cout << "  r_dot_z = " << r_dot_z << endl;
#endif
      }
      res.n_count = k;
      res.resid = sqrt(r_norm2);
      if(FindEvals)
      {
	// Evec Code ------ Before we return --------
	// vs is the current number of vectors stored
	// Neig is the number of eigenvector we want to compute
	normGramSchmidt(vec.vec,0,2*Neig,A.subset());
	normGramSchmidt(vec.vec,0,2*Neig,A.subset());
	Matrix<DComplex> Htmp(2*Neig);
	SubSpaceMatrix(Htmp,A,vec.vec,2*Neig);

	char V = 'V'; char U = 'U';
	multi1d<Double> tt_eval;
	QDPLapack::zheev(V,U,Htmp.mat,tt_eval);
	evec.resize(Neig,Ls);
	eval.resize(Neig);
	for(int i(0);i<Neig;i++){
	  evec[i][A.subset()] = zero;
	  eval[i] = tt_eval[i];
	  for(int j(0);j<2*Neig;j++)
	    evec[i][A.subset()] += Htmp(i,j)*vec[j];
	}
      
	// I will do the checking of eigenvector quality outside this routine
	// -------------------------------------------
#ifdef DEBUG_FINAL
	// CHECK IF vec are eigenvectors...                                   
	{
	  multi1d<T> Av;
	  for(int i(0);i<Neig;i++)
	  {
	    A(Av,evec[i],PLUS);
	    DComplex rq = innerProduct(evec[i],Av,A.subset());
	    Av[A.subset()] -= eval[i]*evec[i];
	    Double tt = sqrt(norm2(Av,A.subset()));
	    QDPIO::cout<<"FINAL: error eigenvector["<<i<<"] = "<<tt<<" ";
	    tt =  sqrt(norm2(evec[i],A.subset()));
	    QDPIO::cout<<"--- rq ="<<real(rq)<<" ";
	    QDPIO::cout<<"--- norm = "<<tt<<endl ;
	  }
	
	}
#endif
      }

      res.n_count = k;
      res.resid   = sqrt(r_norm2);
      swatch.stop();
      QDPIO::cout << "InvEigCG2: k = " << k << endl;
      flopcount.report("InvEigCG2", swatch.getTimeInSeconds());
      END_CODE();
      return res;
    }

    // A  should be Hermitian positive definite
    template<typename T>
    SystemSolverResults_t vecPrecondCG_T(const LinearOperatorArray<T>& A, 
					 multi1d<T>& x, 
					 const multi1d<T>& b, 
					 const multi1d<Double>& eval, 
					 const multi2d<T>& evec, 
					 int startV, int endV,
					 const Real& RsdCG, int MaxCG)
    {
      START_CODE();

      FlopCounter flopcount;
      flopcount.reset();
      StopWatch swatch;
      swatch.reset();
      swatch.start();

      const int Ls = A.size();

      if(endV>eval.size()){
	QDP_error_exit("vPrecondCG: not enought evecs eval.size()=%d",eval.size());
      }
 
      SystemSolverResults_t  res;
    
      multi1d<T> p; 
      multi1d<T> Ap; 
      multi1d<T> r,z;

      Double rsd_sq = (RsdCG * RsdCG) * Real(norm2(b,A.subset()));
      Double alpha,pAp;
      Real beta;
      Double r_dot_z, r_dot_z_old;
      //Complex r_dot_z, r_dot_z_old,beta;
      //Complex alpha,pAp;

      int k = 0;
      A(Ap,x,PLUS);
      for(int l=0; l < Ls; ++l)
	r[l][A.subset()] = b[l] - Ap[l];
      Double r_norm2 = norm2(r,A.subset());

#if 1
      QDPIO::cout << "vecPrecondCG: k = " << k << "  res^2 = " << r_norm2 << endl;
#endif

      // Algorithm from page 529 of Golub and Van Loan
      while(toBool(r_norm2>rsd_sq)){
	/** preconditioning algorithm **/
	for(int l=0; l < Ls; ++l)
	  z[l][A.subset()] = r[l]; // not optimal but do it for now...
	/**/
	for(int i(startV);i<endV;i++){
	  DComplex d = innerProduct(evec[i],r,A.subset());
	  //QDPIO::cout<<"vecPrecondCG: "<< d<<" "<<(1.0/eval[i]-1.0)<<endl;
	  for(int l=0; l < Ls; ++l)
	    z[l][A.subset()] += (1.0/eval[i]-1.0)*d*evec[i][l];
	}
	/**/
	/**/
	r_dot_z = innerProductReal(r,z,A.subset());
	k++;
	if(k==1){
	  for(int l=0; l < Ls; ++l)
	    p[l][A.subset()] = z[l];
	}
	else{
	  beta = r_dot_z/r_dot_z_old;
	  for(int l=0; l < Ls; ++l)
	    p[l][A.subset()] = z[l] + beta*p[l]; 
	}
	//Cheb.Qsq(Ap,p);
	A(Ap,p,PLUS);
	pAp = innerProductReal(p,Ap,A.subset());
      
	alpha = r_dot_z/pAp;
	for(int l=0; l < Ls; ++l)
	  x[l][A.subset()] += alpha*p[l];
	for(int l=0; l < Ls; ++l)
	  r[l][A.subset()] -= alpha*Ap[l];
	r_norm2 =  norm2(r,A.subset());
	r_dot_z_old = r_dot_z;

	if(k>MaxCG){
	  res.n_count = k;
	  QDP_error_exit("too many CG iterations: count = %d", res.n_count);
	  return res;
	}
#if 1
	QDPIO::cout << "vecPrecondCG: k = " << k << "  res^2 = " << r_norm2;
	QDPIO::cout << "  r_dot_z = " << r_dot_z << endl;
#endif
      }
    
      res.n_count = k;
      res.resid   = sqrt(r_norm2); 
      swatch.stop();
      QDPIO::cout << "vPreconfCG: k = " << k << endl;
      flopcount.report("vPrecondCG", swatch.getTimeInSeconds());
      END_CODE();
      return res;
    }
  
    template<typename T>
    void InitGuess_T(const LinearOperatorArray<T>& A, 
		     multi1d<T>& x, 
		     const multi1d<T>& b, 
		     const multi1d<Double>& eval, 
		     const multi2d<T>& evec, 
		     int& n_count)
    {
      int N = evec.size2();
      InitGuess(A,x,b,eval,evec,N,n_count);
    }

    template<typename T>
    void InitGuess_T(const LinearOperatorArray<T>& A, 
		     multi1d<T>& x, 
		     const multi1d<T>& b, 
		     const multi1d<Double>& eval, 
		     const multi2d<T>& evec, 
		     int N, // number of vectors to use
		     int& n_count)
    {
      multi1d<T> p; 
      multi1d<T> Ap; 
      multi1d<T> r;

      const int Ls = A.size();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      A(Ap,x,PLUS);
      for(int l=0; l < Ls; ++l)
	r[l][A.subset()] = b[l] - Ap[l];
      // Double r_norm2 = norm2(r,A.subset());
   
      for(int i(0);i<N;i++)
      {
	DComplex d = innerProduct(evec[i],r,A.subset());
	for(int l=0; l < Ls; ++l)
	  x[l][A.subset()] += (d/eval[i])*evec[i][l];
	//QDPIO::cout<<"InitCG: "<<d<<" "<<eval[i]<<endl;

      }
   
      snoop.stop();
      QDPIO::cout << "InitGuess:  time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      n_count = 1;
    }

    
    //
    // Wrappers
    //
    // LatticeFermionF
    void SubSpaceMatrix(LinAlg::Matrix<DComplex>& H,
			const LinearOperatorArray<LatticeFermionF>& A,
			const multi2d<LatticeFermionF>& evec,
			int Nvecs)
    {
      SubSpaceMatrix_T(H, A, evec, Nvecs);
    }

    SystemSolverResults_t InvEigCG2(const LinearOperatorArray<LatticeFermionF>& A,
				    multi1d<LatticeFermionF>& x, 
				    const multi1d<LatticeFermionF>& b,
				    multi1d<Double>& eval, 
				    multi2d<LatticeFermionF>& evec,
				    int Neig,
				    int Nmax,
				    const Real& RsdCG, int MaxCG)
    {
      return InvEigCG2_T(A, x, b, eval, evec, Neig, Nmax, RsdCG, MaxCG);
    }
  
    SystemSolverResults_t vecPrecondCG(const LinearOperatorArray<LatticeFermionF>& A, 
				       multi1d<LatticeFermionF>& x, 
				       const multi1d<LatticeFermionF>& b, 
				       const multi1d<Double>& eval, 
				       const multi2d<LatticeFermionF>& evec, 
				       int startV, int endV,
				       const Real& RsdCG, int MaxCG)
    {
      return vecPrecondCG_T(A, x, b, eval, evec, startV, endV, RsdCG, MaxCG);
    }

    void InitGuess(const LinearOperatorArray<LatticeFermionF>& A, 
		   multi1d<LatticeFermionF>& x, 
		   const multi1d<LatticeFermionF>& b, 
		   const multi1d<Double>& eval, 
		   const multi2d<LatticeFermionF>& evec, 
		   int& n_count)
    {
      InitGuess_T(A, x, b, eval, evec, n_count);
    }
  
    void InitGuess(const LinearOperatorArray<LatticeFermionF>& A, 
		   multi1d<LatticeFermionF>& x, 
		   const multi1d<LatticeFermionF>& b, 
		   const multi1d<Double>& eval, 
		   const multi2d<LatticeFermionF>& evec, 
		   int N, // number of vectors to use
		   int& n_count)
    {
      InitGuess_T(A, x, b, eval, evec, N, n_count);
    }


    // LatticeFermionD
    void SubSpaceMatrix(LinAlg::Matrix<DComplex>& H,
			const LinearOperatorArray<LatticeFermionD>& A,
			const multi2d<LatticeFermionD>& evec,
			int Nvecs)
    {
      SubSpaceMatrix_T(H, A, evec, Nvecs);
    }

    SystemSolverResults_t InvEigCG2(const LinearOperatorArray<LatticeFermionD>& A,
				    multi1d<LatticeFermionD>& x, 
				    const multi1d<LatticeFermionD>& b,
				    multi1d<Double>& eval, 
				    multi2d<LatticeFermionD>& evec,
				    int Neig,
				    int Nmax,
				    const Real& RsdCG, int MaxCG)
    {
      return InvEigCG2_T(A, x, b, eval, evec, Neig, Nmax, RsdCG, MaxCG);
    }
  
    SystemSolverResults_t vecPrecondCG(const LinearOperatorArray<LatticeFermionD>& A, 
				       multi1d<LatticeFermionD>& x, 
				       const multi1d<LatticeFermionD>& b, 
				       const multi1d<Double>& eval, 
				       const multi2d<LatticeFermionD>& evec, 
				       int startV, int endV,
				       const Real& RsdCG, int MaxCG)
    {
      return vecPrecondCG_T(A, x, b, eval, evec, startV, endV, RsdCG, MaxCG);
    }

    void InitGuess(const LinearOperatorArray<LatticeFermionD>& A, 
		   multi1d<LatticeFermionD>& x, 
		   const multi1d<LatticeFermionD>& b, 
		   const multi1d<Double>& eval, 
		   const multi2d<LatticeFermionD>& evec, 
		   int& n_count)
    {
      InitGuess_T(A, x, b, eval, evec, n_count);
    }
  
    void InitGuess(const LinearOperatorArray<LatticeFermionD>& A, 
		   multi1d<LatticeFermionD>& x, 
		   const multi1d<LatticeFermionD>& b, 
		   const multi1d<Double>& eval, 
		   const multi2d<LatticeFermionD>& evec, 
		   int N, // number of vectors to use
		   int& n_count)
    {
      InitGuess_T(A, x, b, eval, evec, N, n_count);
    }
  } // namespace InvEigCG2ArrayEnv
  
}// End Namespace Chroma

