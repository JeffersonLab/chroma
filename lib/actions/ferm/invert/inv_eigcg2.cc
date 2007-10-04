// -*- C++ -*-
// $Id: inv_eigcg2.cc,v 1.3 2007-10-04 20:39:57 kostas Exp $
#include <sstream>
#include "inv_stathoCG_w.h"
#include "octave_debug.h"

//#include "simpleGramSchmidt.h"
//#include "lapack_wrapper.h"
//#include "containers.h"

//#define DEBUG
#define DEBUG_FINAL


// NEEDS A LOT OF CLEAN UP
namespace Chroma 
{

  using namespace LinAlg ;
  using namespace Octave ;

  //typedef LatticeFermion T ; // save sometyping 

  // needs a little clean up... for possible exceptions if the size of H
  // is smaller than hind 
  void SubSpaceMatrix(Matrix<DComplex>& H,
		      LinearOperator<LatticeFermion>& A,
		      const multi1d<LatticeFermion>& evec,
		      const int Nvecs){
    
    OrderedSubset s(A.subset());
    
    LatticeFermion Ap ;
    H.N = Nvecs ;
    if(Nvecs>evec.size()){
      H.N = evec.size();
    }
    if(H.N>H.size()){
      QDPIO::cerr<<"OOPS! your matrix can't take this many matrix elements\n";
      exit(1);	
    }
    // fill H  with zeros
    H.mat = 0.0 ;
    
    for(int i(0);i<H.N;i++){
      A(Ap,evec[i],PLUS) ;
      for(int j(i);j<H.N;j++){
	H(j,i) = innerProduct(evec[j], Ap, s) ;
	//enforce hermiticity
	H(i,j) = conj(H(j,i));
	if(i==j) H(i,j) = real(H(i,j));
      }
    } 
  }
   
  void SubSpaceMatrix(Matrix<DComplex>& H,
		      LinearOperator<LatticeFermion>& A,
		      const std::vector< RitzPair<LatticeFermion> >& ritz,
		      const int Nvecs){
    
    OrderedSubset s(A.subset());
    
    LatticeFermion Ap ;
    H.N = Nvecs ;
    if(Nvecs>ritz.size()){
      H.N = ritz.size();
    }
    if(H.N>H.size()){
      QDPIO::cerr<<"OOPS! your matrix can't take this many matrix elements\n";
      exit(1);	
    }
    // fill H  with zeros
    H.mat = 0.0 ;
    
    for(int i(0);i<H.N;i++){
      A(Ap,ritz[i].evec,PLUS) ;
      for(int j(i);j<H.N;j++){
	H(j,i) = innerProduct(ritz[j].evec, Ap, s) ;
	//enforce hermiticity
	H(i,j) = conj(H(j,i));
	if(i==j) H(i,j) = real(H(i,j));
      }
    } 
  }
 
  SystemSolverResults_t InvEigCG2(LinearOperator<LatticeFermion>& A,
				  LatticeFermion& x, 
				  const LatticeFermion& b,
				  multi1d<Double>& eval, 
				  multi1d<LatticeFermion>& evec,
				  const int Neig,
				  const int Nmax,
				  const Real RsdCG, const int MaxCG)
  {
    StopWatch snoop;
    snoop.reset();
    snoop.start();
    
    SystemSolverResults_t  res;
    OrderedSubset s(A.subset());
    
    LatticeFermion p ; 
    LatticeFermion Ap; 
    LatticeFermion r,z ;

    Double rsd_sq = (RsdCG * RsdCG) * Real(norm2(b,s));
    Double alphaprev, alpha,pAp;
    Real beta ;
    Double r_dot_z, r_dot_z_old ;
    //Complex r_dot_z, r_dot_z_old,beta ;
    //Complex alpha,pAp ;

    int k = 0 ;
    A(Ap,x,PLUS) ;
    r[s] = b - Ap ;
    Double r_norm2 = norm2(r,s) ;

#if 1
    QDPIO::cout << "StathoCG: Nevecs(input) = " << evec.size() << endl;
    QDPIO::cout << "StathoCG: k = " << k << "  res^2 = " << r_norm2 << endl;
#endif
    Real inorm(Real(1.0/sqrt(r_norm2)));
    bool FindEvals = (Neig>0);
    int tr; // Don't know what this does...
    bool from_restart;
    Matrix<DComplex> H(Nmax) ; // square matrix containing the tridiagonal
    Vectors<LatticeFermion> vec(Nmax) ; // contains the vectors we use...
    //-------- Eigenvalue eigenvector finding code ------
    vec.AddVectors(evec,s) ;
    p[s] = inorm*r ; // this is not needed GramSchmidt will take care of it
    vec.AddOrReplaceVector(p,s);

    SimpleGramSchmidt(vec.vec,vec.N,vec.N+1,s);
    SimpleGramSchmidt(vec.vec,vec.N,vec.N+1,s); //call twice: need to improve...
    SubSpaceMatrix(H,A,vec.vec,vec.N);
    from_restart = true ;
    tr = H.N - 1; // this is just a flag as far as I can tell  
    //-------

    Double betaprev ;
    beta=0.0;
    alpha = 1.0;
    // Algorithm from page 529 of Golub and Van Loan
    // Modified to match the m-code from A. Stathopoulos
    while(toBool(r_norm2>rsd_sq)){
      /** preconditioning algorithm **/
      z[s]=r ; //preconditioning can be added here
      /**/
      r_dot_z = innerProductReal(r,z,s);
      k++ ;
      betaprev = beta; // additional line for Eigenvalue eigenvector code
      if(k==1){
	p[s] = z ;	
	H.N++ ;
      }
      else{
	beta = r_dot_z/r_dot_z_old ;
	p[s] = z + beta*p ; 
	//-------- Eigenvalue eigenvector finding code ------
	// fist block
	if(FindEvals){
	  if(!((from_restart)&&(H.N == tr+1))){
	    if(!from_restart){
	      H(H.N-2,H.N-2) = 1/alpha + betaprev/alphaprev;
	    } 
	    from_restart = false ; 
	    H(H.N-1,H.N-2) = -sqrt(beta)/alpha;
	    H(H.N-2,H.N-1) = -sqrt(beta)/alpha;
	  }
	  H.N++ ;
	}
	//---------------------------------------------------
      }
      A(Ap,p,PLUS) ;
      pAp = innerProductReal(p,Ap,s);
      
      alphaprev = alpha ;// additional line for Eigenvalue eigenvector code
      alpha = r_dot_z/pAp ;
      x[s] += alpha*p ;
      r[s] -= alpha*Ap ;
      r_norm2 =  norm2(r,s) ;
      r_dot_z_old = r_dot_z ;

      
      //-------- Eigenvalue eigenvector finding code ------
      // second block
      if(FindEvals){
	if (vec.N==Nmax){//we already have stored the maximum number of vectors
	  // The magic begins here....
QDPIO::cout<<"MAGIC BEGINS: H.N ="<<H.N<<endl ;
	  H(Nmax-1,Nmax-1) = 1/alpha + beta/alphaprev;

#ifdef DEBUG
	  {
	    stringstream tag ;
	    tag<<"H"<<k ;
	    OctavePrintOut(H.mat,Nmax,tag.str(),"Hmatrix.m");
	  }
	  {
	    Matrix<DComplex> tmp(Nmax) ; 
	    SubSpaceMatrix(tmp,A,vec.vec,vec.N);
	    stringstream tag ;
	    tag<<"H"<<k<<"ex" ;
	    OctavePrintOut(tmp.mat,Nmax,tag.str(),"Hmatrix.m");
	  }
#endif
	  //exit(1);
	  multi2d<DComplex> Hevecs(H.mat) ;
	  multi1d<Double> Heval ;
	  char V = 'V' ; char U = 'U' ;
	  Lapack::zheev(V,U,Hevecs,Heval);
	  multi2d<DComplex> Hevecs_old(H.mat) ;

#ifdef DEBUG
 {
   multi1d<LatticeFermion> tt_vec(vec.size());
	  for(int i(0);i<Nmax;i++){
	    tt_vec[i][s] = zero ;
	    for(int j(0);j<Nmax;j++)
	      tt_vec[i][s] +=Hevecs(i,j)*vec[j] ; // NEED TO CHECK THE INDEXINT
	  }
	  // CHECK IF vec are eigenvectors...
	  
	    LatticeFermion Av ;
	    for(int i(0);i<Nmax;i++){
	      A(Av,tt_vec[i],PLUS) ;
	      DComplex rq = innerProduct(tt_vec[i],Av,s);
	      Av[s] -= Heval[i]*tt_vec[i] ;
	      Double tt = sqrt(norm2(Av,s));
	      QDPIO::cout<<"1 error eigenvector["<<i<<"] = "<<tt<<" " ;
	      tt =  sqrt(norm2(tt_vec[i],s));
	      QDPIO::cout<<" --- eval = "<<Heval[i]<<" ";
	      QDPIO::cout<<" --- rq = "<<real(rq)<<" ";
	      QDPIO::cout<<"--- norm = "<<tt<<endl  ;
	    }
 }
#endif
	  multi1d<Double> Heval_old ;
	  Lapack::zheev(V,U,Nmax-1,Hevecs_old,Heval_old);
	  for(int i(0);i<Nmax;i++)	    
	    Hevecs_old(i,Nmax-1) = Hevecs_old(Nmax-1,i) = 0.0 ;
#ifdef DEBUG
	  {
	    stringstream tag ;
	    tag<<"Hevecs_old"<<k ;
	    OctavePrintOut(Hevecs_old,Nmax,tag.str(),"Hmatrix.m");
	  }
	  
#endif
	  tr = Neig + Neig ; // v_old = Neig optimal choice
	   
	  for(int i(Neig);i<tr;i++)
	    for(int j(0);j<Nmax;j++)	    
	      Hevecs(i,j) = Hevecs_old(i-Neig,j) ;
#ifdef DEBUG
	  {
	    stringstream tag ;
	    tag<<"Hevecs"<<k ;
	    OctavePrintOut(Hevecs,Nmax,tag.str(),"Hmatrix.m");
	  }
#endif

	  // Orthogonalize the last Neig vecs (keep the first Neig fixed)
	  // zgeqrf(Nmax, 2*Neig, Hevecs, Nmax,
	  //    TAU_CMPLX_ofsize_2Neig, Workarr, 2*Neig*Lapackblocksize, info);
	  multi1d<DComplex> TAU ;
	  Lapack::zgeqrf(Nmax,2*Neig,Hevecs,TAU) ;
	  char R = 'R' ; char L = 'L' ; char N ='N' ; char C = 'C' ;
	  multi2d<DComplex> Htmp = H.mat ;
	  Lapack::zunmqr(R,N,Nmax,Nmax,Hevecs,TAU,Htmp);
	  Lapack::zunmqr(L,C,Nmax,2*Neig,Hevecs,TAU,Htmp);
	  // Notice that now H is of size 2Neig x 2Neig, 
	  // but still with LDA = Nmax 
#ifdef DEBUG
	  {
	    stringstream tag ;
	    tag<<"Htmp"<<k ;
	    OctavePrintOut(Htmp,Nmax,tag.str(),"Hmatrix.m");
	  }
#endif
          Lapack::zheev(V,U,2*Neig,Htmp,Heval);
 	  for(int i(Neig); i< 2*Neig;i++ ) 
	     for(int j(2*Neig); j<Nmax; j++)
                Htmp(i,j) =0.0;
#ifdef DEBUG
          {
            stringstream tag ;
            tag<<"evecstmp"<<k ;
            OctavePrintOut(Htmp,Nmax,tag.str(),"Hmatrix.m");
          }
#endif
	  Lapack::zunmqr(L,N,Nmax,2*Neig,Hevecs,TAU,Htmp);
#ifdef DEBUG
          {
            stringstream tag ;
            tag<<"Yevecstmp"<<k ;
            OctavePrintOut(Htmp,Nmax,tag.str(),"Hmatrix.m");
          }
#endif
	  // Here I need the restart_X bit
	  multi1d<LatticeFermion> tt_vec = vec.vec;
	  for(int i(0);i<2*Neig;i++){
	    vec[i][s] = zero ;
	    for(int j(0);j<Nmax;j++)
	      vec[i][s] +=Htmp(i,j)*tt_vec[j] ; // NEED TO CHECK THE INDEXING
	  }

#ifdef DEBUG
	  // CHECK IF vec are eigenvectors...
	  {
	    LatticeFermion Av ;
	    for(int i(0);i<2*Neig;i++){
	      A(Av,vec[i],PLUS) ;
	      DComplex rq = innerProduct(vec[i],Av,s);
	      Av[s] -= Heval[i]*vec[i] ;
	      Double tt = sqrt(norm2(Av,s));
	      QDPIO::cout<<"error eigenvector["<<i<<"] = "<<tt<<" " ;
	      tt =  sqrt(norm2(vec[i],s));
	      QDPIO::cout<<"--- rq ="<<real(rq)<<" ";
	      QDPIO::cout<<"--- norm = "<<tt<<endl  ;
	    }

	  }
#endif
	  H.mat = 0.0 ; // zero out H 
	  for (int i=0;i<2*Neig;i++) H(i,i) = Heval[i];
	  
	  LatticeFermion Ar ;
	  // Need to reorganize so that we aboid this extra matvec
	  A(Ar,r,PLUS) ;
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
	    H(2*Neig,i) = innerProduct(vec[i], Ar, s) ;
	    H(i,2*Neig) = conj(H(2*Neig,i)) ;
	  }
	  H(2*Neig,2*Neig) = innerProduct(r, Ar, s)/sqrt(r_norm2) ;
	  
	  H.N = 2*Neig + 1  ; // why this ?
	  from_restart = true ;
	  vec.N = 2*Neig ; // only keep the lowest Neig. Is this correct???
#ifdef DEBUG
	  {
	    stringstream tag ;
	    tag<<"finalH"<<k ;
	    OctavePrintOut(H.mat,Nmax,tag.str(),"Hmatrix.m");
	  }
#endif
	}
	// Here we add a vector in the list
	Double inorm = 1.0/sqrt(r_norm2) ;
	vec.NormalizeAndAddVector(r,inorm,s) ;
	// Shouldn't be the z vector when a preconditioner is used?
      }
      //---------------------------------------------------

      if(k>MaxCG){
	res.n_count = k;
	res.resid   = sqrt(r_norm2);
	QDP_error_exit("too many CG iterations: count = %d", res.n_count);
	return res;
      }
#if 1
      QDPIO::cout << "StathoCG: k = " << k << "  res^2 = " << r_norm2 ;
      QDPIO::cout << "  r_dot_z = " << r_dot_z << endl;
#endif
    }
    res.n_count = k ;
    snoop.stop();
    Double snoop1_time(0.0);
    if(FindEvals){
      StopWatch snoop1 ;
      snoop1.reset();
      snoop1.start();
      // Evec Code ------ Before we return --------
      // vs is the current number of vectors stored
      // Neig is the number of eigenvector we want to compute
      SimpleGramSchmidt(vec.vec,0,2*Neig,s);
      SimpleGramSchmidt(vec.vec,0,2*Neig,s);
      Matrix<DComplex> Htmp(2*Neig) ;
      SubSpaceMatrix(Htmp,A,vec.vec,2*Neig);

      char V = 'V' ; char U = 'U' ;
      multi1d<Double> tt_eval ;
      Lapack::zheev(V,U,Htmp.mat,tt_eval);
      evec.resize(Neig) ;
      eval.resize(Neig);
      for(int i(0);i<Neig;i++){
	evec[i][s] = zero ;
	eval[i] = tt_eval[i] ;
	for(int j(0);j<2*Neig;j++)
	  evec[i][s] += Htmp(i,j)*vec[j] ;
      }
      
      // I will do the checking of eigenvector quality outside this routine
      // -------------------------------------------
      snoop1.stop(); // do not time the checking...
      snoop1_time = snoop1.getTimeInSeconds() ;
#ifdef DEBUG_FINAL
      // CHECK IF vec are eigenvectors...                                   
      {
	LatticeFermion Av ;
	for(int i(0);i<Neig;i++){
	  A(Av,evec[i],PLUS) ;
	  DComplex rq = innerProduct(evec[i],Av,s);
	  Av[s] -= eval[i]*evec[i] ;
	  Double tt = sqrt(norm2(Av,s));
	  QDPIO::cout<<"FINAL: error eigenvector["<<i<<"] = "<<tt<<" " ;
	  tt =  sqrt(norm2(evec[i],s));
	  QDPIO::cout<<"--- rq ="<<real(rq)<<" ";
	  QDPIO::cout<<"--- norm = "<<tt<<endl  ;
	}
	
      }
#endif
    }

    Double time = (snoop.getTimeInSeconds() + snoop1_time) ;
    QDPIO::cout << "InvStathoCG: time = "
		<< time 
		<< " secs" << endl;
    QDPIO::cout << "InvStathoCG: time per iteration = "
		<< time/res.n_count 
		<< " secs/iter" << endl;

    res.resid   = sqrt(r_norm2);

    return res;
  }

  // A  should be Hermitian positive definite
  SystemSolverResults_t vecPrecondCG(LinearOperator<LatticeFermion>& A, 
				     LatticeFermion& x, 
				     const LatticeFermion& b, 
				     const multi1d<Double>& eval, 
				     const multi1d<LatticeFermion>& evec, 
				     const startV,const endV,
				     const Real RsdCG, const int MaxCG){


    StopWatch snoop;
    snoop.reset();
    snoop.start();
    if(endV>eval.size()){
      QDP_error_exit("vPrecondCG: not enought evecs eval.size()=%d",eval.size());
    }
    OrderedSubset s(A.subset());
 
    LatticeFermion p ; 
    LatticeFermion Ap; 
    LatticeFermion r,z ;

    Double rsd_sq = (RsdCG * RsdCG) * Real(norm2(b,s));
    Double alpha,pAp;
    Real beta ;
    Double r_dot_z, r_dot_z_old ;
    //Complex r_dot_z, r_dot_z_old,beta ;
    //Complex alpha,pAp ;

    int k = 0 ;
    A(Ap,x,PLUS) ;
    r[s] = b - Ap ;
    Double r_norm2 = norm2(r,s) ;

#if 1
    QDPIO::cout << "vecPrecondCG: k = " << k << "  res^2 = " << r_norm2 << endl;
#endif

    // Algorithm from page 529 of Golub and Van Loan
    while(toBool(r_norm2>rsd_sq)){
      /** preconditioning algorithm **/
      z[s]=r ; // not optimal but do it for now...
      /**/
      for(int i(startV);i<endV;i++){
	DComplex d = innerProduct(evec[i],r,s) ;
	//QDPIO::cout<<"vecPrecondCG: "<< d<<" "<<(1.0/eval[i]-1.0)<<endl;
	z[s] += (1.0/eval[i]-1.0)*d*evec[i];
      }
      /**/
      /**/
      r_dot_z = innerProductReal(r,z,s);
      k++ ;
      if(k==1){
	p[s] = z ;	
      }
      else{
	beta = r_dot_z/r_dot_z_old ;
	p[s] = z + beta*p ; 
      }
      //Cheb.Qsq(Ap,p) ;
      A(Ap,p,PLUS) ;
      pAp = innerProductReal(p,Ap,s);
      
      alpha = r_dot_z/pAp ;
      x[s] += alpha*p ;
      r[s] -= alpha*Ap ;
      r_norm2 =  norm2(r,s) ;
      r_dot_z_old = r_dot_z ;

      if(k>MaxCG){
	res.n_count = k ;
	QDP_error_exit("too many CG iterations: count = %d", res.n_count);
	return ;
      }
#if 1
      QDPIO::cout << "vecPrecondCG: k = " << k << "  res^2 = " << r_norm2 ;
      QDPIO::cout << "  r_dot_z = " << r_dot_z << endl;
#endif
    }
    
    res.n_count = k ;
    res.resid   = sqrt(r_norm2); 

    snoop.stop();
    QDPIO::cout << "vPrecondCG: time = "
		<< time 
		<< " secs" << endl;
    QDPIO::cout << "vPrecondCG: time per iteration = "
		<< time/res.n_count 
		<< " secs/iter" << endl;

    return res ;
  }
  
void InitGuess(LinearOperator<LatticeFermion>& A, 
		LatticeFermion& x, 
		const LatticeFermion& b, 
		const multi1d<Double>& eval, 
		const multi1d<LatticeFermion>& evec, 
		int& n_count){

  int N = evec.size();
  InitGuess(A,x,b,eval,evec,N,n_count);
 }

 void InitGuess(LinearOperator<LatticeFermion>& A, 
		LatticeFermion& x, 
		const LatticeFermion& b, 
		const multi1d<Double>& eval, 
		const multi1d<LatticeFermion>& evec, 
		const int N, // number of vectors to use
		int& n_count){

   OrderedSubset s(A.subset());

   LatticeFermion p ; 
   LatticeFermion Ap; 
   LatticeFermion r ;

   StopWatch snoop;
   snoop.reset();
   snoop.start();

   A(Ap,x,PLUS) ;
   r[s] = b - Ap ;
   // Double r_norm2 = norm2(r,s) ;
   
   for(int i(0);i<N;i++){
     DComplex d = innerProduct(evec[i],r,s) ;
     x[s] += (d/eval[i])*evec[i];
     //QDPIO::cout<<"InitCG: "<<d<<" "<<eval[i]<<endl ;

   }
   
   snoop.stop();
   QDPIO::cout << "InitGuess:  time = "
              << snoop.getTimeInSeconds() 
              << " secs" << endl;

   n_count = 1 ;
 }

  
}// End Namespace Chroma

