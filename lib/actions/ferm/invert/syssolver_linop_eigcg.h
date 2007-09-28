// -*- C++ -*-
// $Id: syssolver_linop_eigcg.h,v 1.2 2007-09-28 02:02:47 kostas Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system by CG2
 */

#ifndef __syssolver_linop_eigcg_h__
#define __syssolver_linop_eigcg_h__
#include "chroma_config.h"
#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "actions/ferm/invert/syssolver_linop.h"
#include "actions/ferm/invert/syssolver_eigcg_params.h"
#include "actions/ferm/invert/inv_eigcg2.h"


namespace Chroma
{

  //! Eigenvector accelerated CG system solver namespace
  namespace LinOpSysSolverEigCGEnv
  {
    //! Name to be used
    extern const std::string name;

    //! Register the syssolver
    bool registerAll();
  }


  //! Solve a M*psi=chi linear system by CG2 with eigenvectors
  /*! \ingroup invert
   */
  template<typename T>
  class LinOpSysSolverEigCG : public LinOpSystemSolver<T>
  {
  public:
    //! Constructor
    /*!
     * \param M_        Linear operator ( Read )
     * \param invParam  inverter parameters ( Read )
     */
    LinOpSysSolverEigCG(Handle< LinearOperator<T> > A_,
		 const SysSolverEigCGParams& invParam_) : 
      A(A_), invParam(invParam_) 
      {
	if(invParam.Neig_max>0){
	  GoodEval.vec.resize(invParams.Neig_max) ;
	  GoodEvec.vec.resize(invParams.Neig_max) ;
	  GoodEval.N = 0 ;
	  GoodEvec.N = 0 ;
	}
	else{
	  GoodEval.vec.resize(invParams.Neig) ;
	  GoodEvec.vec.resize(invParams.Neig) ;
	  GoodEval.N = 0 ;
	  GoodEvec.N = 0 ;
	}
	// NEED to grab the eignvectors from the named buffer here
      }

    //! Destructor is automatic
    ~LinOpSysSolverEigCG() {}

    //! Return the subset on which the operator acts
    const Subset& subset() const {return A->subset();}

    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     */
    SystemSolverResults_t operator() (T& psi, const T& chi) const
      {
	START_CODE();

	T chi_tmp;
	(*A)(chi_tmp, chi, MINUS);

	multi1d<Double> lambda ; //the eigenvalues
	multi1d<T> evec(0); // The eigenvectors  
	SystemSolverResults_t res;  // initialized by a constructor
	int n_CG(0);
	// Need to pass the appropriate operators here...
	InitGuess(*A,psi,chi_tmp,GoodEval.vec,GoodEvec.vec,GoodEvec.N,n_CG);
	if((GoodEvec.N+invParams.Neig)<=invParams.Neig_max){// if there is space for new
	  evec.resize(0);// get in there with no evecs so that it computes new
	  res = InvEigCG2(*A, psi, chi_tmp, eval, evec, 
			  invParam.Neig, invParam.Nmax, 
			  invParam.RsdCG, invParam.MaxCG);
	  res.n_count += n_CG ;
	  
	  snoop.start();
	  GoodEvec.AddVectors(evec,  s);
	  GoodEval.AddVectors(lambda,s);
	  snoop.stop();
	  Time = snoop.getTimeInSeconds() ;
	  QDPIO::cout<<"GoodEval.N= "<<GoodEval.N<<endl;
	  
	  snoop.start();
	  SimpleGramSchmidt(GoodEvec.vec,GoodEvec.N-params.Neig,GoodEvec.N,s);
	  SimpleGramSchmidt(GoodEvec.vec,GoodEvec.N-params.Neig,GoodEvec.N,s);
	  snoop.stop();
	  Time += snoop.getTimeInSeconds() ;
	  
	  snoop.start();
	  Matrix<DComplex> Htmp(GoodEvec.N) ;
	  SubSpaceMatrix(Htmp,*MM,GoodEvec.vec,GoodEvec.N);
	  //OctavePrintOut(Htmp.mat,Htmp.N,tag("H",i),"RayleighRich.m") ;
	  char V = 'V' ; char U = 'U' ;
	  Lapack::zheev(V,U,Htmp.mat,lambda);
	  evec.resize(GoodEvec.N) ;
	  //OctavePrintOut(Htmp.mat,Htmp.N,tag("Htmp",i),"RayleighRich.m") ;
	  for(int k(0);k<GoodEvec.N;k++){
	    GoodEval[k] = lambda[k];
	    evec[k][s] = zero ;
	    //cout<<k<<endl ;
	    for(int j(0);j<GoodEvec.N;j++)
	      evec[k][s] += conj(Htmp(k,j))*GoodEvec[j] ;
	    //QDPIO::cout<<"norm(evec["<<k<<"])=";
	    //QDPIO::cout<< sqrt(norm2(GoodEvec[k],s))<<endl;
	  }
	  snoop.stop();
	  Time += snoop.getTimeInSeconds() ;
	  QDPIO::cout << "Evec_Refinement: time = "
		      << Time
		      << " secs" << endl;
	  
	  QDPIO::cout<<"GoodEval.N= "<<GoodEval.N<<endl;
	  QDPIO::cout<<"GoodEvec.N= "<<GoodEvec.N<<endl;
	  for(int k(0);k<GoodEvec.N;k++)
	    GoodEvec[k][s]  = evec[k] ;
	  
	  //Check the quality of eigenvectors
#ifdef TEST_ALGORITHM
	  {
	    LatticeFermion Av ;
	    for(int k(0);k<GoodEvec.N;k++){
	      (*MM)(Av,GoodEvec[k],PLUS) ;
	      DComplex rq = innerProduct(GoodEvec[k],Av,s);
	      Av[s] -= GoodEval[k]*GoodEvec[k] ;
	      Double tt = sqrt(norm2(Av,s));
	      QDPIO::cout<<"REFINE: error evec["<<k<<"] = "<<tt<<" " ;
	      QDPIO::cout<<"--- eval ="<<GoodEval[k]<<" ";
	      tt =  sqrt(norm2(GoodEvec[k],s));
	      QDPIO::cout<<"--- rq ="<<real(rq)<<" ";
	      QDPIO::cout<<"--- norm = "<<tt<<endl  ;
	    } 
	  }
#endif
	}// if there is space
	else{// call CG but ask it not to compute vectors
	  evec.resize(0);
	  n_CG = res.n_count ;
	  if(invParams.vPrecCGvecs ==0){
	    res =InvEigCG2(*A, 
			   psi,
			   chi_tmp,
			   lambda, 
			   evec, 
			   0, //Eigenvectors to keep
			   invParams.Nmax,  // Max vectors to work with
			   invParams.RsdCGRestart, // CG residual.... Empirical restart need a param here...
			   invParams.MaxCG // Max CG itterations
			   ); 
	  }
	  else
	    res= vecPrecondCG(*A,psi,chi_tmp,GoodEval.vec,GoodEvec.vec,
			      invParams.vPrecCGvecStart,
			      invParams.vPrecCGvecStart+invParams.vPrecCGvecs,
                              invParams.RsdCG, // CG residual       
                              invParams.MaxCG // Max CG itterations    
			      );
	  res.n_count += n_CG ;
	  n_CG=0 ;
	  QDPIO::cout<<"Restart1: "<<endl ;
	  InitGuess(*A,psi,chi_tmp,GoodEval.vec,GoodEvec.vec,GoodEvec.N, n_CG);
	  res.n_count += n_CG ;
	  n_CG = res.n_count ;
	  if(invParams.vPrecCGvecs ==0)
	    res =InvEigCG2(*A,
                           psi,
                           chi_tmp,
                           lambda,
                           evec,
                           0, //Eigenvectors to keep                             
                           invParams.Nmax,  // Max vectors to work with           
                           invParams.RsdCG, // CG residual....
                           invParams.MaxCG // Max CG itterations             
                           );
	  else
	    res= vecPrecondCG(*A,psi,chi_tmp,GoodEval.vec,GoodEvec.vec,
			      invParams.vPrecCGvecStart,
                              invParams.vPrecCGvecStart+invParams.vPrecCGvecs,
			      invParams.RsdCG, // CG residual             
			      invParams.MaxCG // Max CG itterations          
                              );
	  res.n_count += n_CG ;
	}

	END_CODE();

	return res;
      }

    // IN THE DISTRACTOR I need to store back to the name buffer the computed eigenvectors

  private:
    // Hide default constructor
    LinOpSysSolverEigCG() {}

    Handle< LinearOperator<T> > A;
    SysSolverEigCGParams invParam;

    Vectors<T> GoodEvec; //
    Vectors<Double> GoodEval; //

  };

} // End namespace

#endif 

