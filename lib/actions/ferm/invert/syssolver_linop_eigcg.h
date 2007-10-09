// -*- C++ -*-
// $Id: syssolver_linop_eigcg.h,v 1.5 2007-10-09 18:07:14 edwards Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system by CG2
 */

#ifndef __syssolver_linop_eigcg_h__
#define __syssolver_linop_eigcg_h__

#include "handle.h"
#include "syssolver.h"
#include "linearop.h"
#include "named_obj.h"
#include "meas/inline/io/named_objmap.h"

#include "actions/ferm/invert/syssolver_linop.h"
#include "actions/ferm/invert/syssolver_eigcg_params.h"
#include "actions/ferm/invert/inv_eigcg2.h"

#include "actions/ferm/invert/containers.h"
#include "actions/ferm/invert/lapack_wrapper.h"


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
     * \param M_         Linear operator ( Read )
     * \param invParam_  inverter parameters ( Read )
     */
    LinOpSysSolverEigCG(Handle< LinearOperator<T> > A_,
			const SysSolverEigCGParams& invParam_) : 
      A(A_), invParam(invParam_) 
      {
	// NEED to grab the eignvectors from the named buffer here
	if (! TheNamedObjMap::Instance().check(invParam.eigen_id))
	{
	  TheNamedObjMap::Instance().create< LinAlg::RitzPairs<T> >(invParam.eigen_id);
	  LinAlg::RitzPairs<T>& GoodEvecs = 
	    TheNamedObjMap::Instance().getData< LinAlg::RitzPairs<T> >(invParam.eigen_id);

	  if(invParam.Neig_max>0 ){
	    GoodEvecs.init(invParam.Neig_max);
	  }
	  else{
	    GoodEvecs.init(invParam.Neig);
	  }
	}
      }

    //! Destructor is automatic
    ~LinOpSysSolverEigCG()
      {
	if (invParam.cleanUpEvecs)
	{
	  TheNamedObjMap::Instance().erase(invParam.eigen_id);
	}
      }

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

	LinAlg::RitzPairs<T>& GoodEvecs = TheNamedObjMap::Instance().getData< LinAlg::RitzPairs<T> >(invParam.eigen_id);

	multi1d<Double> lambda ; //the eigenvalues
	multi1d<T> evec(0); // The eigenvectors  
	SystemSolverResults_t res;  // initialized by a constructor
	int n_CG(0);
	// Need to pass the appropriate operators here...
	InitGuess(*A,psi,chi_tmp,GoodEvecs.eval.vec,GoodEvecs.evec.vec,GoodEvecs.Neig,n_CG);
	if((GoodEvecs.Neig+invParam.Neig)<=invParam.Neig_max) // if there is space for new
	{
	  StopWatch snoop;
	  snoop.reset();
	  snoop.start();

	  evec.resize(0);// get in there with no evecs so that it computes new
	  res = InvEigCG2(*A, psi, chi_tmp, lambda, evec, 
			  invParam.Neig, invParam.Nmax, 
			  invParam.RsdCG, invParam.MaxCG);
	  res.n_count += n_CG ;
	  
	  snoop.start();
	  GoodEvecs.evec.AddVectors(evec,  subset());
	  GoodEvecs.eval.AddVectors(lambda,subset());
	  snoop.stop();
	  double Time = snoop.getTimeInSeconds() ;
	  QDPIO::cout<<"GoodEvecs.Neig= "<<GoodEvecs.Neig<<endl;
	  
	  snoop.start();
//????	  SimpleGramSchmidt(GoodEvecs.evec,GoodEvecs.Neig-invParam.Neig,GoodEvecs.Neig,subset());
	  SimpleGramSchmidt(GoodEvecs.evec.vec,GoodEvecs.Neig-invParam.Neig,GoodEvecs.Neig,subset());
	  snoop.stop();
	  Time += snoop.getTimeInSeconds() ;
	  
	  snoop.start();
	  LinAlg::Matrix<DComplex> Htmp(GoodEvecs.Neig) ;
	  SubSpaceMatrix(Htmp,*A,GoodEvecs.evec.vec,GoodEvecs.Neig);
	  //OctavePrintOut(Htmp.mat,Htmp.N,tag("H",i),"RayleighRich.m") ;
	  char V = 'V' ; char U = 'U' ;
	  Lapack::zheev(V,U,Htmp.mat,lambda);
	  evec.resize(GoodEvecs.Neig) ;
	  //OctavePrintOut(Htmp.mat,Htmp.N,tag("Htmp",i),"RayleighRich.m") ;
	  for(int k(0);k<GoodEvecs.Neig;k++){
	    GoodEvecs.eval[k] = lambda[k];
	    GoodEvecs.evec[k][subset()] = zero ;
	    //cout<<k<<endl ;
	    for(int j(0);j<GoodEvecs.Neig;j++)
	      evec[k][subset()] += conj(Htmp(k,j))*GoodEvecs.evec[j] ;
	    //QDPIO::cout<<"norm(evec["<<k<<"])=";
	    //QDPIO::cout<< sqrt(norm2(GoodEvecs.evec[k],subset()))<<endl;
	  }
	  snoop.stop();
	  Time += snoop.getTimeInSeconds() ;
	  QDPIO::cout << "Evec_Refinement: time = "
		      << Time
		      << " secs" << endl;
	  
	  QDPIO::cout<<"GoodEvecs.Neig= "<<GoodEvecs.Neig<<endl;
	  for(int k(0);k<GoodEvecs.Neig;k++)
	    GoodEvecs.evec[k][subset()]  = evec[k] ;
	  
	  //Check the quality of eigenvectors
#ifdef TEST_ALGORITHM
	  {
	    T Av ;
	    for(int k(0);k<GoodEvecs.Neig;k++)
	    {
	      (A*)(Av,GoodEvecs.evec[k],PLUS) ;
	      DComplex rq = innerProduct(GoodEvec[k],Av,subset());
	      Av[subset()] -= GoodEvecs.eval[k]*GoodEvecs.evec[k] ;
	      Double tt = sqrt(norm2(Av,subset()));
	      QDPIO::cout<<"REFINE: error evec["<<k<<"] = "<<tt<<" " ;
	      QDPIO::cout<<"--- eval ="<<GoodEvecs.eval[k]<<" ";
	      tt =  sqrt(norm2(GoodEvecs.evec[k],subset()));
	      QDPIO::cout<<"--- rq ="<<real(rq)<<" ";
	      QDPIO::cout<<"--- norm = "<<tt<<endl  ;
	    } 
	  }
#endif
	}// if there is space
	else // call CG but ask it not to compute vectors
	{
	  evec.resize(0);
	  n_CG = res.n_count ;
	  if(invParam.vPrecCGvecs ==0){
	    res =InvEigCG2(*A, 
			   psi,
			   chi_tmp,
			   lambda, 
			   evec, 
			   0, //Eigenvectors to keep
			   invParam.Nmax,  // Max vectors to work with
			   invParam.RsdCGRestart, // CG residual.... Empirical restart need a param here...
			   invParam.MaxCG // Max CG itterations
			   ); 
	  }
	  else
	    res= vecPrecondCG(*A,psi,chi_tmp,GoodEvecs.eval.vec,GoodEvecs.evec.vec,
			      invParam.vPrecCGvecStart,
			      invParam.vPrecCGvecStart+invParam.vPrecCGvecs,
                              invParam.RsdCG, // CG residual       
                              invParam.MaxCG // Max CG itterations    
			      );
	  res.n_count += n_CG ;
	  n_CG=0 ;
	  QDPIO::cout<<"Restart1: "<<endl ;
	  InitGuess(*A,psi,chi_tmp,GoodEvecs.eval.vec,GoodEvecs.evec.vec,GoodEvecs.Neig, n_CG);
	  res.n_count += n_CG ;
	  n_CG = res.n_count ;
	  if(invParam.vPrecCGvecs ==0)
	    res =InvEigCG2(*A,
                           psi,
                           chi_tmp,
                           lambda,
                           evec,
                           0, //Eigenvectors to keep                             
                           invParam.Nmax,  // Max vectors to work with           
                           invParam.RsdCG, // CG residual....
                           invParam.MaxCG // Max CG itterations             
                           );
	  else
	    res= vecPrecondCG(*A,psi,chi_tmp,GoodEvecs.eval.vec,GoodEvecs.evec.vec,
			      invParam.vPrecCGvecStart,
                              invParam.vPrecCGvecStart+invParam.vPrecCGvecs,
			      invParam.RsdCG, // CG residual             
			      invParam.MaxCG // Max CG itterations          
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
  };

} // End namespace

#endif 

