// $Id: syssolver_linop_eigcg.cc,v 1.5 2007-10-24 01:50:15 edwards Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"

#include "actions/ferm/invert/syssolver_linop_eigcg.h"

namespace Chroma
{

  //! CG1 system solver namespace
  namespace LinOpSysSolverEigCGEnv
  {
    //! Callback function
    LinOpSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new LinOpSysSolverEigCG<LatticeFermion>(A, SysSolverEigCGParams(xml_in, path));
    }

#if 0
    //! Callback function
    LinOpSystemSolver<LatticeStaggeredFermion>* createStagFerm(
      XMLReader& xml_in,
      const std::string& path,
      Handle< LinearOperator<LatticeStaggeredFermion> > A)
    {
      return new LinOpSysSolverEigCG<LatticeStaggeredFermion>(A, SysSolverEigCGParams(xml_in, path));
    }
#endif

    //! Name to be used
    const std::string name("EIG_CG_INVERTER");

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheLinOpFermSystemSolverFactory::Instance().registerObject(name, createFerm);
//	success &= Chroma::TheLinOpStagFermSystemSolverFactory::Instance().registerObject(name, createStagFerm);
	registered = true;
      }
      return success;
    }
  }


  //! Anonymous namespace holding method
  namespace
  {
    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     */
    template<typename T>
    SystemSolverResults_t sysSolver(T& psi, const T& chi, 
				    const LinearOperator<T>& A, 
				    const LinearOperator<T>& MdagM, 
				    const SysSolverEigCGParams& invParam)
    {
      START_CODE();

      T chi_tmp;
      A(chi_tmp, chi, MINUS);

      LinAlg::RitzPairs<T>& GoodEvecs = TheNamedObjMap::Instance().getData< LinAlg::RitzPairs<T> >(invParam.eigen_id);

      multi1d<Double> lambda ; //the eigenvalues
      multi1d<T> evec(0); // The eigenvectors  
      SystemSolverResults_t res;  // initialized by a constructor
      int n_CG(0);
      // Need to pass the appropriate operators here...
      InvEigCG2Env::InitGuess(MdagM,psi,chi_tmp,GoodEvecs.eval.vec,GoodEvecs.evec.vec,GoodEvecs.Neig,n_CG);
      if((GoodEvecs.Neig+invParam.Neig)<=invParam.Neig_max) // if there is space for new
      {
	StopWatch snoop;
	snoop.reset();
	snoop.start();

	evec.resize(0);// get in there with no evecs so that it computes new
	res = InvEigCG2Env::InvEigCG2(MdagM, psi, chi_tmp, lambda, evec, 
				      invParam.Neig, invParam.Nmax, 
				      invParam.RsdCG, invParam.MaxCG);
	res.n_count += n_CG ;
	  
	snoop.start();
	GoodEvecs.evec.AddVectors(evec,  MdagM.subset());
	GoodEvecs.eval.AddVectors(lambda,MdagM.subset());
	snoop.stop();
	double Time = snoop.getTimeInSeconds() ;
	QDPIO::cout<<"GoodEvecs.Neig= "<<GoodEvecs.Neig<<endl;
	  
	snoop.start();
//????	  normGramSchmidt(GoodEvecs.evec,GoodEvecs.Neig-invParam.Neig,GoodEvecs.Neig,MdagM.subset());
	normGramSchmidt(GoodEvecs.evec.vec,GoodEvecs.Neig-invParam.Neig,GoodEvecs.Neig,MdagM.subset());
	snoop.stop();
	Time += snoop.getTimeInSeconds() ;
	  
	snoop.start();
	LinAlg::Matrix<DComplex> Htmp(GoodEvecs.Neig) ;
	InvEigCG2Env::SubSpaceMatrix(Htmp,MdagM,GoodEvecs.evec.vec,GoodEvecs.Neig);
	//OctavePrintOut(Htmp.mat,Htmp.N,tag("H",i),"RayleighRich.m") ;
	char V = 'V' ; char U = 'U' ;
	Lapack::zheev(V,U,Htmp.mat,lambda);
	evec.resize(GoodEvecs.Neig) ;
	//OctavePrintOut(Htmp.mat,Htmp.N,tag("Htmp",i),"RayleighRich.m") ;
	for(int k(0);k<GoodEvecs.Neig;k++){
	  GoodEvecs.eval[k] = lambda[k];
	  GoodEvecs.evec[k][MdagM.subset()] = zero ;
	  //cout<<k<<endl ;
	  for(int j(0);j<GoodEvecs.Neig;j++)
	    evec[k][MdagM.subset()] += conj(Htmp(k,j))*GoodEvecs.evec[j] ;
	  //QDPIO::cout<<"norm(evec["<<k<<"])=";
	  //QDPIO::cout<< sqrt(norm2(GoodEvecs.evec[k],MdagM.subset()))<<endl;
	}
	snoop.stop();
	Time += snoop.getTimeInSeconds() ;
	QDPIO::cout << "Evec_Refinement: time = "
		    << Time
		    << " secs" << endl;
	  
	QDPIO::cout<<"GoodEvecs.Neig= "<<GoodEvecs.Neig<<endl;
	for(int k(0);k<GoodEvecs.Neig;k++)
	  GoodEvecs.evec[k][MdagM.subset()]  = evec[k] ;
	  
	//Check the quality of eigenvectors
#ifdef TEST_ALGORITHM
	{
	  T Av ;
	  for(int k(0);k<GoodEvecs.Neig;k++)
	  {
	    (A*)(Av,GoodEvecs.evec[k],PLUS) ;
	    DComplex rq = innerProduct(GoodEvec[k],Av,MdagM.subset());
	    Av[MdagM.subset()] -= GoodEvecs.eval[k]*GoodEvecs.evec[k] ;
	    Double tt = sqrt(norm2(Av,MdagM.subset()));
	    QDPIO::cout<<"REFINE: error evec["<<k<<"] = "<<tt<<" " ;
	    QDPIO::cout<<"--- eval ="<<GoodEvecs.eval[k]<<" ";
	    tt =  sqrt(norm2(GoodEvecs.evec[k],MdagM.subset()));
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
	  res = InvEigCG2Env::InvEigCG2(MdagM, 
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
	  res = InvEigCG2Env::vecPrecondCG(MdagM,psi,chi_tmp,GoodEvecs.eval.vec,GoodEvecs.evec.vec,
					   invParam.vPrecCGvecStart,
					   invParam.vPrecCGvecStart+invParam.vPrecCGvecs,
					   invParam.RsdCG, // CG residual       
					   invParam.MaxCG // Max CG itterations    
	    );
	res.n_count += n_CG ;
	n_CG=0 ;
	QDPIO::cout<<"Restart1: "<<endl ;
	InvEigCG2Env::InitGuess(MdagM,psi,chi_tmp,GoodEvecs.eval.vec,GoodEvecs.evec.vec,GoodEvecs.Neig, n_CG);
	res.n_count += n_CG ;
	n_CG = res.n_count ;
	if(invParam.vPrecCGvecs ==0)
	  res = InvEigCG2Env::InvEigCG2(MdagM,
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
	  res= InvEigCG2Env::vecPrecondCG(MdagM,psi,chi_tmp,GoodEvecs.eval.vec,GoodEvecs.evec.vec,
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

  } // anonymous namespace


  //
  // Wrappers
  //

  // LatticeFermionF
  template<>
  SystemSolverResults_t
  LinOpSysSolverEigCG<LatticeFermionF>::operator()(LatticeFermionF& psi, const LatticeFermionF& chi) const
  {
    return sysSolver(psi, chi, *A, *MdagM, invParam);
  }

  // LatticeFermionD
  template<>
  SystemSolverResults_t
  LinOpSysSolverEigCG<LatticeFermionD>::operator()(LatticeFermionD& psi, const LatticeFermionD& chi) const
  {
    return sysSolver(psi, chi, *A, *MdagM, invParam);
  }

#if 0
  // Not quite ready yet for these - almost there
  // LatticeStaggeredFermionF
  template<>
  SystemSolverResults_t
  LinOpSysSolverEigCG<LatticeStaggeredFermionF>::operator()(LatticeStaggeredFermionF& psi, const LatticeStaggeredFermionF& chi) const
  {
    return sysSolver(psi, chi, *A, *MdagM, invParam);
  }

  // LatticeStaggeredFermionD
  template<>
  SystemSolverResults_t
  LinOpSysSolverEigCG<LatticeStaggeredFermionD>::operator()(LatticeStaggeredFermionD& psi, const LatticeStaggeredFermionD& chi) const
  {
    return sysSolver(psi, chi, *A, *MdagM, invParam);
  }
#endif

}
