// $Id: syssolver_mdagm_eigcg_qdp.cc,v 3.2 2009-04-17 02:05:32 bjoo Exp $
/*! \file
 *  \brief Solve a M^dag*M*psi=chi linear system by EigCG
 */

#include <qdp-lapack.h>

#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"

#include "actions/ferm/invert/syssolver_mdagm_eigcg_qdp.h"
#include "actions/ferm/invert/inv_eigcg2.h"
#include "actions/ferm/invert/norm_gram_schm.h"
#include "actions/ferm/invert/invcg2.h"

//for debugging
//#include "octave.h"
#define TEST_ALGORITHM
namespace Chroma
{

  //! CG1 system solver namespace
  namespace MdagMSysSolverQDPEigCGEnv
  {
    //! Callback function
    MdagMSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 

						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new MdagMSysSolverQDPEigCG<LatticeFermion>(A, SysSolverEigCGParams(xml_in, path));
    }

#if 0
    //! Callback function
    MdagMSystemSolver<LatticeStaggeredFermion>* createStagFerm(
      XMLReader& xml_in,
      const std::string& path,
      Handle< LinearOperator<LatticeStaggeredFermion> > A)
    {
      return new MdagMSysSolverQDPEigCG<LatticeStaggeredFermion>(A, SysSolverEigCGParams(xml_in, path));
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
	success &= Chroma::TheMdagMFermSystemSolverFactory::Instance().registerObject(name, createFerm);
//	success &= Chroma::TheMdagMStagFermSystemSolverFactory::Instance().registerObject(name, createStagFerm);
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

      LinAlg::RitzPairs<T>& GoodEvecs = TheNamedObjMap::Instance().getData< LinAlg::RitzPairs<T> >(invParam.eigen_id);

      multi1d<Double> lambda ; //the eigenvalues
      multi1d<T> evec(0); // The eigenvectors  
      SystemSolverResults_t res;  // initialized by a constructor
      int n_CG(0);
      int flag(-1);//first time through
      int restart(0);
      Real restartTol = invParam.restartTol ;
      StopWatch snoop;
      while((flag==-1)||flag==3){
	flag=0 ;
	if(invParam.PrintLevel>0)
	  QDPIO::cout<<"GoodEvecs.Neig= "<<GoodEvecs.Neig<<endl;
	if(GoodEvecs.Neig>0){//deflate if there avectors to deflate
	  if(invParam.PrintLevel>0){
	          snoop.reset();
		  snoop.start();
	  }
	  InvEigCG2Env::InitGuess(MdagM,psi,chi,GoodEvecs.eval.vec,GoodEvecs.evec.vec,GoodEvecs.Neig,n_CG);
	  if(invParam.PrintLevel>0) snoop.stop();
	  if(invParam.PrintLevel>0)
	    QDPIO::cout << "InitGuess:  time = "
			<< snoop.getTimeInSeconds() 
			<< " secs" << endl;
	}
	//if there is space for new
	if((GoodEvecs.Neig)<GoodEvecs.evec.vec.size())
	  {

	    evec.resize(0);//get in there with no evecs so that it computes new
	    res = InvEigCG2Env::InvEigCG2(MdagM, psi, chi, lambda, evec, 
					  invParam.Neig, invParam.Nmax, 
					  invParam.RsdCG, invParam.MaxCG,
					  invParam.PrintLevel);
	    res.n_count += n_CG ;

	    snoop.reset();
	    snoop.start();
	    GoodEvecs.AddVectors(lambda, evec, MdagM.subset());
	    snoop.stop();
	    double Time = snoop.getTimeInSeconds() ;
	  
	    snoop.start();
	    normGramSchmidt(GoodEvecs.evec.vec,GoodEvecs.Neig-invParam.Neig,GoodEvecs.Neig,MdagM.subset());
	    normGramSchmidt(GoodEvecs.evec.vec,GoodEvecs.Neig-invParam.Neig,GoodEvecs.Neig,MdagM.subset());
	    snoop.stop();
	    Time += snoop.getTimeInSeconds() ;
	  
	    snoop.start();
	    LinAlg::Matrix<DComplex> Htmp(GoodEvecs.Neig) ;
	    //Only do the matrix elements for the new vectors
	    //InvEigCG2Env::SubSpaceMatrix(Htmp,MdagM,GoodEvecs.evec.vec,GoodEvecs.Neig);
	    /**/
	    InvEigCG2Env::SubSpaceMatrix(Htmp,MdagM,
					 GoodEvecs.evec.vec,
					 GoodEvecs.eval.vec,
					 GoodEvecs.Neig,
					 GoodEvecs.Neig-invParam.Neig);
	    /**/
	    char V = 'V' ; char U = 'U' ;
	    QDPLapack::zheev(V,U,Htmp.mat,lambda);
	    evec.resize(GoodEvecs.Neig) ;

	    for(int k(0);k<GoodEvecs.Neig;k++){
	      GoodEvecs.eval[k] = lambda[k];
	      evec[k][MdagM.subset()] = zero ;
	      for(int j(0);j<GoodEvecs.Neig;j++)
		evec[k][MdagM.subset()] += conj(Htmp(k,j))*GoodEvecs.evec[j] ;
	    }
	    snoop.stop();
	    Time += snoop.getTimeInSeconds() ;
	    if(invParam.PrintLevel>0){
	      QDPIO::cout << "Evec_Refinement: time = "
			  << Time
			  << " secs" << endl;
	      
	      //QDPIO::cout<<"GoodEvecs.Neig= "<<GoodEvecs.Neig<<endl;
	    }
	    for(int k(0);k<GoodEvecs.Neig;k++)
	      GoodEvecs.evec[k][MdagM.subset()]  = evec[k] ;
	  
	    //Check the quality of eigenvectors
	    if(invParam.PrintLevel>4){
	      T Av ;
	      for(int k(0);k<GoodEvecs.Neig;k++){
		MdagM(Av,GoodEvecs.evec[k],PLUS) ;
		DComplex rq = innerProduct(GoodEvecs.evec[k],Av,MdagM.subset());
		Av[MdagM.subset()] -= GoodEvecs.eval[k]*GoodEvecs.evec[k] ;
		Double tt = sqrt(norm2(Av,MdagM.subset()));
		QDPIO::cout<<"REFINE: error evec["<<k<<"] = "<<tt<<" " ;
		QDPIO::cout<<"--- eval ="<<GoodEvecs.eval[k]<<" ";
		tt =  sqrt(norm2(GoodEvecs.evec[k],MdagM.subset()));
		QDPIO::cout<<"--- rq ="<<real(rq)<<" ";
		QDPIO::cout<<"--- norm = "<<tt<<endl  ;
	      } 
	    }
	  }// if there is space
	else // call CG but ask it not to compute vectors
	  {
	    evec.resize(0);
	    n_CG = res.n_count ;
	    if(invParam.vPrecCGvecs ==0){
	      if(invParam.PrintLevel<2)// Call the CHROMA CG 
		res = InvCG2(A, chi, psi, restartTol, invParam.MaxCG);
	      else
		res = InvEigCG2Env::InvEigCG2(MdagM, 
					      psi,
					      chi,
					      lambda, 
					      evec, 
					      0, //Eigenvectors to keep
					      invParam.Nmax,  // Max vectors to work with
					      restartTol, // CG residual...
					      invParam.MaxCG, // Max CG itterations
					      invParam.PrintLevel
					      );
	    }
	    else
	      res = InvEigCG2Env::vecPrecondCG(MdagM,psi,chi,
					       GoodEvecs.eval.vec,
					       GoodEvecs.evec.vec,
					       invParam.vPrecCGvecStart,
					       invParam.vPrecCGvecStart+invParam.vPrecCGvecs,
					       restartTol, // CG residual      
					       invParam.MaxCG // Max CG itterations    
					       );
	    res.n_count += n_CG ;
	    if(toBool(restartTol!=invParam.RsdCG)){
	      restart++;//count the number of restarts
	      if(invParam.PrintLevel>0)
		QDPIO::cout<<"Restart: "<<restart<<endl ;
	      flag=3 ; //restart
	    }
	    else{
	      flag=0; //stop restarting
	    }
	    restartTol *=restartTol ;
	    if(toBool(restartTol < invParam.RsdCG)){
	      restartTol = invParam.RsdCG;
	    }
	  }
      }//while

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
  MdagMSysSolverQDPEigCG<LatticeFermionF>::operator()(LatticeFermionF& psi, const LatticeFermionF& chi) const
  {
    return sysSolver(psi, chi, *A, *MdagM, invParam);
  }

  // LatticeFermionD
  template<>
  SystemSolverResults_t
  MdagMSysSolverQDPEigCG<LatticeFermionD>::operator()(LatticeFermionD& psi, const LatticeFermionD& chi) const
  {
    return sysSolver(psi, chi, *A, *MdagM, invParam);
  }

#if 0
  // Not quite ready yet for these - almost there
  // LatticeStaggeredFermionF
  template<>
  SystemSolverResults_t
  MdagMSysSolverQDPEigCG<LatticeStaggeredFermionF>::operator()(LatticeStaggeredFermionF& psi, const LatticeStaggeredFermionF& chi) const
  {
    return sysSolver(psi, chi, *A, *MdagM, invParam);
  }

  // LatticeStaggeredFermionD
  template<>
  SystemSolverResults_t
  MdagMSysSolverQDPEigCG<LatticeStaggeredFermionD>::operator()(LatticeStaggeredFermionD& psi, const LatticeStaggeredFermionD& chi) const
  {
    return sysSolver(psi, chi, *A, *MdagM, invParam);
  }
#endif

}
