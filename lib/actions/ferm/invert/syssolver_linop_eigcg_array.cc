// $Id: syssolver_linop_eigcg_array.cc,v 1.4 2008-04-09 04:49:23 kostas Exp $
/*! \file
 *  \brief Solve a M*psi=chi linear system array by EigCG2
 */

#include <qdp-lapack.h>

#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"

#include "actions/ferm/invert/syssolver_linop_eigcg_array.h"
#include "actions/ferm/invert/inv_eigcg2_array.h"
#include "actions/ferm/invert/norm_gram_schm.h"

//for debugging
//#include "octave.h"
#define TEST_ALGORITHM
namespace Chroma
{

  //! CG1 system solver namespace
  namespace LinOpSysSolverEigCGArrayEnv
  {
    //! Callback function
    LinOpSystemSolverArray<LatticeFermion>* createFerm(XMLReader& xml_in,
						       const std::string& path,
						       Handle< LinearOperatorArray<LatticeFermion> > A)
    {
      return new LinOpSysSolverEigCGArray<LatticeFermion>(A, SysSolverEigCGParams(xml_in, path));
    }

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
	success &= Chroma::TheLinOpFermSystemSolverArrayFactory::Instance().registerObject(name, createFerm);
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
    SystemSolverResults_t sysSolver(multi1d<T>& psi, const multi1d<T>& chi, 
				    const LinearOperatorArray<T>& A, 
				    const LinearOperatorArray<T>& MdagM, 
				    const SysSolverEigCGParams& invParam)
    {
      START_CODE();
#ifdef _WORKING_
      multi1d<T> chi_tmp;
      A(chi_tmp, chi, MINUS);

      int Ls = chi.size();

      LinAlg::RitzPairsArray<T>& GoodEvecs = 
	TheNamedObjMap::Instance().getData< LinAlg::RitzPairsArray<T> >(invParam.eigen_id);

      multi1d<Double> lambda; //the eigenvalues
      multi2d<T>      evec; // The eigenvectors  
      SystemSolverResults_t res;  // initialized by a constructor
      int n_CG(0);
      // Need to pass the appropriate operators here...
      InvEigCG2ArrayEnv::InitGuess(MdagM,psi,chi_tmp,GoodEvecs.eval.vec,GoodEvecs.evec.vec,GoodEvecs.Neig,n_CG);
      if((GoodEvecs.Neig+invParam.Neig)<=invParam.Neig_max) // if there is space for new
      {
	StopWatch snoop;
	snoop.reset();
	snoop.start();

	evec.resize(0,0);// get in there with no evecs so that it computes new
	res = InvEigCG2ArrayEnv::InvEigCG2(MdagM, psi, chi_tmp, lambda, evec, 
					   invParam.Neig, invParam.Nmax, 
					   invParam.RsdCG, invParam.MaxCG);
	res.n_count += n_CG;
	  
	snoop.start();
	//GoodEvecs.evec.AddVectors(evec,  MdagM.subset());
	//GoodEvecs.eval.AddVectors(lambda,MdagM.subset());
	GoodEvecs.AddVectors(lambda, evec, MdagM.subset());
	snoop.stop();
	double Time = snoop.getTimeInSeconds();
	QDPIO::cout<<"GoodEvecs.Neig= "<<GoodEvecs.Neig<<endl;
	  
	snoop.start();
        normGramSchmidt(GoodEvecs.evec.vec,GoodEvecs.Neig-invParam.Neig,GoodEvecs.Neig,MdagM.subset());
	normGramSchmidt(GoodEvecs.evec.vec,GoodEvecs.Neig-invParam.Neig,GoodEvecs.Neig,MdagM.subset());
	snoop.stop();
	Time += snoop.getTimeInSeconds();
	  
	snoop.start();
	LinAlg::Matrix<DComplex> Htmp(GoodEvecs.Neig);
	InvEigCG2ArrayEnv::SubSpaceMatrix(Htmp,MdagM,GoodEvecs.evec.vec,GoodEvecs.Neig);
	//Octave::PrintOut(Htmp.mat,Htmp.N,Octave::tag("H"),"RayleighRich.m");
	char V = 'V'; char U = 'U';
	QDPLapack::zheev(V,U,Htmp.mat,lambda);
	evec.resize(GoodEvecs.Neig,Ls);
	//Octave::PrintOut(Htmp.mat,Htmp.N,Octave::tag("Hevec"),"RayleighRich.m");
	//evec.resize(GoodEvecs.Neig);
	for(int k(0);k<GoodEvecs.Neig;k++){
	  GoodEvecs.eval[k] = lambda[k];
	  for(int s=0; s < Ls; ++s)
	    evec[k][s][MdagM.subset()] = zero;
	  for(int j(0);j<GoodEvecs.Neig;j++)
	    for(int s=0; s < Ls; ++s)
	      evec[k][s][MdagM.subset()] += conj(Htmp(k,j))*GoodEvecs.evec[j][s];
	  //QDPIO::cout<<"norm(evec["<<k<<"])=";
	  //QDPIO::cout<< sqrt(norm2(GoodEvecs.evec[k],MdagM.subset()))<<endl;
	}
	snoop.stop();
	Time += snoop.getTimeInSeconds();
	QDPIO::cout << "Evec_Refinement: time = "
		    << Time
		    << " secs" << endl;
	  
	QDPIO::cout<<"GoodEvecs.Neig= "<<GoodEvecs.Neig<<endl;
	for(int k(0);k<GoodEvecs.Neig;k++)
	  for(int s=0; s < Ls; ++s)
	    GoodEvecs.evec[k][s][MdagM.subset()]  = evec[k][s];
	  
	//Check the quality of eigenvectors
#ifdef TEST_ALGORITHM
	{
	  multi1d<T> Av;
	  for(int k(0);k<GoodEvecs.Neig;k++)
	  {
	    MdagM(Av,GoodEvecs.evec[k],PLUS);
	    DComplex rq = innerProduct(GoodEvecs.evec[k],Av,MdagM.subset());
	    for(int s=0; s < Ls; ++s)
	      Av[s][MdagM.subset()] -= GoodEvecs.eval[k]*GoodEvecs.evec[k][s];
	    Double tt = sqrt(norm2(Av,MdagM.subset()));
	    QDPIO::cout<<"REFINE: error evec["<<k<<"] = "<<tt<<" ";
	    QDPIO::cout<<"--- eval ="<<GoodEvecs.eval[k]<<" ";
	    tt =  sqrt(norm2(GoodEvecs.evec[k],MdagM.subset()));
	    QDPIO::cout<<"--- rq ="<<real(rq)<<" ";
	    QDPIO::cout<<"--- norm = "<<tt<<endl;
	  } 
	}
#endif
      }// if there is space
      else // call CG but ask it not to compute vectors
      {
	evec.resize(0,0);
	n_CG = res.n_count;
	if(invParam.vPrecCGvecs ==0){
	  res = InvEigCG2ArrayEnv::InvEigCG2(MdagM, 
					     psi,
					     chi_tmp,
					     lambda, 
					     evec, 
					     0, //Eigenvectors to keep
					     invParam.Nmax,  // Max vectors to work with
					     invParam.RsdCGRestart[0], // CG residual.... Empirical restart need a param here...
					     invParam.MaxCG // Max CG itterations
	    ); 
	}
	else
	  res = InvEigCG2ArrayEnv::vecPrecondCG(MdagM,psi,chi_tmp,GoodEvecs.eval.vec,GoodEvecs.evec.vec,
						invParam.vPrecCGvecStart,
						invParam.vPrecCGvecStart+invParam.vPrecCGvecs,
						invParam.RsdCG, // CG residual       
						invParam.MaxCG // Max CG itterations    
	    );
	res.n_count += n_CG;
	n_CG=0;
	QDPIO::cout<<"Restart1: "<<endl;
	InvEigCG2ArrayEnv::InitGuess(MdagM,psi,chi_tmp,GoodEvecs.eval.vec,GoodEvecs.evec.vec,GoodEvecs.Neig, n_CG);
	res.n_count += n_CG;
	n_CG = res.n_count;
	if(invParam.vPrecCGvecs ==0)
	  res = InvEigCG2ArrayEnv::InvEigCG2(MdagM,
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
	  res= InvEigCG2ArrayEnv::vecPrecondCG(MdagM,psi,chi_tmp,GoodEvecs.eval.vec,GoodEvecs.evec.vec,
					       invParam.vPrecCGvecStart,
					       invParam.vPrecCGvecStart+invParam.vPrecCGvecs,
					       invParam.RsdCG, // CG residual             
					       invParam.MaxCG // Max CG itterations          
	    );
	res.n_count += n_CG;
      }
#else
      SystemSolverResults_t res;  // initialized by a constructor
      QDPIO::cerr<<" OOOPS! EigCG does not work with DWF\n";
      exit(1);
#endif
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
  LinOpSysSolverEigCGArray<LatticeFermionF>::operator()(multi1d<LatticeFermionF>& psi, const multi1d<LatticeFermionF>& chi) const
  {
    return sysSolver(psi, chi, *A, *MdagM, invParam);
  }

  // LatticeFermionD
  template<>
  SystemSolverResults_t
  LinOpSysSolverEigCGArray<LatticeFermionD>::operator()(multi1d<LatticeFermionD>& psi, const multi1d<LatticeFermionD>& chi) const
  {
    return sysSolver(psi, chi, *A, *MdagM, invParam);
  }


}
