// $Id: syssolver_mdagm_OPTeigcg.cc,v 3.4 2009-07-31 14:09:31 bjoo Exp $
/*! \file
 *  \brief Solve a M^dag*M*psi=chi linear system by EigCG
 */


#include "qdp-lapack_Complex.h"
#include "qdp-lapack_eigpcg.h"  
#include "qdp-lapack_IncrEigpcg.h"

#include "actions/ferm/invert/syssolver_mdagm_factory.h"
#include "actions/ferm/invert/syssolver_mdagm_aggregate.h"

#include "actions/ferm/invert/syssolver_mdagm_OPTeigcg.h"
#include "containers.h"


namespace Chroma
{

  //! CG1 system solver namespace
  namespace MdagMSysSolverOptEigCGEnv
  {
    //! Callback function
    MdagMSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new MdagMSysSolverOptEigCG<LatticeFermion>(A, SysSolverOptEigCGParams(xml_in, path));
    }

#if 0
    //! Callback function
    MdagMSystemSolver<LatticeStaggeredFermion>* createStagFerm(
      XMLReader& xml_in,
      const std::string& path,
      Handle< LinearOperator<LatticeStaggeredFermion> > A)
    {
      return new MdagMSysSolverOptEigCG<LatticeStaggeredFermion>(A, SysSolverOptEigCGParams(xml_in, path));
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
     template<typename T>
     struct MatVecArg{
       T XX ;
       T YY ;
       Handle< LinearOperator<T> > MdagM;
    };

    template<typename T>
    void MatrixMatvec(void *x, void *y, void *params) {
      
      MatVecArg<T> &arg = *((MatVecArg<T> *) params) ;
      
      //Works only in single precision CHROMA
      RComplex<float> *px = (RComplex<float> *) x;
      RComplex<float> *py = (RComplex<float> *) y;

      //XX.getF() = (Complex *) x ; //AN DOULEPSEI AUTO NA ME FTUSEIS
      //YY.getF() = (Complex *) y ; //AN DOULEPSEI AUTO NA ME FTUSEIS

      //Alliws kanoume copy
      //copy x into XX
      Subset s = arg.MdagM->subset() ; 
      if(s.hasOrderedRep()){
	/**
	int one(1);
        int VecSize = s.numSiteTable()*Nc*Ns ;
	BLAS_DCOPY(&VecSize, 
	       (double *)&arg.XX.elem(s.start()).elem(0).elem(0).real(),
	       &one,
	       (double *)x, &one);
	**/
	/**/
	int count=0 ;
	//can be done with ccopy for speed...
	for(int i=s.start(); i <= s.end(); i++)
	  for(int ss(0);ss<Ns;ss++)
	    for(int c(0);c<Nc;c++){
	      arg.XX.elem(i).elem(ss).elem(c)  = *(px+count);
	      count++;
	    }
	/**/
      }
      else{
	int i ;
	const int *tab = s.siteTable().slice();
	int count=0;
	for(int x=0; x < s.numSiteTable(); ++x){
	  i = tab[x] ;
	  for(int ss(0);ss<Ns;ss++)
	    for(int c(0);c<Nc;c++){
	      arg.XX.elem(i).elem(ss).elem(c) = *(px+count);
	      count++;
	    }
	}
      }
      

      (*arg.MdagM)(arg.YY,arg.XX,PLUS) ;
      //T foo,boo;
      //*(arg.MdagM)(boo,foo,PLUS) ;

      //copy back..
      if(s.hasOrderedRep()){
	int count=0 ;
	//can be done with ccopy for speed...
	for(int i=s.start(); i <= s.end(); i++)
	  for(int ss(0);ss<Ns;ss++)
	    for(int c(0);c<Nc;c++){
	      *(py+count) = arg.YY.elem(i).elem(ss).elem(c) ;
	      count++;
	    }
      }
      else{
	int i ;
	int count=0;
	const int *tab = s.siteTable().slice();
	for(int x=0; x < s.numSiteTable(); ++x){
	  i = tab[x] ;
	  for(int ss(0);ss<Ns;ss++)
	    for(int c(0);c<Nc;c++){
	      *(py+count) = arg.YY.elem(i).elem(ss).elem(c) ;
	      count++;
	    }
	}
      }
    }
    

    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     */
    template<typename T>
    SystemSolverResults_t sysSolver(T& psi, const T& chi, 
				    const LinearOperator<T>& A, 
				    Handle< LinearOperator<T> > MdagM, 
				    const SysSolverOptEigCGParams& invParam)
    {
#ifndef QDP_IS_QDPJIT
      START_CODE();

      SystemSolverResults_t res;  // initialized by a constructor

      LinAlg::OptEigInfo& EigInfo = TheNamedObjMap::Instance().getData< LinAlg::OptEigInfo >(invParam.eigen_id);

      QDPIO::cout<<"EigInfo.N= "<<EigInfo.N<<endl ;
      QDPIO::cout<<"EigInfo.lde= "<<EigInfo.lde<<endl ;
      QDPIO::cout<<"EigInfo.ldh= "<<EigInfo.evals.size()<<endl ;
      QDPIO::cout<<"EigInfo.ncurEvals= "<<EigInfo.ncurEvals<<endl ;
      QDPIO::cout<<"EigInfo.restartTol= "<<EigInfo.restartTol<<endl ;

      Subset s = A.subset() ;

      Complex_C *X ; 
      Complex_C *B ; 
      Complex_C *work=NULL  ;
      Complex_C *V=NULL     ;
      Complex_C *ework=NULL ;

      // OK, for now Weirdly enough psi and chi have to be
      // explicit single precision, so downcast those...
      // Most generically the single prec type of say a 
      // type T should be given by the QDP++ Type trait:
      // SinglePrecType<T>::Type_t
      typename SinglePrecType<T>::Type_t psif = psi;
      typename SinglePrecType<T>::Type_t chif = chi;

      if(s.hasOrderedRep()){
	X = (Complex_C *) &psif.elem(s.start()).elem(0).elem(0).real();
	B = (Complex_C *) &chif.elem(s.start()).elem(0).elem(0).real();
      }
      else{//need to copy
	//X = allocate space for them
	//B =  allocate space for them...
	QDPIO::cout<<"OPPS! I have not implemented OPT_EigCG for Linops with non contigius subset\n";
	exit(1);
      }

      Complex_C *evecs = (Complex_C *) &EigInfo.evecs[0] ;
      float *evals = (float *) &EigInfo.evals[0].elem() ;
      Complex_C *H  = (Complex_C *) &EigInfo.H[0] ;
      Complex_C *HU = (Complex_C *) &EigInfo.HU[0] ;
      MatVecArg<T> arg ;
      arg.MdagM = MdagM ;
      int esize = invParam.esize*Layout::sitesOnNode()*Nc*Ns ;

      QDPIO::cout<<"OPT_EIGCG_SYSSOLVER= "<<esize<<endl ;
      //multi1d<Complex_C> ework(esize);
      float resid = (float) invParam.RsdCG.elem().elem().elem().elem();
      float AnormEst = invParam.NormAest.elem().elem().elem().elem();

      float restartTol;
      if (EigInfo.ncurEvals < EigInfo.evals.size()) 
         restartTol = 0.0; 		  // Do not restart the first phase 
      else {
        // set it to the user parameter:
	restartTol = invParam.restartTol.elem().elem().elem().elem();

        // restartTol = EigInfo.restartTol; //or restart with tol as computed 
      }

      IncrEigpcg(EigInfo.N, EigInfo.lde, 1, X, B, 
		 &EigInfo.ncurEvals, EigInfo.evals.size(), 
		 evecs, evals, H, HU, 
		 MatrixMatvec<T>, NULL, (void *)&arg, work, V, 
		 ework, esize, 
		 resid, &restartTol,
		 AnormEst, invParam.updateRestartTol, 
		 invParam.MaxCG, invParam.PrintLevel, 
		 invParam.Neig, invParam.Nmax, stdout);

      /* Update the restartTol in the EigInfo function */
      EigInfo.restartTol = restartTol;

      // Copy result back
      psi = psif;

      T tt;
      (*MdagM)(tt,psi,PLUS);
      QDPIO::cout<<"OPT_EICG_SYSSOLVER: True residual after solution : "<<sqrt(norm2(tt-chi,s))<<endl ;
      QDPIO::cout<<"OPT_EICG_SYSSOLVER: norm of  solution            : "<<sqrt(norm2(psi,s))<<endl ;
      QDPIO::cout<<"OPT_EICG_SYSSOLVER: norm of rhs                  : "<<sqrt(norm2(chi,s))<<endl ;
      


      if(!s.hasOrderedRep()){
	QDPIO::cout<<"OPPS! I have no implemented OPT_EigCG for Linops with non contigius subset\n";
      }
      END_CODE();

      return res;
#endif
    }

  } // anonymous namespace


  //
  // Wrappers
  //

  // LatticeFermionF
  template<>
  SystemSolverResults_t
  MdagMSysSolverOptEigCG<LatticeFermionF>::operator()(LatticeFermionF& psi, const LatticeFermionF& chi) const
  {
    return sysSolver(psi, chi, *A, MdagM, invParam);
  }

  // LatticeFermionD
  template<>
  SystemSolverResults_t
  MdagMSysSolverOptEigCG<LatticeFermionD>::operator()(LatticeFermionD& psi, const LatticeFermionD& chi) const
  {
    return sysSolver(psi, chi, *A, MdagM, invParam);
  }

#if 1
 

  // Not quite ready yet for these - almost there
  // LatticeStaggeredFermionF
  template<>
  SystemSolverResults_t
  MdagMSysSolverOptEigCG<LatticeStaggeredFermionF>::operator()(LatticeStaggeredFermionF& psi, const LatticeStaggeredFermionF& chi) const
  {
    return sysSolver(psi, chi, *A, MdagM, invParam);
  }

  // LatticeStaggeredFermionD
  template<>
  SystemSolverResults_t
  MdagMSysSolverOptEigCG<LatticeStaggeredFermionD>::operator()(LatticeStaggeredFermionD& psi, const LatticeStaggeredFermionD& chi) const
  {
    return sysSolver(psi, chi, *A, MdagM, invParam);
  }
#endif

}
