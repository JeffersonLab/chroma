// $Id: syssolver_linop_OPTeigbicg.cc,v 3.1 2009-11-02 21:34:08 kostas Exp $
/*! \file
 *  \brief Solve a A*psi=chi linear system by EigBiCG
 */


#include "qdp-lapack_Complex.h"
#include "qdp-lapack_eigbicg.h"  
#include "qdp-lapack_IncrEigbicg.h"

#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"

#include "actions/ferm/invert/syssolver_linop_OPTeigbicg.h"
#include "containers.h"


namespace Chroma
{

  //! CG1 system solver namespace
  namespace LinOpSysSolverOptEigBiCGEnv
  {
    //! Callback function
    LinOpSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new LinOpSysSolverOptEigBiCG<LatticeFermion>(A, SysSolverOptEigBiCGParams(xml_in, path));
    }

#if 0
    //! Callback function
    LinOpSystemSolver<LatticeStaggeredFermion>* createStagFerm(
      XMLReader& xml_in,
      const std::string& path,
      Handle< LinearOperator<LatticeStaggeredFermion> > A)
    {
      return new LinOpSysSolverOptEigBiCG<LatticeStaggeredFermion>(A, SysSolverOptEigBiCGParams(xml_in, path));
    }
#endif

    //! Name to be used
    const std::string name("EIG_BiCG_INVERTER");

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
     template<typename T>
     struct MatVecArg{
       T XX ;
       T YY ;
       Handle< LinearOperator<T> > LinOp;
    };

    //PLUS is the mutrix and MINUS is the dagger
    template<typename R, typename T, PlusMinus isign>
    void MatrixMatvec(void *x, void *y, void *params){
      
      MatVecArg<T> &arg = *((MatVecArg<T> *) params) ;
      
      //Works only in single precision CHROMA
      RComplex<R> *px = (RComplex<R> *) x;
      RComplex<R> *py = (RComplex<R> *) y;

      //copy x into XX
      LinAlg::CopyToLatFerm<T,R>(arg.XX,px,arg.LinOp->subset());
      (*arg.LinOp)(arg.YY,arg.XX,isign) ;

      //copy back..
      LinAlg::CopyFromLatFerm<T,R>(py,arg.YY,arg.LinOp->subset());
    }


    void IncrEigBiCG(int n, int lde,int nrhs, Complex_C* X, Complex_C* B, 
		     int* ncurEvals, int ldh, Complex_C* evecsl, 
		     Complex_C* evecsr, Complex_C* evals,       
		     Complex_C* H,void (*matvec) (void*, void*, void*), 
		     void (*mathvec)(void*, void*, void*),void* params, 
		     float* AnormEst, Complex_C* work, Complex_C* VL,
		     int ldvl,Complex_C* VR, int ldvr,       
		     Complex_C* ework, int esize,float tol,float* restartTol,  
		     int maxit, char SRT_OPT, float epsi, int ConvTestOpt,
		     int plvl,int nev, int v_max,FILE* outputFile){
      QDPIO::cout<<"IncrEigbicg_C will be called"<<endl ;
      return IncrEigbicg_C(n, lde, nrhs, X, B, 
			   ncurEvals, ldh, 
			   evecsl, evecsr, evals, H,
			   matvec,mathvec,
			   params, AnormEst,
			   work, VL, ldvl, VR, 
			   ldvr, ework, esize, 
			   tol, restartTol,   
			   maxit, SRT_OPT, epsi, 
			   ConvTestOpt, 
			   plvl, nev, v_max, 
			   outputFile);
    }

    void IncrEigBiCG(int n, int lde,int nrhs, Complex_Z* X, Complex_Z* B, 
		     int* ncurEvals, int ldh, Complex_Z* evecsl, 
		     Complex_Z* evecsr, Complex_Z* evals,      
		     Complex_Z* H, void (*matvec) (void*, void*, void*), 
		     void (*mathvec)(void*, void*, void*), void* params, 
		     double *AnormEst, Complex_Z* work, Complex_Z* VL, 
		     int ldvl, Complex_Z* VR, int ldvr, Complex_Z* ework, 
		     int esize, double tol, double* restartTol,  
		     int maxit, char SRT_OPT, double epsi, int ConvTestOpt,
		     int plvl, int nev, int v_max,FILE *outputFile){
      QDPIO::cout<<"IncrEigbicg_Z will be called"<<endl ;
      return IncrEigbicg_Z(n, lde, nrhs, X, B, 
			   ncurEvals, ldh, 
			   evecsl, evecsr, evals, H,
			   matvec,mathvec,
			   params, AnormEst,
			   work, VL, ldvl, VR, 
			   ldvr, ework, esize, 
			   tol, restartTol,   
			   maxit, SRT_OPT, epsi, 
			   ConvTestOpt, 
			   plvl, nev, v_max, 
			   outputFile);
  }

    //! Solver the linear system
    /*!
     * \param psi      solution ( Modify )
     * \param chi      source ( Read )
     * \return syssolver results
     */
    template<typename T, typename F, typename C>
    SystemSolverResults_t sysSolver(T& psi, const T& chi, 
				    Handle< LinearOperator<T> > A, 
				    const SysSolverOptEigBiCGParams& invParam)
    {
      START_CODE();

      SystemSolverResults_t res;  // initialized by a constructor

      LinAlg::OptEigBiInfo<REAL>& EigBiInfo = TheNamedObjMap::Instance().getData< LinAlg::OptEigBiInfo<REAL> >(invParam.eigen_id);

      QDPIO::cout<<"EigBiInfo.N= "<<EigBiInfo.N<<endl ;
      QDPIO::cout<<"EigBiInfo.lde= "<<EigBiInfo.lde<<endl ;
      QDPIO::cout<<"EigBiInfo.ldh= "<<EigBiInfo.evals.size()<<endl ;
      QDPIO::cout<<"EigBiInfo.ncurEvals= "<<EigBiInfo.ncurEvals<<endl ;
      QDPIO::cout<<"EigBiInfo.restartTol= "<<EigBiInfo.restartTol<<endl ;

      Subset s = A->subset() ;

      C *X ; 
      C *B ; 
      C *work=NULL  ;
      C *VL=NULL    ;
      C *VR=NULL    ;
      C *ework=NULL ;

      /***** NO NEED FOR THIS ANY MORE
      // OK, for now Weirdly enough psi and chi have to be
      // explicit single precision, so downcast those...
      // Most generically the single prec type of say a 
      // type T should be given by the QDP++ Type trait:
      // SinglePrecType<T>::Type_t
      
      typename SinglePrecType<T>::Type_t psif = psi;
      typename SinglePrecType<T>::Type_t chif = chi;
      *****/

      if(s.hasOrderedRep()){
	X = (C *) &psi.elem(s.start()).elem(0).elem(0).real();
	B = (C *) &chi.elem(s.start()).elem(0).elem(0).real();
      }
      else{//need to copy
	//X = allocate space for them
	//B =  allocate space for them...
	QDPIO::cout<<"OPPS! I have not implemented OPT_EigBiCG for Linops with non contigius subset\n";
	exit(1);
      }

      C *evecsL = (C *) &EigBiInfo.evecsL[0] ;
      C *evecsR = (C *) &EigBiInfo.evecsR[0] ;
      C *evals      = (C *) &EigBiInfo.evals[0] ;
      C *H      = (C *) &EigBiInfo.H[0] ;
      MatVecArg<T> arg ;
      arg.LinOp = A ;
      int esize = invParam.esize*Layout::sitesOnNode()*Nc*Ns ;

      QDPIO::cout<<"OPT_EIGBICG_SYSSOLVER= "<<esize<<endl ;
      //multi1d<C> ework(esize);
      F resid = (F) invParam.RsdCG.elem().elem().elem().elem();
      F AnormEst = invParam.NormAest.elem().elem().elem().elem();

      F restartTol;
      if (EigBiInfo.ncurEvals < EigBiInfo.evals.size()) 
         restartTol = 0.0; 		  // Do not restart the first phase 
      else {
        // set it to the user parameter:
	restartTol = invParam.restartTol.elem().elem().elem().elem();

        // restartTol = EigInfo.restartTol; //or restart with tol as computed 
      }
      F epsi = invParam.epsi.elem().elem().elem().elem() ;

      IncrEigBiCG(EigBiInfo.N, EigBiInfo.lde, 1, X, B, 
		  &EigBiInfo.ncurEvals, EigBiInfo.evals.size(), 
		  evecsL, evecsR, evals, H, 
		  MatrixMatvec<F,T,PLUS>,  // mat vec
		  MatrixMatvec<F,T,MINUS> ,// hermitian conjugate mat vec
		  (void *)&arg, 
		  &AnormEst,
		  work, 
		  VL, EigBiInfo.lde, VR, EigBiInfo.lde,
		  ework, esize, 
		  resid, &restartTol,
		  invParam.MaxCG, 
		  invParam.sort_option.c_str()[0],
		  epsi,
		  invParam.ConvTestOpt,
		  invParam.PrintLevel, 
		  invParam.Neig, invParam.Nmax, stdout);
      
      // NO NEED OF THIS
      // Copy result back
      // psi = psif;

      T tt;
      (*A)(tt,psi,PLUS);
      QDPIO::cout<<"OPT_EIGBiCG_SYSSOLVER: True residual after solution : "<<sqrt(norm2(tt-chi,s))<<endl ;
      QDPIO::cout<<"OPT_EIGBiCG_SYSSOLVER: norm of  solution            : "<<sqrt(norm2(psi,s))<<endl ;
      QDPIO::cout<<"OPT_EIGBiCG_SYSSOLVER: norm of rhs                  : "<<sqrt(norm2(chi,s))<<endl ;
      


      if(!s.hasOrderedRep()){
	QDPIO::cout<<"OPPS! I have no implemented OPT_EigBiCG for Linops with non contigius subset\n";
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
  LinOpSysSolverOptEigBiCG<LatticeFermionF>::operator()(LatticeFermionF& psi, const LatticeFermionF& chi) const
  {
    return sysSolver<LatticeFermionF,float,Complex_C>(psi, chi, A, invParam);
  }

  // LatticeFermionD
  template<>
  SystemSolverResults_t
  LinOpSysSolverOptEigBiCG<LatticeFermionD>::operator()(LatticeFermionD& psi, const LatticeFermionD& chi) const
  {
    return sysSolver<LatticeFermionD, double, Complex_Z>(psi, chi, A, invParam);
  }

#if 0
 

  // Not quite ready yet for these - almost there
  // LatticeStaggeredFermionF
  template<>
  SystemSolverResults_t
  LinOpSysSolverOptEigBiCG<LatticeStaggeredFermionF>::operator()(LatticeStaggeredFermionF& psi, const LatticeStaggeredFermionF& chi) const
  {
    return sysSolver(psi, chi, A, invParam);
  }

  // LatticeStaggeredFermionD
  template<>
  SystemSolverResults_t
  LinOpSysSolverOptEigBiCG<LatticeStaggeredFermionD>::operator()(LatticeStaggeredFermionD& psi, const LatticeStaggeredFermionD& chi) const
  {
    return sysSolver(psi, chi, A, invParam);
  }
#endif

}
