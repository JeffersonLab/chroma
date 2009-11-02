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
      CopyToLatFerm<R>(XX,px,arg.LinOp->subset());
      (*arg.LinOp)(arg.YY,arg.XX,isign) ;

      //copy back..
      CopyFromLatFerm<R>(py,YY,arg.LinOp->subset());
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
				    const SysSolverOptEigBiCGParams& invParam)
    {
      START_CODE();

      SystemSolverResults_t res;  // initialized by a constructor

      LinAlg::OptEigInfo& EigInfo = TheNamedObjMap::Instance().getData< LinAlg::OptEigInfo >(invParam.eigen_id);

      QDPIO::cout<<"EigBiInfo.N= "<<EigBiInfo.N<<endl ;
      QDPIO::cout<<"EigBiInfo.lde= "<<EigBiInfo.lde<<endl ;
      QDPIO::cout<<"EigBiInfo.ldh= "<<EigBiInfo.evals.size()<<endl ;
      QDPIO::cout<<"EigBiInfo.ncurEvals= "<<EigBiInfo.ncurEvals<<endl ;
      QDPIO::cout<<"EigBiInfo.restartTol= "<<EigBiInfo.restartTol<<endl ;

      Subset s = A.subset() ;

      Complex_C *X ; 
      Complex_C *B ; 
      Complex_C *work=NULL  ;
      Complex_C *VL=NULL    ;
      Complex_C *VR=NULL    ;
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
	QDPIO::cout<<"OPPS! I have not implemented OPT_EigBiCG for Linops with non contigius subset\n";
	exit(1);
      }

      Complex_C *evecsL = (Complex_C *) &EigBiInfo.evecsL[0] ;
      Complex_C *evecsR = (Complex_C *) &EigBiInfo.evecsR[0] ;
      Complex_C *evals      = (Complex_C *) &EigInfo.evals[0].elem() ;
      Complex_C *H      = (Complex_C *) &EigBiInfo.H[0] ;
      MatVecArg<T> arg ;
      arg.LinOp = A ;
      int esize = invParam.esize*Layout::sitesOnNode()*Nc*Ns ;

      QDPIO::cout<<"OPT_EIGBICG_SYSSOLVER= "<<esize<<endl ;
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

      IncrEigbicgC(EigInfo.N, EigInfo.lde, 1, X, B, 
		   &EigInfo.ncurEvals, EigInfo.evals.size(), 
		   evecsL, evecsR, evals, H, 
		   MatrixMatvec<float,T,PLUS>, // mat vec
		   MatrixMatvec<float,T,MINUS> , // hermitian conjugate mat vec
		   (void *)&arg, 
		   AnormEst,
		   work, 
		   VL, EigInfo.lde, VR, EigInfo.lde,
		   ework, esize, 
		   resid, &restartTol,
		   invParam.MaxCG, 
		   invParam.sort_option.c_str(),
		   invParam.epsi,
		   invParam.ConfTestOpt,
		   invParam.PrintLevel, 
		   invParam.Neig, invParam.Nmax, stdout);
      
      // Copy result back
      psi = psif;

      T tt;
      (*LinOp)(tt,psi,PLUS);
      QDPIO::cout<<"OPT_EICG_SYSSOLVER: True residual after solution : "<<sqrt(norm2(tt-chi,s))<<endl ;
      QDPIO::cout<<"OPT_EICG_SYSSOLVER: norm of  solution            : "<<sqrt(norm2(psi,s))<<endl ;
      QDPIO::cout<<"OPT_EICG_SYSSOLVER: norm of rhs                  : "<<sqrt(norm2(chi,s))<<endl ;
      


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
    return sysSolver(psi, chi, *A, LinOp, invParam);
  }

  // LatticeFermionD
  template<>
  SystemSolverResults_t
  LinOpSysSolverOptEigBiCG<LatticeFermionD>::operator()(LatticeFermionD& psi, const LatticeFermionD& chi) const
  {
    return sysSolver(psi, chi, *A, LinOp, invParam);
  }

#if 1
 

  // Not quite ready yet for these - almost there
  // LatticeStaggeredFermionF
  template<>
  SystemSolverResults_t
  LinOpSysSolverOptEigBiCG<LatticeStaggeredFermionF>::operator()(LatticeStaggeredFermionF& psi, const LatticeStaggeredFermionF& chi) const
  {
    return sysSolver(psi, chi, *A, LinOp, invParam);
  }

  // LatticeStaggeredFermionD
  template<>
  SystemSolverResults_t
  LinOpSysSolverOptEigBiCG<LatticeStaggeredFermionD>::operator()(LatticeStaggeredFermionD& psi, const LatticeStaggeredFermionD& chi) const
  {
    return sysSolver(psi, chi, *A, LinOp, invParam);
  }
#endif

}
