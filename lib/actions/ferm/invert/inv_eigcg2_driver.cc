// $Id: inv_eigcg2_driver.cc,v 1.2 2007-09-28 02:05:08 kostas Exp $
/*! \file
 * \brief Eig CG driver
 *
 * Driver

 THIS IS USELESS... it is here to help me remember what I need to do....
 Will be removed from the repository as soon as I am done...
 */

#include "simpleGramSchmidt.h"
#include "lapack_wrapper.h"
#include "containers.h"
#include "inv_eigcg2.h"
#include "octave_debug.h"


using namespace LinAlg ;
using namespace Octave ;

namespace Chroma 
{ 

  /**
 
    TheNamedObjMap::Instance().create<EigenInfo>(params.named_obj.eigen_id);
    EigenInfo& eigenvec_val = 
      TheNamedObjMap::Instance().getData<EigenInfo>(params.named_obj.eigen_id);

    
     TheNamedObjMap::Instance().get(params.named_obj.eigen_id).setFileXML(file_xml);
    TheNamedObjMap::Instance().get(params.named_obj.eigen_id).setRecordXML(record_xml);
  **/





  void EigCGdriver(Handle< LinearOperator<LatticeFermion> >& MM,
		   Handle< LinearOperator<LatticeFermion> >& H,
		   const StathoCGParams_t& params,
		   XMLWriter& xml_out,
		   const multi1d<LatticeFermion>& rhss,
		   multi1d<LatticeFermion>& sol,
		   EigenInfo& eigenvec_val)
  {
    StopWatch snoop;
    snoop.reset();
    double Time ;

    // Try and get lowest eigenvalue of MM
    const OrderedSubset& s = MM->subset();
  
    multi1d<Double> lambda ; //the eigenvalues
    multi1d<LatticeFermion> evec(0); // The eigenvectors  
    // eigenvectors will be computed here. The array will be resized by StathoCG
    //Logic needs to be added through out the code to bring eigenvectors
    // in from outside
      
    Vectors<LatticeFermion> GoodEvec; //
    Vectors<Double> GoodEval; //
    if(params.Neig_max>0){
      GoodEval.vec.resize(params.Neig_max) ;
      GoodEvec.vec.resize(params.Neig_max) ;
      GoodEval.N = 0 ;
      GoodEvec.N = 0 ;
    }
    else{
      GoodEval.vec.resize(params.Neig) ;
      GoodEvec.vec.resize(params.Neig) ;
      GoodEval.N = 0 ;
      GoodEvec.N = 0 ;
    }
    /**
    for(int i =0; i < params.Neig; i++){
      psi[i] = zero;
      gaussian(psi[i],s);
      lambda[i] = Real(1);
    }
    **/
    
    int n_CG_count;
    int tt_n_CG;
    
    //OctavePrintClear("RayleighRich.m") ;

    //might want to give the option to do preconditioned CG
    //need different flag here ex. param.AlgType==vPrecCG
    if(params.Neig_max == 0 ){// Do the preconditioned CG
      vecPrecondCG(*MM,sol,rhss,lambda,evec,	
		   params.RsdCG, // CG residual
		   params.MaxCG, // Max CG itterations
		   tt_n_CG ); n_CG_count += tt_n_CG ;
    }
    /**
      else if(params.Neig_max == -1 ){// Do the InitCG
	InitCG(*MM,sol[i],rhss[i],lambda,evec,	
	       params.RsdCG, // CG residual
	       params.MaxCG, // Max CG itterations
	       n_CG_count );n_CG_count += tt_n_CG ;
	write(xml_out, "n_CG", tt_n_CG); 
      }
    **/
    else{ //Do the eigenvalue refinement 
      InitGuess(*MM,sol,rhss,GoodEval.vec,GoodEvec.vec,GoodEvec.N,
		tt_n_CG);
      n_CG_count += tt_n_CG ;
      if((GoodEvec.N+params.Neig)<=params.Neig_max){// if there is space for new
	evec.resize(0); // get in there with no evecs so that it computes new
	InvStathoCG(*MM, 
		    sol,
		    rhss,
		    lambda, 
		    evec, 
		    params.Neig, //Eigenvectors to keep
		    params.Nmax,  // Max vectors to work with
		    params.RsdCG, // CG residual
		    params.MaxCG, // Max CG itterations
		    tt_n_CG  // CG iteration count
		    ); n_CG_count += tt_n_CG ;
	write(xml_out, "n_CG", tt_n_CG);
	
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
      } // if there is space
      else{ // call stathoCG but ask it not to computed evecs
	evec.resize(0);
	InvStathoCG(*MM, 
		    sol,
		    rhss,
		    lambda, 
		    evec, 
		    0, //Eigenvectors to keep
		    params.Nmax,  // Max vectors to work with
		    params.RsdCG*33, // CG residual.... Empirical restart need a param here...
		    params.MaxCG, // Max CG itterations
		    tt_n_CG  // CG iteration count
		    ); n_CG_count += tt_n_CG ;
	int tt ;
	QDPIO::cout<<"Restart1: "<<endl ;
	InitGuess(*MM,sol,rhss,GoodEval.vec,GoodEvec.vec,GoodEvec.N, tt);
	InvStathoCG(*MM, 
		    sol,
		    rhss,
		    lambda, 
		    evec, 
		    0, //Eigenvectors to keep
		    params.Nmax,  // Max vectors to work with
		    params.RsdCG, // CG residual
		    params.MaxCG, // Max CG itterations
		    tt  // CG iteration count
		    ); 
	n_CG_count += tt; tt_n_CG += tt ; 
	//n_CG_count += tt; tt_n_CG += tt ; 
	write(xml_out, "n_CG", tt_n_CG);
      }
    } // else
    /**/
    // Dump output
    // xml_out << eig_spec_xml;
    write(xml_out, "n_CG_count", n_CG_count); 
    write(xml_out, "lambda_Msq", lambda); 
    
    /** This is just checking how good the eigenvectors are... Need not be done
	every time. But some how the option to do it should be there for testing
	reasons...
   // Check norms
  multi1d<Double> check_norm(GoodEvec.N);
  multi1d<Double> check_norm_rel(GoodEvec.N);
  for(int i=0; i < GoodEvec.N; i++) { 
    LatticeFermion Me;
    LatticeFermion lambda_e;
    (*MM)(Me, GoodEvec[i], PLUS);
    Real ll = GoodEval[i] ;
    lambda_e[s] = ll*GoodEvec[i];
    
    LatticeFermion r_norm;
    r_norm[s] = Me - lambda_e;
    
    check_norm[i] = norm2(r_norm,s);
    check_norm[i] = sqrt(check_norm[i]);
  }
  write(xml_out, "check_norm", check_norm);
  **/

    // Here we need to store back the vectors ....
    for(int i=0; i < GoodEvec.N; i++) {
      check_norm_rel[i] = check_norm[i]/fabs(GoodEval[i]);
    }
    write(xml_out, "check_norm_rel", check_norm_rel);
    
    for(int i=0; i < GoodEvec.N; i++){
      QDPIO::cout<<"eval["<<i<<"]="<<GoodEval[i]<<" ";
      QDPIO::cout<<"check_norm["<<i<<"]="<<check_norm[i]<<" ";
      QDPIO::cout<<"check_norm_rel["<<i<<"]="<<check_norm_rel[i]<<endl;
    }
    
    
    // this stores the eigenvectors in the eigenvector structure
    /** will reinstate once I know what the algorithm is doing **/
    eigenvec_val.getEvalues().resize(params.Neig);
    eigenvec_val.getEvectors().resize(params.Neig);
    for (int i=0; i<params.Neig; i++)
      eigenvec_val.getEvalues()[i]=GoodEval[i];
    eigenvec_val.getLargest()= -10000.00 ;// Just a random number//lambda_high_aux[0];
    for (int i=0; i<params.Neig; i++)
      eigenvec_val.getEvectors()[i]=GoodEvec[i];
    /**/
    
    //   QDPIO::cout << "Writing low eigenvalues/vectors" << endl;
    //   writeEigen(input, lambda, psi, lambda_high_aux[0], QDPIO_SERIAL);
    
  }

}
