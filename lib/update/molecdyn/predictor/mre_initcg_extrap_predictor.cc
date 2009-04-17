#include "chromabase.h"
#include "meas/inline/io/named_objmap.h"
#include "update/molecdyn/predictor/mre_initcg_extrap_predictor.h"
#include "meas/eig/gramschm.h"
#include "meas/eig/gramschm_array.h"
#include "update/molecdyn/predictor/lu_solve.h"
#include "actions/ferm/invert/containers.h"
#include "meas/eig/sn_jacob.h"

namespace Chroma 
{ 
  
  namespace MREInitCG4DChronoPredictorEnv 
  {
    namespace
    {
      // Create a new 4D Zero Guess Predictor
      // No params to read -- but preserve form
      AbsChronologicalPredictor4D<LatticeFermion>* createPredictor(XMLReader& xml,
								   const std::string& path) 
      {
	unsigned int max_chrono = 1;
	std::string opt_eigen_id;
	int nevec;
	int max_evec;
	try 
	{
	  XMLReader paramtop(xml, path);
	  read( paramtop, "./MaxChrono", max_chrono);
	  read( paramtop, "./MaxEvec", max_evec);
	  read( paramtop, "./opt_eigen_id", opt_eigen_id);
	}
	catch( const std::string& e ) { 
	  QDPIO::cerr << "Caught exception reading XML: " << e << endl;
	  QDP_abort(1);
	}
      
	return new MREInitCG4DChronoPredictor(max_chrono, opt_eigen_id, max_evec);
      }
    
      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "MRE_INITCG_4D_PREDICTOR";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= The4DChronologicalPredictorFactory::Instance().registerObject(name, createPredictor);
  	registered = true;
      }
      return success;
    }
  }


  void MREInitCG4DChronoPredictor::operator()(
					      LatticeFermion& psi,
					      const LinearOperator<LatticeFermion>& M,
					      const LatticeFermion& chi) 
  {
    START_CODE();

    // Always do this...
    find_extrap_solution(psi, M, chi);

    
    END_CODE();
  }


  void 
  MREInitCG4DChronoPredictor::find_extrap_solution(
						   LatticeFermion& psi,
						   const LinearOperator<LatticeFermion>& A,
						   const LatticeFermion& chi) 
  {
    START_CODE();

    const Subset& s= A.subset();


  

    int Nchrono = chrono_buf->size();

    QDPIO::cout << "MREInitCG Predictor: Got " << Nchrono << " chrono vecs" << endl;


    psi = zero;

    if (Nchrono > 0 ) { 

      if( Nchrono == 1 ) { 

	// If only one chrono vector exists, give that.
	LatticeFermion tmpvec;
	chrono_buf->get(0, tmpvec);
	psi[s] = tmpvec;

      }

      else {
	
	// Otherwise do minimum norm extrapolation
	multi1d<LatticeFermion> v(Nchrono);
	QDPIO::cout << "MREInitCG Predictor: Orthonormalizing the " << Nchrono << " chrono vecs" << endl;
	
	// Orthogonalize current against the last N....
	for(int i=0; i < Nchrono; i++) { 
	  
	  // Zero out the non subsetted part
	  v[i] = zero;
	  
	  // Grab the relevant vector from the chronobuf
	  // The way the circular buffer works is that i=0 is the
	  // most recent... as i increases vectors become less recent.
	  LatticeFermion tmpvec;
	  chrono_buf->get(i, tmpvec);
	  
	  if( i == 0 ) { 
	    // First vector we just take
	    v[i][s] = tmpvec;
	  }
	  else { 
	    // i-th vector. Orthogonalise against i-1 previous
	    // vector, but i is an index running from 0. So I need
	    // to pass i+1-1=i as the number of vectors to orthog against
	    //
	    // This is a very dumb GramSchmidt process and possibly
	    // unstable, however apparently (according to the paper)
	    // this is OK.
	    GramSchm(tmpvec, v, i, s);
	    GramSchm(tmpvec, v, i, s);
	    v[i][s] = tmpvec;
	  }
	  // QDPIO::cout << "Norm v[i] = " << norm2(v[i],s) << endl;
	  // Normalise v[i]
	  Double norm = sqrt(norm2(v[i], s));
	  v[i][s] /= norm;
	}

	// Now I need to form G_n m = v_[n]^{dag} A v[m]
	multi2d<DComplex> G(Nchrono,Nchrono);
	
	for(int m = 0 ; m < Nchrono; m++) { 
	  LatticeFermion tmpvec;
	  A(tmpvec, v[m], PLUS);
	
	  for(int n = 0; n < Nchrono; n++) { 
	    G(n,m) = innerProduct(v[n], tmpvec, s);
	  }
	}

	// Now I need to form b_n = v[n]^{dag} chi
	multi1d<DComplex> b(Nchrono);
	
	for(int n = 0; n < Nchrono; n++) { 
	  b[n] = innerProduct(v[n], chi, s);
	}

	// Solve G_nm a_m = b_n:
	
	// First LU decompose G in place 
	multi1d<DComplex> a(Nchrono);
	
	LUSolve(a, G, b);
	
	// Create teh lnear combination
	psi[s] = Complex(a[0])*v[0];
	for(int n=1; n < Nchrono; n++) { 
	  psi[s] += Complex(a[n])*v[n];
	}
      }
    }

    QDPIO::cout << "MRE InitCG predictor: psi prepared " << endl;

    // Strategy: Rayleigh Ritz all the EVs, put back Neig of them.
    LinAlg::OptEigInfo& eiginfo = TheNamedObjMap::Instance().getData< LinAlg::OptEigInfo >(opt_eigen_id);

    int Nevec = eiginfo.ncurEvals;
    QDPIO::cout << "There are " << Nevec << " current EVs" << endl;

    if( Nevec > 0 ) {
      LatticeFermion tmpvec;
      multi1d<LatticeFermion> ev(Nevec); 
      for(int i=0; i < Nevec; i++) { 

	// The EigCG rather conveniently Orthonormalizes 
	// evs for me, so I don't need to. Just copy them into ev
	ev[i] = zero;
	eiginfo.CvToLatFerm(ev[i],s,i);
	
      }
    
      QDPIO::cout << "MREInitCG: Making up little matrix" << endl;

      multi1d<Real> lambda(Nevec);
      for(int i=0; i < Nevec;i++) { 
	A(tmpvec, ev[i], PLUS);
	lambda[i] = innerProductReal(ev[i], tmpvec, s);
      }


      multi1d<Complex> offd(Nevec*(Nevec-1)/2);
      int ij=0;
      for(int i=0; i < Nevec; i++) { 
	A(tmpvec, ev[i], PLUS);
	for(int j=0; j < i;  j++) { 
	  offd[ij] = innerProduct(ev[j], tmpvec);
	  ij++;
	}
      }
      
      QDPIO::cout << "MREInitCG: Calling Jacobi Diagonalizer " << endl;
      Real rsda=Real(1.0e-5);
      int n_jacob = SN_Jacob(ev, Nevec, lambda, offd,
			     rsda, 50, s);
      
	
      int vecs_to_copy=Nevec;  // Copy back all Nevec
      if( Nevec > Neig ) { 
	vecs_to_copy = Neig;  // If we have more vecs than max, only copy Max back
      }

      QDPIO::cout << "MREInitCG: Copying " << vecs_to_copy << " vecs into EigInfo " << endl;
      // EigCG Will Recompute necessary vectors & Matrices...
      // Copy the Nvec chrono vectors into the buffer 
      for(int i=0; i < vecs_to_copy; i++) { 
	eiginfo.CvToEigCGvec(ev[i],s,i);
      }

      QDPIO::cout << "MREInitCG: Resetting nCurVal" << endl;
      eiginfo.ncurEvals = vecs_to_copy;
      
    }
  

    END_CODE();
  }




    
}
