#include "chromabase.h"
#include "update/molecdyn/predictor/mre_extrap_predictor.h"
#include "meas/eig/gramschm.h"
#include "meas/eig/gramschm_array.h"
#include "update/molecdyn/predictor/lu_solve.h"


namespace Chroma 
{ 
  
  namespace MinimalResidualExtrapolation4DChronoPredictorEnv 
  {
    namespace
    {
      // Create a new 4D Zero Guess Predictor
      // No params to read -- but preserve form
      AbsChronologicalPredictor4D<LatticeFermion>* createPredictor(XMLReader& xml,
								   const std::string& path) 
      {
	unsigned int max_chrono = 1;
	
	try 
	{
	  XMLReader paramtop(xml, path);
	  read( paramtop, "./MaxChrono", max_chrono);
	}
	catch( const std::string& e ) { 
	  QDPIO::cerr << "Caught exception reading XML: " << e << endl;
	  QDP_abort(1);
	}
      
	return new MinimalResidualExtrapolation4DChronoPredictor(max_chrono);
      }
    
      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "MINIMAL_RESIDUAL_EXTRAPOLATION_4D_PREDICTOR";

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


  void MinimalResidualExtrapolation4DChronoPredictor::predictX(LatticeFermion& X,
							       const LinearOperator<LatticeFermion>& M,
							       const LatticeFermion& chi) 
  {
    START_CODE();

    int Nvec = chrono_bufX->size();
    switch(Nvec) { 
    case 0:
      {
	QDPIO::cout << "MRE Predictor: Zero vectors stored. Giving you zero guess" << endl;
	X= zero;
      }
      break;
    case 1:
      {
	QDPIO::cout << "MRE Predictor: Only 1 vector stored. Giving you last solution " << endl;
	chrono_bufX->get(0,X);
      }
      break;
    default:
      {
	QDPIO::cout << "MRE Predictor: Finding X extrapolation with "<< Nvec << " vectors" << endl;
	
	// Expect M is actually M^\dagger M here 
	find_extrap_solution(X, M, chi, chrono_bufX, PLUS);
      }
      break;
    }
    
    END_CODE();
  }

  void MinimalResidualExtrapolation4DChronoPredictor::predictY(LatticeFermion& Y,
							       const LinearOperator<LatticeFermion>& M,
							       const LatticeFermion& chi) 
  {
    START_CODE();

    int Nvec = chrono_bufY->size();
    switch(Nvec) { 
    case 0:
      {
	QDPIO::cout << "MRE Predictor: Zero vectors stored. Giving you zero guess" << endl;
	Y = zero;
      }
      break;
    case 1:
      {
	QDPIO::cout << "MRE Predictor: Only 1 vector stored. Giving you last solution " << endl;
	chrono_bufY->get(0,Y);
      }
      break;
    default:
      {
	QDPIO::cout << "MRE Predictor: Finding Y extrapolation with "<< Nvec << " vectors" << endl;
	// Should have M as just M (not M^\dagger M) here.
	find_extrap_solution(Y, M, chi, chrono_bufY, MINUS);
      }
      break;
    }
    
    END_CODE();
  }


  void 
  MinimalResidualExtrapolation4DChronoPredictor::find_extrap_solution(
					  LatticeFermion& psi,
					  const LinearOperator<LatticeFermion>& A,
					  const LatticeFermion& chi,
					  const Handle<CircularBuffer<LatticeFermion> >& chrono_buf, 
					  enum PlusMinus isign) 
  {
    START_CODE();

    const Subset& s= A.subset();

    
    int Nvec = chrono_buf->size();


    // Construct an orthonormal basis from the 
    // vectors in the buffer. Stick to notation of paper and call these
    // v
    multi1d<LatticeFermion> v(Nvec);
    
    for(int i=0; i < Nvec; i++) { 
      // Zero out the non subsetted part
      v[i] = zero;

      // Grab the relevant vector from the chronobuf
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
	v[i][s] = tmpvec;
      }
      // QDPIO::cout << "Norm v[i] = " << norm2(v[i],s) << endl;
      // Normalise v[i]
      Double norm = sqrt(norm2(v[i], s));
      v[i][s] /= norm;
			 
    }

    // Now I need to form G_n m = v_[n]^{dag} A v[m]
    multi2d<DComplex> G(Nvec,Nvec);

    for(int m = 0 ; m < Nvec; m++) { 
      LatticeFermion tmpvec;
      A(tmpvec, v[m], isign);

      for(int n = 0; n < Nvec; n++) { 
	G(n,m) = innerProduct(v[n], tmpvec, s);
      }
    }

    // Now I need to form b_n = v[n]^{dag} chi
    multi1d<DComplex> b(Nvec);

    for(int n = 0; n < Nvec; n++) { 
      b[n] = innerProduct(v[n], chi, s);
    }

    // Solve G_nm a_m = b_n:

    // First LU decompose G in place 
    multi1d<DComplex> a(Nvec);

    LUSolve(a, G, b);

    // Check solution 
    multi1d<DComplex> Ga(Nvec);

    for(int i=0; i < Nvec; i++) { 
      Ga[i] = G(i,0)*a[0];
      for(int j=1; j < Nvec; j++) { 
	Ga[i] += G(i,j)*a[j];
      }
    }

    multi1d<DComplex> r(Nvec);
    for(int i=0; i < Nvec; i++) { 
      r[i] = b[i]-Ga[i];
    }

#if 0
    QDPIO::cout << "Constraint Eq Solution Check" << endl;
    for(int i=0; i < Nvec; i++) { 
      QDPIO::cout << "   r[ " << i << "] = " << r[i] << endl;
    }
#endif

    // Create teh lnear combination
    psi[s] = Complex(a[0])*v[0];
    for(int n=1; n < Nvec; n++) { 
      psi[s] += Complex(a[n])*v[n];
    }

    END_CODE();
  }





  namespace MinimalResidualExtrapolation5DChronoPredictorEnv 
  {
    namespace
    {
      // Create a new 5D Zero Guess Predictor
      // No params to read 
      AbsChronologicalPredictor5D<LatticeFermion>* createPredictor(const int N5,
								   XMLReader& xml,
								   const std::string& path) 
      {
	unsigned int max_chrono = 1;

	try 
	{
	  XMLReader paramtop(xml, path);
	  read( paramtop, "./MaxChrono", max_chrono);
	}
	catch( const std::string& e ) { 
	  QDPIO::cerr << "Caught exception reading XML: " << e << endl;
	  QDP_abort(1);
	}
	
	return new MinimalResidualExtrapolation5DChronoPredictor(N5,max_chrono);
      }
      
      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "MINIMAL_RESIDUAL_EXTRAPOLATION_5D_PREDICTOR";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= The5DChronologicalPredictorFactory::Instance().registerObject(name, createPredictor);
      	registered = true;
      }
      return success;
    }
  }



  void MinimalResidualExtrapolation5DChronoPredictor::operator()(
		    multi1d<LatticeFermion>& psi,
		    const LinearOperatorArray<LatticeFermion>& M,
		    const multi1d<LatticeFermion>& chi) 
  {
    START_CODE();

    int Nvec = chrono_buf->size();
    switch(Nvec) { 
    case 0:
      {
	QDPIO::cout << "MRE Predictor: Zero vectors stored. Giving you zero guess" << endl;
	psi = zero;
      }
      break;
    case 1:
      {
	QDPIO::cout << "MRE Predictor: Only 1 vector stored. Giving you last solution " << endl;
	chrono_buf->get(0,psi);
      }
      break;
    default:
      {
	QDPIO::cout << "MRE Predictor: Finding  extrapolation with "<< Nvec << " vectors" << endl;
	find_extrap_solution(psi, M, chi);
      }
      break;
    }
    
    END_CODE();
  }


  void 
  MinimalResidualExtrapolation5DChronoPredictor::find_extrap_solution(
		       	  multi1d<LatticeFermion>& psi,
			  const LinearOperatorArray<LatticeFermion>& A,
			  const multi1d<LatticeFermion>& chi) 
  {
    START_CODE();

    const Subset& s= A.subset();
    
    int Nvec = chrono_buf->size();


    // Construct an orthonormal basis from the 
    // vectors in the buffer. Stick to notation of paper and call these
    // v
    multi2d<LatticeFermion> v(Nvec, N5);
    
    for(int i=0; i < Nvec; i++) { 
      // Zero out the non subsetted part
      for(int d5=0; d5 < N5; d5++) { 
	v[i][d5] = zero;
      }

      // Grab the relevant vector from the chronobuf
      multi1d<LatticeFermion> tmpvec(N5);
      chrono_buf->get(i, tmpvec);

      if( i == 0 ) { 
	// First vector we just take
	// Loop over 5th dim
	for(int d5=0; d5 < N5; d5++) { 
	  v[i][d5][s] = tmpvec[d5];
	}
      }
      else { 
	// i-th vector. Orthogonalise against i-1 previous
	// vector, but i is an index running from 0. So I need
	// to pass i+1-1=i as the number of vectors to orthog against
	//
	// This is a very dumb GramSchmidt process and possibly
	// unstable, however apparently (according to the paper)
	// this is OK.
	GramSchmArray(tmpvec, v, i, s);

	// Loop over 5th dim
	for(int d5=0; d5 < N5; d5++) {
	  v[i][d5][s] = tmpvec[d5];
	}
      }

      // Normalise v[i]
      Double norm = norm2(v[i][0], s);

      // Loop over 5th dim and accumulate
      for(int d5=1; d5 < N5; d5++) { 
	norm += norm2(v[i][d5], s);
      }

      // Square root
      Double root_norm = sqrt(norm);

      // Normalise -- Loop over 5th dim
      for(int d5=0; d5 < N5; d5++) { 
	v[i][d5][s] /= root_norm;
      }

			 
    }

    // Now I need to form G_n m = v_[n]^{dag} A v[m]
    multi2d<DComplex> G(Nvec,Nvec);

    for(int m = 0 ; m < Nvec; m++) { 
      multi1d<LatticeFermion> tmpvec(N5);

      // 5D Matrix application
      A(tmpvec, v[m], PLUS);


      for(int n = 0; n < Nvec; n++) { 

	// Compute 5D ioner product
	G(n,m) = innerProduct(v[n][0], tmpvec[0], s);
	for(int d5=1; d5 < N5; d5++) {
	  G(n,m) += innerProduct(v[n][d5], tmpvec[d5], s);
	}

      }
    }

    // Now I need to form b_n = v[n]^{dag} chi
    multi1d<DComplex> b(Nvec);

    for(int n = 0; n < Nvec; n++) { 

      // Loop over 5th Dim
      b[n] = innerProduct(v[n][0], chi[0], s);
      for(int d5=1; d5 < N5; d5++) { 
	b[n] += innerProduct(v[n][d5], chi[d5], s);
      }
    }

    // Solve G_nm a_m = b_n:

    // First LU decompose G in place 
    multi1d<DComplex> a(Nvec);

    LUSolve(a, G, b);

#if 0
    // Check solution 
    multi1d<DComplex> Ga(Nvec);


    for(int i=0; i < Nvec; i++) { 
      Ga[i] = G(i,0)*a[0];
      for(int j=1; j < Nvec; j++) { 
	Ga[i] += G(i,j)*a[j];
      }
    }

    multi1d<DComplex> r(Nvec);
    for(int i=0; i < Nvec; i++) { 
      r[i] = b[i]-Ga[i];
    }

    QDPIO::cout << "Constraint Eq Solution Check" << endl;
    for(int i=0; i < Nvec; i++) { 
      QDPIO::cout << "   r[ " << i << "] = " << r[i] << endl;
    }
#endif

    // Create the lnear combination

    // First piece
    for(int d5=0; d5 < N5; d5++) { 
      psi[d5][s] = Complex(a[0])*v[0][d5];
    

      for(int n=1; n < Nvec; n++) { 
	psi[d5][s] += Complex(a[n])*v[n][d5];
      }

    }
    
    END_CODE();
  }
    
}
