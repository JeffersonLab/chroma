#include "chromabase.h"
#include "update/molecdyn/predictor/mre_extrap_predictor.h"
#include "meas/eig/gramschm.h"
#include "meas/eig/gramschm_array.h"
#include "update/molecdyn/predictor/lu_solve.h"


namespace Chroma { 
  
  namespace MinimalResidualExtrapolation4DChronoPredictorEnv {

    // Create a new 4D Zero Guess Predictor
    // No params to read -- but preserve form
    AbsChronologicalPredictor4D<LatticeFermion>* createPredictor(XMLReader& xml,
								 const std::string& path) {

      unsigned int max_chrono=1;

      try { 
	XMLReader paramtop(xml, path);
	read( paramtop, "./MaxChrono", max_chrono);
      }
      catch( const std::string& e ) { 
	QDPIO::cerr << "Caught exception reading XML: " << e << endl;
	QDP_abort(1);
      }
      
      return new MinimalResidualExtrapolation4DChronoPredictor(max_chrono);
    }
    
    const std::string name = "MINIMAL_RESIDUAL_EXTRAPOLATION_4D_PREDICTOR";

    // Register it
    const bool registered = The4DChronologicalPredictorFactory::Instance().registerObject(name, createPredictor);
  
  };

  void MinimalResidualExtrapolation4DChronoPredictor::operator()(
								 LatticeFermion& psi,
								 const LinearOperator<LatticeFermion>& M,
								 const LatticeFermion& chi) 
  {

    switch(chrono_buf->size()) { 
    case 0:
      {
	QDPIO::cout << "Minimal Residual Extrapolation: zero vectors stored. Giving you zero guess" << endl;
	psi = zero;
      }
      break;
    case 1:
      {
	QDPIO::cout << "Minimal Residual Extrapolation: Only 1 vector stored. Giving you last solution " << endl;
	chrono_buf->get(0,psi);
      }
      break;
    default:
      {
	QDPIO::cout << "Minimal Residual Extrapolation: Trying to find you the minimal residual extrapolation " << endl;
	find_extrap_solution(psi, M, chi);
      }
      break;
    }
  }


  void 
  MinimalResidualExtrapolation4DChronoPredictor::find_extrap_solution(
					  LatticeFermion& psi,
					  const LinearOperator<LatticeFermion>& A,
					  const LatticeFermion& chi) 
  {
    const OrderedSubset& s= A.subset();

    
    int Nvec = chrono_buf->size();

    QDPIO::cout << "Number of chrono vectors: " << Nvec << endl;

    // Construct an orthogonal (but not orthonormal basis from the 
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
    }

    // Now I need to form G_n m = v_[n]^{dag} A v[m]
    multi2d<DComplex> G(Nvec,Nvec);

    for(int m = 0 ; m < Nvec; m++) { 
      LatticeFermion tmpvec;
      A(tmpvec, v[m], PLUS);

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

    QDPIO::cout << "Constraint Eq Solution Check" << endl;
    for(int i=0; i < Nvec; i++) { 
      QDPIO::cout << "   r[ " << i << "] = " << r[i] << endl;
    }

    // Create teh lnear combination
    psi[s] = Complex(a[0])*v[0];
    for(int n=1; n < Nvec; n++) { 
      psi[s] += Complex(a[n])*v[n];
    }

  }




    /*
      namespace MinimalResidualExtrapolation5DChronoPredictorEnv {
      // Create a new 5D Zero Guess Predictor
      // No params to read 
      AbsChronologicalPredictor5D<LatticeFermion>* createPredictor(const int N5,
      XMLReader& xml,
      const std::string& path) {
      return new MinimalResidualExtrapolation5DChronoPredictor(N5);
      };
      
      const std::string name = "LINEAR_EXTRAPOLATION_5D_PREDICTOR";
      
      const bool registered = The5DChronologicalPredictorFactory::Instance().registerObject(name, createPredictor);
      
      };
    */
};
