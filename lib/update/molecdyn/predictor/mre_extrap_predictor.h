// -*- C++ -*-
/*! \file
 * \brief Minimal residual predictor
 *
 * Predictors for HMC
 */

#ifndef __mre_extrap_predictor_h__
#define __mre_extrap_predictor_h__

#include "chromabase.h"
#include "handle.h"
#include "update/molecdyn/predictor/chrono_predictor.h"
#include "update/molecdyn/predictor/chrono_predictor_factory.h"
#include "update/molecdyn/predictor/circular_buffer.h"
#include "update/molecdyn/predictor/lu_solve.h"
#include "meas/eig/gramschm.h"

namespace Chroma 
{ 
  
  /*! @ingroup predictor */
  namespace MinimalResidualExtrapolation4DChronoPredictorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Minimal residual predictor
  /*! @ingroup predictor */
  template<typename T>
  class MinimalResidualExtrapolation4DChronoPredictor  
    : public AbsTwoStepChronologicalPredictor4D<T> 
  {
  private:
    Handle< CircularBuffer<T> > chrono_bufX;
    Handle< CircularBuffer<T> > chrono_bufY;

  void 
  find_extrap_solution(
		       T& psi,
		       const LinearOperator<T>& M,
		       const T& chi,
		       const Handle<CircularBuffer<T> >& chrono_buf, 
		       enum PlusMinus isign) 
    {
      START_CODE();
      
      const Subset& s= M.subset();
      
      
#if 0
      {
	T r;
	r[s] = chi;
	T tmp;
	M(tmp, psi, isign);
	r[s] -= tmp;
	Double norm_r = sqrt(norm2(r,s));
	Double norm_chi = sqrt(norm2(chi,s));
	QDPIO::cout << "MRE Predictor: before prediction || r || / || b || =" << norm_r/norm_chi << endl;
      }
#endif
      
      int Nvec = chrono_buf->size();
      
      
      // Construct an orthonormal basis from the 
      // vectors in the buffer. Stick to notation of paper and call these
      // v
      multi1d<T> v(Nvec);
      
      for(int i=0; i < Nvec; i++) { 
	// Zero out the non subsetted part
	v[i] = zero;
	
	// Grab the relevant vector from the chronobuf
	T tmpvec;
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
	T tmpvec;
	M(tmpvec, v[m], isign);
	
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
      
      
#if 0
      {
	T r;
	r[s] = chi;
	T tmp;
	M(tmp, psi, isign);
	r[s] -= tmp;
	Double norm_r = sqrt(norm2(r,s));
	Double norm_chi = sqrt(norm2(chi,s));
	QDPIO::cout << "MRE Predictor: after prediction || r || / || b || =" << norm_r/norm_chi << endl;
      }
#endif
      
      END_CODE();
    }



  public:
    
    MinimalResidualExtrapolation4DChronoPredictor(unsigned int max_chrono) : 
      chrono_bufX(new CircularBuffer<T>(max_chrono)),
      chrono_bufY(new CircularBuffer<T>(max_chrono)) {}
    
    // Destructor is automagic
    ~MinimalResidualExtrapolation4DChronoPredictor(void) {}

    void predictX(T& X,
		  const LinearOperator<T>& M,
		  const T& chi) 
    {
      START_CODE();
      StopWatch swatch;
      swatch.reset();
      swatch.start();
      
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
	  
	  // Expect M is either  MdagM if we use chi
	  // or                   M    if we minimize against Y
	  find_extrap_solution(X, M, chi, chrono_bufX, PLUS);
	}
	break;
      }
      
      swatch.stop();
      QDPIO::cout << "MRE_PREDICT_X_TIME = " << swatch.getTimeInSeconds() << " s" << endl;
      
      END_CODE();
    }
    
    void predictY(T& Y,
		  const LinearOperator<T>& M,
		  const T& chi) 
    {
      START_CODE();
      StopWatch swatch;
      swatch.reset();
      swatch.start();
      
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
      
      swatch.stop();
      QDPIO::cout << "MRE_PREDICT_Y_TIME = " << swatch.getTimeInSeconds() << " s" << endl;
      END_CODE();
    }
    
    // No internal state so reset is a nop
    void reset(void) {
      chrono_bufX->reset();
      chrono_bufY->reset();
    }


    void newXVector(const T& X) 
    {
      START_CODE();

      QDPIO::cout << "MREPredictor: registering new X solution. " << endl;
      chrono_bufX->push(X);
      QDPIO::cout << "MREPredictor: number of X vectors stored is = " << chrono_bufX->size() << endl;
    
      END_CODE();
    }


    void newYVector(const T& Y) 
    {
      START_CODE();

      QDPIO::cout << "MREPredictor: registering new Y solution. " << endl;
      chrono_bufY->push(Y);
      QDPIO::cout << "MREPredictor: number of Y vectors stored is = " << chrono_bufY->size() << endl;
    
      END_CODE();
    }

    void replaceXHead(const T& v)
    {
      chrono_bufX->replaceHead(v);
    }

    void replaceYHead(const T& v)
    {
      chrono_bufY->replaceHead(v);
    }


  };

  

  /*! @ingroup predictor */
  namespace MinimalResidualExtrapolation5DChronoPredictorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  
  //! Minimal residual predictor
  /*! @ingroup predictor */
  class MinimalResidualExtrapolation5DChronoPredictor :
    public AbsChronologicalPredictor5D<LatticeFermion> {
    
  private: 
    Handle< CircularBufferArray<LatticeFermion>  > chrono_buf;
    const int N5;

    void find_extrap_solution(multi1d<LatticeFermion>& psi, 
			      const LinearOperatorArray<LatticeFermion>& A,
			      const multi1d<LatticeFermion>& chi);

  public:
    MinimalResidualExtrapolation5DChronoPredictor(const int N5_, const unsigned int max_chrono) : chrono_buf( new CircularBufferArray<LatticeFermion>(max_chrono, N5_) ), N5(N5_) {}

      
    ~MinimalResidualExtrapolation5DChronoPredictor(void) {}
    
    // Copying
    MinimalResidualExtrapolation5DChronoPredictor(const MinimalResidualExtrapolation5DChronoPredictor& p) : chrono_buf(p.chrono_buf), N5(p.N5) {}

    // Zero out psi -- it is a zero guess after all
    void operator()(multi1d<LatticeFermion>& psi,
		    const LinearOperatorArray<LatticeFermion>& A,
		    const multi1d<LatticeFermion>& chi); 

    // No internal state so reset is a Nop
    void reset(void) {
      chrono_buf->reset();

    }

    // Ignore new vector
    // Ignore new vector
    void newVector(const multi1d<LatticeFermion>& psi) 
    {
      START_CODE();

      QDPIO::cout << "MRE Predictor: registering new solution. " << endl;
      chrono_buf->push(psi);
      QDPIO::cout << "MRE Predictor: number of vectors stored is = " << chrono_buf->size() << endl;
    
      END_CODE();
    }

  };
  
} // End Namespace Chroma

#endif 
