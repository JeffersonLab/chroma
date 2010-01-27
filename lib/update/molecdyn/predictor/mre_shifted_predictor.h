// -*- C++ -*-
// $Id: mre_extrap_predictor.h,v 3.5 2009-06-01 20:39:37 bjoo Exp $
/*! \file
 * \brief Minimal residual predictor
 *
 * Predictors for HMC
 */

#ifndef __mre_shifted_predictor_h__
#define __mre_shifted_predictor_h__

#include "chromabase.h"
#include "handle.h"
#include "update/molecdyn/predictor/chrono_predictor.h"
#include "update/molecdyn/predictor/chrono_predictor_factory.h"
#include "update/molecdyn/predictor/circular_buffer.h"
#include "update/molecdyn/predictor/lu_solve.h"
namespace Chroma 
{ 
  
  /*! @ingroup predictor */
  namespace MinimalResidualExtrapolationShifted4DChronoPredictorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }

  //! Minimal residual predictor
  /*! @ingroup predictor */
  template<typename T, typename R>
  class MinimalResidualExtrapolationShifted4DChronoPredictor  
  {
  private:
    Handle< CircularBuffer<T> > chrono_buf;
    Handle< CircularBuffer<T> > chrono_bufM;
    const LinearOperator<T>& M;
    

  void 
  find_extrap_solution(
		       T& psi,
		       const T& chi,
		       const R& shift,
		       enum PlusMinus isign) 
    {
      START_CODE();
      
      const Subset& s= M.subset();
      
      
#if 1
      T rvec;
      {
	
	rvec[s] = chi;
	T tmp;
	M(tmp, psi, isign);
	tmp[s] += shift*psi;
	rvec[s] -= tmp;
	// Double norm_r = sqrt(norm2(rvec,s));
	// Double norm_chi = sqrt(norm2(chi,s));
	// QDPIO::cout << "MRE Predictor: before prediction || r || / || b || =" << norm_r/norm_chi << endl;
      }
#endif
      
      int Nvec = chrono_buf->size();
      
      
      // Construct an orthonormal basis from the 
      // vectors in the buffer. Stick to notation of paper and call these
      // v
      //      multi1d<T> v(Nvec);
      
      
      // Now I need to form G_n m = v_[n]^{dag} A v[m]
      multi2d<DComplex> G(Nvec,Nvec);
      multi1d<DComplex> b(Nvec);
      
      for(int m = 0 ; m < Nvec; m++) { 
	T v_m;
	T Mv_m;
	chrono_buf->get(m, v_m);
	chrono_bufM->get(m, Mv_m);
	//	M(Mv_m, v_m, isign);
	b[m] = innerProduct(v_m, rvec, s);
	for(int n = 0; n < Nvec; n++) { 
	  T v_n;
	  chrono_buf->get(n, v_n);
	  G(n,m) = innerProduct(v_n, Mv_m, s);
	  DComplex dcshift(shift);
	  G(n,m) += dcshift*innerProduct(v_n, v_m,s);
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
      
      // Create teh lnear combination
      {
	T v;
	chrono_buf->get(0,v);
	psi[s] += Complex(a[0])*v;
	for(int n=1; n < Nvec; n++) { 
	  chrono_buf->get(n, v);
	  psi[s] += Complex(a[n])*v;
	}
      }
      
#if 0
      {
	//	T rvec;
	rvec[s] = chi;
	T tmp;
	M(tmp, psi, isign);
	tmp[s] += shift*psi;

	rvec[s] -= tmp;
	Double norm_r = sqrt(norm2(rvec,s));
	Double norm_chi = sqrt(norm2(chi,s));
	QDPIO::cout << "MRE Predictor: after prediction || r || / || b || =" << norm_r/norm_chi << endl;
      }
#endif
      
      END_CODE();
    }



  public:
    
    MinimalResidualExtrapolationShifted4DChronoPredictor(unsigned int max_chrono, const LinearOperator<T>& M_) : 
      chrono_buf(new CircularBuffer<T>(max_chrono)),
      chrono_bufM(new CircularBuffer<T>(max_chrono)),
      M(M_) {}
    
    // Destructor is automagic
    ~MinimalResidualExtrapolationShifted4DChronoPredictor(void) {}

    void predictX(T& X,
		  const R& shift,
		  const T& chi) 
    {
      START_CODE();
      StopWatch swatch;
      swatch.reset();
      swatch.start();
      
      int Nvec = chrono_buf->size();
      switch(Nvec) { 
      case 0:
	{
	  return;
	}
	break;
#if 0
      case 1:
	{
	  QDPIO::cout << "MRE Predictor: Only 1 vector stored. Giving you last solution " << endl;
	  chrono_buf->get(0,X);
	}
	break;
#endif
      default:
	{
	  QDPIO::cout << "MRE Predictor: Finding X extrapolation with "<< Nvec << " vectors" << endl;
	  
	  // Expect M is either  MdagM if we use chi
	  // or                   M    if we minimize against Y
	  find_extrap_solution(X, chi, shift, PLUS);
	}
	break;
      }
      
      swatch.stop();
      QDPIO::cout << "MRE_PREDICT_X_TIME = " << swatch.getTimeInSeconds() << " s" << endl;
      
      END_CODE();
    }
    

    
    // No internal state so reset is a nop
    void reset(void) {
      chrono_buf->reset();
    }


    void newXVector(const T& X) 
    {
      START_CODE();
      const Subset& s = M.subset();
      chrono_buf->push(X);
      T Mv;
      M(Mv,X,PLUS);
      chrono_bufM->push(Mv);


      QDPIO::cout << "MREPredictor: number of X vectors stored is = " << chrono_buf->size() << endl;

      END_CODE();
    }



    void replaceXHead(const T& v_)
    {
      START_CODE();
      const Subset& s = M.subset();

      chrono_buf->replaceHead(v_);

      T Mv;
      M(Mv, v_, PLUS);
      chrono_bufM->replaceHead(Mv);


      QDPIO::cout << "MREPredictor: number of X vectors stored is = " << chrono_buf->size() << endl;

      END_CODE();
    }


  };

    
} // End Namespace Chroma

#endif 
