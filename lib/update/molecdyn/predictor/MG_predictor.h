// MG_predictor.h
// Chrono predictor which will be used for incorporating MG into HMC.
/*
 * 
 *
 * Predictors for HMC
 */

#ifndef __MG_predictor_h__
#define __MG_predictor_h__

#include "chromabase.h"
#include "update/molecdyn/predictor/abs_MG_chrono_predictor.h"
#include "update/molecdyn/predictor/chrono_predictor_factory.h"


namespace Chroma 
{ 
  
  namespace MG4DChronoPredictorEnv 
  {
    extern const std::string name;
    bool registerAll();
  }


  // MG predictor
  template<typename T>
  class MG4DChronoPredictor : 
    public AbsMGChronologicalPredictor4D<LatticeFermion> 
  {
   
  public:

    // Destructor is automagic
    ~MG4DChronoPredictor(void) {}

    MG4DChronoPredictor(void) : subspace_name("Default"),
					  refresh_rate(0){}

    MG4DChronoPredictor(const MG4DChronoPredictor& p) :
      subspace_name(p.subspace_name),
      refresh_rate(p.refresh_rate){}

    MG4DChronoPredictor(const std::string subspace_name_, const int refresh_rate_) :
      subspace_name(subspace_name_),
      refresh_rate(refresh_rate_){}

    void getSubspace() 
    {
      START_CODE();

      QDPIO::cout << "MG4DChronoPredictor - Subspace Name: "<<subspace_name<<std::endl;
    
      END_CODE();
    }

    void resetSubspace(int counter) 
    { 
      START_CODE();

      QDPIO::cout << "MG4DChronoPredictor - Current Counter: "<<counter<<std::endl;
      QDPIO::cout << "MG4DChronoPredictor - Refresh Rate: "<<refresh_rate<<std::endl;
    
      END_CODE();
    }
    
    void reset(void) 
    {

      QDPIO::cout<<"MG4DChronoPredictor - Do nothing!"<<std::endl;

    }

    void operator()(T& psi, 
	            const LinearOperator<T>& A, 
	            const T& chi)
    {
    
     QDPIO::cout<<"MG4DChronoPredictor - Do nothing!"<<std::endl;

    }

    virtual void newVector(const T& psi)
    { 

     QDPIO::cout<<"MG4DChronoPredictor - Do nothing!"<<std::endl;

    }


  private:

    std::string subspace_name;
    int refresh_rate;

  };

  
} // End Namespace Chroma

#endif 
