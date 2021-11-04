// -*- C++ -*-
/*! \file
 *  \brief Stout utilities
 */

#include "chroma_config.h"
#include "chromabase.h"
#include "util/gauge/hyp_utils.h"
#if QDP_NC != 3
#include "util/gauge/expm12.h"
#endif


namespace Chroma 
{ 

  //! Timings
  /*! \ingroup gauge */
  namespace HypLinkTimings { 
    static double smearing_secs = 0;
    double getSmearingTime() { 
      return smearing_secs;
    }

    static double force_secs = 0;
    double getForceTime() { 
      return force_secs;
    }

    static double functions_secs = 0;
    double getFunctionsTime() { 
      return functions_secs;
    }
  }

  //! Utilities
  /*! \ingroup gauge */
  namespace Hyping 
  {

    //! Do the force recursion from level i+1, to level i
    void deriv_recurse(multi1d<LatticeColorMatrix>&  F,
		       const multi1d<bool>& smear_in_this_dirP,
                       const Real alpha1,
                       const Real alpha2,
                       const Real alpha3,
                       const int BlkMax,
                       const Real BlkAccu,
                       const multi1d<LatticeColorMatrix>& u)
    {
      START_CODE();

      END_CODE();
    }

    /*! \ingroup gauge */
    void smear_links(const multi1d<LatticeColorMatrix>& current, 
                     multi1d<LatticeColorMatrix>& next,
                     const multi1d<bool>& smear_in_this_dirP,
                     const Real alpha1,
                     const Real alpha2,
                     const Real alpha3,
                     const int BlkMax,
                     const Real BlkAccu)
    {
      START_CODE();
      
      for(int mu = 0; mu < Nd; mu++) {
        if( smear_in_this_dirP[mu] ) {
          LatticeColorMatrix Q;
          
#if QDP_NC == 3
          next = current[mu];
#else 
          // Q is in hermitian form. We can simply exponentiate with a 
          // Taylor series for Nc != 3 builds
          //expim20(Q);
          next = current[mu];
#endif
        }
        else { 
          next[mu] = current[mu];  // Unsmeared
        }	
      }
      END_CODE();
    }
  } // End Namespace Hyping
} // End Namespace Chroma
	
    

