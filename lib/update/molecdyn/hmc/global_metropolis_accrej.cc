// $Id: global_metropolis_accrej.cc,v 2.1 2006-02-02 18:44:20 edwards Exp $
/*! \file
 * \brief Simple metropolis accept/reject
 */

#include "chromabase.h"
#include "update/molecdyn/hmc/global_metropolis_accrej.h"

namespace Chroma { 

  bool globalMetropolisAcceptReject(const Double& DeltaH) { 
    // If deltaH is negative then always accept
    bool ret_val;
    
    
    if ( toBool( DeltaH <= Double(0)) ) {
      ret_val = true;
    }
    else {
      Double AccProb = exp(-DeltaH);
      Double uni_dev;
      random(uni_dev);
      
      
      if( toBool( uni_dev <= AccProb ) ) { 
	
	ret_val = true;
	
      }
      else {    
	
	ret_val = false;
      }
    }
    
    return ret_val;
  }  

}; // End namespace
