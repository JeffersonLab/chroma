#include "chromabase.h"
#include "update/molecdyn/global_metropolis_accrej.h"


bool globalMetropolisAcceptReject(const Double& DeltaH, XMLWriter& monitor) { 
  // If deltaH is negative then always accept
  bool ret_val;
  
  push(monitor, "MetropolisAcceptRejectTest");
  if ( toBool( DeltaH <= Double(0)) ) {
    ret_val = true;
    write(monitor, "AccProb", Double(1));
    
  }
  else {
    Double AccProb = exp(-DeltaH);
    Double uni_dev;
    random(uni_dev);
    
    write(monitor, "AccProb", AccProb);
    write(monitor, "random", uni_dev);
    
    if( toBool( uni_dev <= AccProb ) ) { 
      
      ret_val = true;
      
    }
    else {    
      
      ret_val = false;
    }
  }
  
  write(monitor, "AcceptState", ret_val);
  pop(monitor);
  
  return ret_val;
}  
