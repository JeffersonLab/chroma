#ifndef TIME_SLICES_H
#define TIME_SLICES_H

#include "qdp.h" 

namespace Chroma {
  namespace LaphEnv {

// *********************************************************

  // This function can be used to define Timeslices on the lattice.
  //
  // Usage:
  //     Set timeslices;
  //     timeslices.make(TimeSliceFunc(Tdir)); 
  //     Subset slice = timeslices[t];  // how to get a time slice
  //     int nslices = timeslices.numSubsets();  


class TimeSliceFunc : public SetFunc 
{
public:

  TimeSliceFunc(int dir): mu(dir) {}

  // Simply return the mu'th coordinate as the subset index
  // for a 4D coord
  int operator()(const multi1d<int>& coord) const
  {return coord[mu];}

  // The number of subsets is the length of the lattice
  // in direction mu
  int numSubsets() const {return Layout::lattSize()[mu];}
  
private: 
  int mu;// Time direction
};

// *********************************************************
  }
}
#endif

