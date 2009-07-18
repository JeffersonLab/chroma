#ifndef LAPH_NOISE_H
#define LAPH_NOISE_H

#include "qdp.h"
#include "chromabase.h"
#include "laph.h"

namespace Chroma {
  namespace LaphEnv {


// *******************************************************************
// *                                                                 *
// *   An object of class "LaphNoiseInfo" stores identifying info    *
// *   about a quark source noise in the Laph subspace.  The         *
// *   notion of a noise identity extends to an entire ensemble.     *
// *   The XML input must have the format                            *
// *                                                                 *
// *         <laph_noise>                                            *
// *            <seed0> 1234567890 </seed0>                          *
// *            <seed1> 2345678901 </seed1>                          *
// *            <seed2> 3456789012 </seed2>                          *
// *            <seed3> 4012345678 </seed3>                          *
// *            <traj_offset> 1024 </traj_offset>                    *
// *            <nhits> 8 </nhits>                                   *
// *         </laph_noise>                                           *
// *                                                                 *
// *   Each seed is an unsigned integer (0..2^32-1).  Basically,     *
// *   a noise on a configuration is constructed using the random    *
// *   number generator with a seed.  The information contained      *
// *   in a "LaphNoiseInfo" object is used to compute a seed for     *
// *   a configuration identified by a given RHMC trajectory         *
// *   number.  The four integers are combined to compute a QDP++    *
// *   Seed, which is used for the configuration having trajectory   *
// *   number equal to "traj_offset".  For configurations having     *
// *   other trajectory numbers, a simple linear congruential        *
// *   RNG is called some prescribed number of times to obtain       *
// *   the Seed corresponding to that trajectory number.  For        *
// *   trajectory number "k", the number of applications of the      *
// *   RNG is                                                        *
// *      "nhits"*("k"-"traj_offset")  for "k">"traj_offset"         *
// *   or                                                            *
// *     ("nhits"+1)*("traj_offset"-"k") for "k"<"traj_offset"       *
// *                                                                 *
// *                                                                 *
// *   Example usage:                                                *
// *                                                                 *
// *     XMLReader xml_in(...);                                      *
// *     LaphNoiseInfo rho(xml_in);                                  *
// *                                                                 *
// *     LaphNoiseInfo rho2(....);                                   *
// *     rho.check(rho2);        // aborts if rho2 != rho            *
// *     if (rho==rho2) ...      // returns boolean                  *
// *                                                                 *
// *     GaugeConfigurationInfo G(....)                              *
// *     Seed s = rho.getSeed(G);                                    *
// *                                                                 *
// *     int ival = getNumberOfHits();                               *
// *     int jval = getTrajectoryOffset();                           *
// *                                                                 *
// *     string str = rho.output();   // xml string                  *
// *                                                                 *
// *******************************************************************

class LaphNoiseInfo
{

  unsigned int s0,s1,s2,s3;
  int trajOffset;
  unsigned int nHits;
  static const unsigned long long int ran_mult=1103515245ULL;
  static const unsigned long long int ran_shift=12345ULL;
  static const unsigned long long int ran_mod=4294967296ULL;

  mutable int currTraj;
  mutable Seed rngseed;

 public:  

  LaphNoiseInfo(XMLReader& xml_in);

  LaphNoiseInfo(const LaphNoiseInfo& in) 
     :  s0(in.s0), s1(in.s1), s2(in.s2), s3(in.s3),
        trajOffset(in.trajOffset), nHits(in.nHits),
        currTraj(in.currTraj) 
   {
    rngseed=in.rngseed;
   }

  LaphNoiseInfo& operator=(const LaphNoiseInfo& in)
   {
    s0=in.s0; s1=in.s1; s2=in.s2; s3=in.s3;
    trajOffset=in.trajOffset;
    nHits=in.nHits;
    currTraj=in.currTraj;
    rngseed=in.rngseed;
    return *this;
   }

  void check(const LaphNoiseInfo& in) const
   {
    if  ((s0!=in.s0)||(s1!=in.s1)||(s2!=in.s2)||(s3!=in.s3)
       ||(trajOffset!=in.trajOffset)||(nHits!=in.nHits)){
       QDPIO::cerr << "LaphNoiseInfo does not check...abort"<<endl;
       QDP_abort(1);}
   }

  bool operator==(const LaphNoiseInfo& in) const
   {
    return ((s0==in.s0)&&(s1==in.s1)&&(s2==in.s2)&&(s3==in.s3)
       &&(trajOffset==in.trajOffset)&&(nHits==in.nHits));
   }

    // output functions

  Seed getSeed(const GaugeConfigurationInfo& G) const;
  int getNumberOfHits() const { return nHits; }
  int getTrajectoryOffset() const { return trajOffset; }

  std::string output() const;


 private:

  void lcongr(unsigned int& ss0, unsigned int& ss1,
              unsigned int& ss2, unsigned int& ss3) const;


};



// **************************************************
  }
}
#endif
