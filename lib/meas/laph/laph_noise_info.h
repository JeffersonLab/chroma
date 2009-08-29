#ifndef LAPH_NOISE_H
#define LAPH_NOISE_H

#include "qdp.h"
#include "chromabase.h"
#include "gauge_configuration_info.h"

namespace Chroma {
  namespace LaphEnv {


// *******************************************************************
// *                                                                 *
// *   An object of class "LaphNoiseInfo" stores identifying info    *
// *   about a quark source noise in the Laph subspace.  The         *
// *   notion of a noise identity extends to an entire ensemble.     *
// *   The XML input must have the format                            *
// *                                                                 *
// *         <LaphNoise>                                             *
// *            <ZNGroup> 4 </ZNGroup>                               *
// *            <seed0> 3156 </seed0>                                *
// *            <seed1> 2981 </seed1>                                *
// *            <seed2> 4013 </seed2>                                *
// *            <seed3> 1132 </seed3>                                *
// *         </LaphNoise>                                            *
// *                                                                 *
// *   The group Z(N) is used, so the value of "N" must be specified *
// *   in the tag named "ZNGroup".  This must be an integer >=2.     *
// *   The value of 4 is recommended.                                *
// *                                                                 *
// *   "seed0", "seed1", "seed2", "seed3" are integers; each must be *
// *   given a value in the range 0 to 4095 (2^12-1), except "seed3" *
// *   whose range is from 0 to 2047 (2^11-1).                       *
// *                                                                 *
// *   A noise vector for a configuration is constructed using a     *
// *   random number generator with a seed.  The four seed integers  *
// *   in a "LaphNoiseInfo" object are used to compute a seed for    *
// *   a configuration identified by a given RHMC trajectory         *
// *   number.  The four integers are combined to compute a QDP++    *
// *   Seed, which is used for the configuration having trajectory   *
// *   number equal to "trajOffset" (a pre-set number, usually 1000).*
// *   For configurations having other trajectory numbers, a simple  *
// *   linear congruential RNG is called some prescribed number of   *
// *   times to obtain the Seed corresponding to that trajectory     *
// *   number.  For trajectory number "k", the number of             *
// *   applications of the RNG is                                    *
// *      "nhits"*("k"-"traj_offset")  for "k">"traj_offset"         *
// *   or                                                            *
// *     ("nhits"+1)*("traj_offset"-"k") for "k"<"traj_offset"       *
// *   where "nhits" has a pre-set value (usually, 1).               *
// *                                                                 *
// *   Example usage:                                                *
// *                                                                 *
// *     XMLReader xml_in(...);                                      *
// *     LaphNoiseInfo rho(xml_in);                                  *
// *                                                                 *
// *     LaphNoiseInfo rho2(....);                                   *
// *     rho.checkEqual(rho2);   // throws exception if rho2 != rho  *
// *     if (rho==rho2) ...      // returns boolean                  *
// *                                                                 *
// *     GaugeConfigurationInfo G(....)                              *
// *     Seed s = rho.getSeed(G);                                    *
// *                                                                 *
// *     string str = rho.output();   // xml string                  *
// *                                                                 *
// *******************************************************************

class LaphNoiseInfo
{

  int znGroup;
  int s0,s1,s2,s3;

  mutable int currTraj;
  mutable Seed rngseed;

 public:  

  LaphNoiseInfo(XMLReader& xml_in);

  LaphNoiseInfo(const std::string& header);

  LaphNoiseInfo(const LaphNoiseInfo& in);

  LaphNoiseInfo& operator=(const LaphNoiseInfo& in);

  ~LaphNoiseInfo(){}

  void checkEqual(const LaphNoiseInfo& in) const;

  bool operator==(const LaphNoiseInfo& in) const;

  bool operator<(const LaphNoiseInfo& in) const;


    // output functions

  Seed getSeed(const GaugeConfigurationInfo& G) const;

  int getZNGroup() const { return znGroup; }
  int getSeed0() const { return s0; }
  int getSeed1() const { return s1; }
  int getSeed2() const { return s2; }
  int getSeed3() const { return s3; }

  std::string output(int indent = 0) const;

  void output(XMLWriter& xmlout) const;

  std::string getHeader() const { return output(0);}

  void getHeader(XMLWriter& xmlout) const { output(xmlout);}

  void binaryWrite(BinaryWriter& out) const;
  void binaryRead(BinaryReader& in);

 private:

  void extract_info_from_reader(XMLReader& xml_in);
  void check_assignment();
  void calc_seed() const;
  void lcongr(int& i0, int& i1, int& i2, int& i3) const;


};



// **************************************************
  }
}
#endif
