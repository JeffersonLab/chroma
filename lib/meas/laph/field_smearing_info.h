#ifndef FIELD_SMEARING_H
#define FIELD_SMEARING_H

#include "qdp.h"
#include "chromabase.h"

namespace Chroma {
  namespace LaphEnv {


// *******************************************************************
// *                                                                 *
// *   Objects of class "FieldSmearingInfo" store identifying info   *
// *   about both the link and quark field smearing used.  Stout     *
// *   link smearing and quark Laph smearing are enforced.  The XML  *
// *   input must have the format                                    *
// *                                                                 *
// *       <stout_laph_smearing>                                     *
// *           <link_iterations> 4 </link_iterations>                *
// *           <link_staple_weight>  0.25 </link_staple_weight>      *
// *           <laph_sigma_cutoff> 0.76 </laph_sigma_cutoff>         *
// *           <number_laph_eigvecs> 32 </number_laph_eigvecs>       *
// *       </stout_laph_smearing>                                    *
// *                                                                 *
// *   Example usage:                                                *
// *                                                                 *
// *     XMLReader xml_in(...);                                      *
// *     FieldSmearingInfo smear(xml_in);                            *
// *                                                                 *
// *     FieldSmearingInfo smear2(....);                             *
// *     smear.checkEqual(smear2);  // throws string exception       *
// *                                //   if smear2 != smear          *
// *     if (smear==smear2) ...     // returns boolean               *
// *                                                                 *
// *     int ival = smear.getNumberOfLinkIterations();               *
// *     double dval = smear.getLinkStapleWeight();                  *
// *     int jval = getNumberOfLaplacianEigenvectors();              *
// *     double dval = smear.getLaphSigmaCutoff();                   *
// *     string sval = getLinkSmearType();                           *
// *     string ssval = getQuarkSmearType();                         *
// *                                                                 *
// *     string out = smear.output();    // xml output               *
// *     string out = smear.output(2);   // indented xml output      *
// *                                                                 *
// *******************************************************************


class FieldSmearingInfo
{

  int linkIterations;
  double linkStapleWeight;
  int laphNumEigvecs;
  double laphSigma;

 public:  

  FieldSmearingInfo(XMLReader& xml_in);

  FieldSmearingInfo(const FieldSmearingInfo& in);

  FieldSmearingInfo& operator=(const FieldSmearingInfo& in);

  void assign(int link_it, double link_wt, int quark_nvecs, double laph_sigma);

  ~FieldSmearingInfo(){}

  void checkEqual(const FieldSmearingInfo& in) const;

  bool operator==(const FieldSmearingInfo& in) const;


    // output functions

  int getNumberOfLinkIterations() const { return linkIterations; }

  double getLinkStapleWeight() const { return linkStapleWeight; }

  int getNumberOfLaplacianEigenvectors() const {return laphNumEigvecs;}
  
  double getLaphSigmaCutoff() const { return laphSigma; }

  std::string getLinkSmearType() const { return "STOUT"; }

  std::string getQuarkSmearType() const { return "LAPH"; }

  std::string output(int indent = 0) const;


};



// **************************************************
  }
}
#endif
