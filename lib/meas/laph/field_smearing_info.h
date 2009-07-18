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
// *       <smearing_scheme>                                         *
// *           <link_smear_type> STOUT </link_smear_type>            *
// *           <link_iterations> 4 </link_iterations>                *
// *           <link_staple_weight>  0.25 </link_staple_weight>      *
// *           <quark_smear_type> LAPH </quark_smear_type>           *
// *           <number_laph_eigvecs> 32 </number_laph_eigvecs>       *
// *       </smearing_scheme>                                        *
// *                                                                 *
// *   Example usage:                                                *
// *                                                                 *
// *     XMLReader xml_in(...);                                      *
// *     FieldSmearingInfo smear(xml_in);                            *
// *                                                                 *
// *     FieldSmearingInfo smear2(....);                             *
// *     smear.check(smear2);  // aborts if smear2 != smear          *
// *     if (smear==smear2) ...  // returns boolean                  *
// *                                                                 *
// *     int ival = smear.getNumberOfLinkIterations();               *
// *     double dval = smear.getLinkStapleWeight();                  *
// *     int jval = getNumberOfLaplacianEigenvectors();              *
// *     string sval = getLinkSmearType();                           *
// *     string ssval = getQuarkSmearType();                         *
// *                                                                 *
// *     string out = smear.output();   // xml output                *
// *                                                                 *
// *******************************************************************

class FieldSmearingInfo
{

  int linkIterations;
  double linkStapleWeight;
  int laphNumEigvecs;

 public:  

  FieldSmearingInfo(XMLReader& xml_in);

  FieldSmearingInfo(const FieldSmearingInfo& in) 
     : linkIterations(in.linkIterations),
       linkStapleWeight(in.linkStapleWeight),
       laphNumEigvecs(in.laphNumEigvecs) {}

  FieldSmearingInfo& operator=(const FieldSmearingInfo& in)
   {
    linkIterations=in.linkIterations;
    linkStapleWeight=in.linkStapleWeight;
    laphNumEigvecs=in.laphNumEigvecs;
    return *this;
   }

  void check(const FieldSmearingInfo& in) const
   {
    if  ((linkIterations!=in.linkIterations)
       ||(abs(linkStapleWeight-in.linkStapleWeight)>1e-12)
       ||(laphNumEigvecs!=in.laphNumEigvecs)){
       QDPIO::cerr << "FieldSmearingInfo does not check...abort"<<endl;
       QDP_abort(1);}
   }

  bool operator==(const FieldSmearingInfo& in) const
   {
    return ((linkIterations==in.linkIterations)
          &&(abs(linkStapleWeight-in.linkStapleWeight)<1e-12)
          &&(laphNumEigvecs==in.laphNumEigvecs));
   }


    // output functions

  int getNumberOfLinkIterations() const { return linkIterations; }
  Double getLinkStapleWeight() const { return linkStapleWeight; }
  int getNumberOfLaplacianEigenvectors() const {return laphNumEigvecs;}
  std::string getLinkSmearType() const { return "STOUT"; }
  std::string getQuarkSmearType() const { return "LAPH"; }

  std::string output() const;

};



// **************************************************
  }
}
#endif
