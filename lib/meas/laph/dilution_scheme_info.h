#ifndef LAPH_DILUTION_SCHEME_H
#define LAPH_DILUTION_SCHEME_H

#include "qdp.h"
#include "chromabase.h"
#include "xml_help.h"
#include "field_smearing_info.h"
#include <vector>

namespace Chroma {
  namespace LaphEnv {


// *******************************************************************
// *                                                                 *
// *   Objects of class "DilutionSchemeInfo" store identifying info  *
// *   about a Laph dilution scheme.  Full time dilution is          *
// *   enforced.  XML input format for a Laph dilution scheme:       *
// *                                                                 *
// *     <laph_dilution_scheme>                                      *
// *        <time_dilution>                                          *
// *           <dilution_type> full </dilution_type>                 *
// *        </time_dilution>                                         *
// *        <spin_dilution>                                          *
// *           <dilution_type> none </dilution_type>                 *
// *        </spin_dilution>                                         *
// *        <eigvec_dilution>                                        *
// *           <dilution_type> block </dilution_type>                *
// *           <number_projectors> 4 </number_projectors>            *
// *        </eigvec_dilution>                                       *
// *     </laph_dilution_scheme>                                     *
// *                                                                 *
// *   Time dilution must be input as "full".  Spin and Laph         *
// *   eigenvector dilution can be "full", "none", "block", or       *
// *   "interlace".  If "block" or "interlace" is specified, a tag   *
// *   specifying the number of projectors must be included (except  *
// *   for spin dilution, since 2 will be assumed).                  *
// *                                                                 *
// *   The members most likely to be used are                        *
// *                                                                 *
// *     XMLReader xlm_in(...);                                      *
// *     DilutionSchemeInfo dil(xml_in);  // construct from xml      *
// *                                                                 *
// *     DilutionSchemeInfo dil2(...);                               *
// *     dil.checkEqual(dil2);  // if dil != dil2 throw exception    *
// *                                                                 *
// *     string xml_out = dil.output();     // xml output            *
// *                                                                 *
// *     FieldSmearingInfo S(...);                                   *
// *     vector<DilutionSchemeInfo::Projector> dilProjs;             *
// *     dil.getProjectors(S, dilProjs);                             *
// *                                                                 *
// *     nProjs = dilProjs.size();  // number of projectors          *
// *                                                                 *
// *     dilProjs[k].onSpinIndices();   // list of "on" spin indices *
// *                                    //   for projector "k"       *
// *     dilProjs[k].onEigvecIndices(); // list of "on" eigenvector  *
// *                                    //   indices for proj "k"    *
// *                                                                 *
// *                                                                 *
// *******************************************************************



class DilutionSchemeInfo
{
                           // full time dilution is enforced!!

   int spinDilutionType;     //  0 = none, 1 = full
   int eigvecDilutionType;   //  x  (x>=2) = block with no. projectors x
                             // -y  (y>=2) = interlace with no. projectors y 

   mutable vector<list<int> > spinProjs;     // holds the spin projectors
   mutable vector<list<int> > eigvecProjs;   // holds the eigvec projectors


 public:  

   class Projector
   {
      list<int>* spin_components_on;
      list<int>* eigvecs_on;

    public:

      const list<int>& onSpinIndices() const { return *spin_components_on;}
      const list<int>& onEigvecIndices() const { return *eigvecs_on; }

      friend class DilutionSchemeInfo;
   };


  DilutionSchemeInfo(XMLReader& xml_in);

  DilutionSchemeInfo(const DilutionSchemeInfo& in);

  DilutionSchemeInfo& operator=(const DilutionSchemeInfo& in);

  ~DilutionSchemeInfo(){}

  void assign(int spin_dil_type, int eigvec_dil_type, int time_dil_type = 1);

  void checkEqual(const DilutionSchemeInfo& in) const;

  bool operator==(const DilutionSchemeInfo& in) const;


    // output functions

  void getProjectors(const FieldSmearingInfo& S, vector<Projector>& dilprojs) const;
  
  int getNumberOfProjectors(const FieldSmearingInfo& S) const;

  std::string output(int indent = 0) const;   // XML output


 private:

  void dil_in(XMLReader& xml_in, const std::string& path, int& DilType);

  std::string dil_out(int indent, int DilType, bool out_nproj = false) const;

  void setProjectorMasks(vector<list<int> >& projs, int dil_type, int nBasis) const;

  int findProjectorNumber(int dil_type, int nBasis) const;

};



// **************************************************
  }
}
#endif
