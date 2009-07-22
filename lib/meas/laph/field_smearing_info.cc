#include "field_smearing_info.h"
#include "xml_help.h"
using namespace std;

namespace Chroma {
  namespace LaphEnv {


   // constructor extracts needed info from XMLReader

FieldSmearingInfo::FieldSmearingInfo(const FieldSmearingInfo& in) 
            : linkIterations(in.linkIterations),
              linkStapleWeight(in.linkStapleWeight),
              laphNumEigvecs(in.laphNumEigvecs) {}

FieldSmearingInfo& FieldSmearingInfo::operator=(const FieldSmearingInfo& in)
{
 linkIterations=in.linkIterations;
 linkStapleWeight=in.linkStapleWeight;
 laphNumEigvecs=in.laphNumEigvecs;
 return *this;
}

void FieldSmearingInfo::checkEqual(const FieldSmearingInfo& in) const
{
 if  ((linkIterations!=in.linkIterations)
    ||(abs(linkStapleWeight-in.linkStapleWeight)>1e-12)
    ||(laphNumEigvecs!=in.laphNumEigvecs)){
    QDPIO::cerr << "FieldSmearingInfo does not check...abort"<<endl;
     QDP_abort(1);}
}

bool FieldSmearingInfo::operator==(const FieldSmearingInfo& in) const
{
 return ((linkIterations==in.linkIterations)
       &&(abs(linkStapleWeight-in.linkStapleWeight)<1e-12)
       &&(laphNumEigvecs==in.laphNumEigvecs));
}


   // fatally aborting if not found

FieldSmearingInfo::FieldSmearingInfo(XMLReader& xml_rdr)
{
 int link_it,quark_nvecs;
 double link_wt;
 try{
    XMLReader xml_in(xml_rdr, "//stout_laph_smearing");
    read(xml_in,"./link_iterations", link_it );
    read(xml_in,"./link_staple_weight", link_wt );
    read(xml_in,"./number_laph_eigvecs", quark_nvecs );
    }
 catch(const string& err){
    QDPIO::cerr << "could not initialize FieldSmearingInfo from XML input"<<endl;
    QDP_abort(1);}
 assign(link_it,link_wt,quark_nvecs);
}

void FieldSmearingInfo::assign(int link_it, double link_wt, int quark_nvecs)
{
 linkIterations=link_it;
 linkStapleWeight=link_wt;
 laphNumEigvecs=quark_nvecs;
 if ((linkIterations<0)||(linkStapleWeight<0.0)||(laphNumEigvecs<1)){
    QDPIO::cerr << "invalid smearing scheme"<<endl;
    QDP_abort(1);}
}

string FieldSmearingInfo::output() const
{
 ostringstream oss;
 oss << "<stout_laph_smearing>"<<endl;
 oss << "  <link_iterations> " << linkIterations 
     << " </link_iterations>"<<endl;
 oss << "  <link_staple_weight> " << linkStapleWeight 
     << " </link_staple_weight>"<<endl;
 oss << "  <number_laph_eigvecs> " << laphNumEigvecs 
     << " </number_laph_eigvecs>"<<endl;
 oss << "</stout_laph_smearing>"<<endl;
 return oss.str();
}


// *************************************************************
  }
}
