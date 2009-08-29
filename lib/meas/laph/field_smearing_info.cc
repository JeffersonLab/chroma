#include "field_smearing_info.h"
#include "xml_help.h"
using namespace std;

namespace Chroma {
  namespace LaphEnv {


// *************************************************************

   // XMLReader constructor

FieldSmearingInfo::FieldSmearingInfo(XMLReader& xml_in)
{
 if (xml_tag_count(xml_in,"StoutLaphSmearing")!=1){
    QDPIO::cerr << "Bad XML input to FieldSmearingInfo"<<endl;
    QDPIO::cerr << "Expected one <StoutLaphSmearing> tag"<<endl;
    QDP_abort(1);}
 XMLReader xmlr(xml_in, "./descendant-or-self::StoutLaphSmearing");
 extract_info_from_reader(xmlr);
}


void FieldSmearingInfo::extract_info_from_reader(XMLReader& xml_in)
{
 xmlread(xml_in,"LinkIterations", linkIterations, "FieldSmearingInfo");
 xmlread(xml_in,"LinkStapleWeight", linkStapleWeight, "FieldSmearingInfo");
 xmlread(xml_in,"LaphSigmaCutoff", laphSigma, "FieldSmearingInfo");
 xmlread(xml_in,"NumberLaphEigvecs", laphNumEigvecs, "FieldSmearingInfo");
 if ((linkIterations<0)||(linkStapleWeight<0.0)||(laphNumEigvecs<1)
    ||(laphSigma<=0)){
    QDPIO::cerr << "invalid smearing scheme parameters in FieldSmearingInfo"<<endl;
    QDP_abort(1);}
}


 // *************************************************************

   // This version of the constructor assumes that header information
   // from a quark_source_sink file, for example, is passed in.

FieldSmearingInfo::FieldSmearingInfo(const string& header)
{
 string smearing_header;
 extract_xml_element(header,"StoutLaphSmearing",smearing_header,
                     "FieldSmearingInfo");
 stringstream tmp;
 tmp << smearing_header;
 XMLReader xmlr0(tmp);
 XMLReader xmlr(xmlr0,"/StoutLaphSmearing");  
 extract_info_from_reader(xmlr);
}


  // ************************************************************

    // copy constructor

FieldSmearingInfo::FieldSmearingInfo(const FieldSmearingInfo& in) 
            : linkIterations(in.linkIterations),
              linkStapleWeight(in.linkStapleWeight),
              laphNumEigvecs(in.laphNumEigvecs),
              laphSigma(in.laphSigma) {}

FieldSmearingInfo& FieldSmearingInfo::operator=(const FieldSmearingInfo& in)
{
 linkIterations=in.linkIterations;
 linkStapleWeight=in.linkStapleWeight;
 laphNumEigvecs=in.laphNumEigvecs;
 laphSigma=in.laphSigma;
 return *this;
}

void FieldSmearingInfo::checkEqual(const FieldSmearingInfo& in) const
{
 if  ((linkIterations!=in.linkIterations)
    ||(abs(linkStapleWeight-in.linkStapleWeight)>1e-12)
    ||(abs(laphSigma-in.laphSigma)>1e-12)
    ||(laphNumEigvecs!=in.laphNumEigvecs))
    throw string("FieldSmearingInfo checkEqual failed...");
}

bool FieldSmearingInfo::operator==(const FieldSmearingInfo& in) const
{
 return ((linkIterations==in.linkIterations)
       &&(abs(linkStapleWeight-in.linkStapleWeight)<1e-12)
       &&(abs(laphSigma-in.laphSigma)<1e-12)
       &&(laphNumEigvecs==in.laphNumEigvecs));
}


string FieldSmearingInfo::output(int indent) const
{
 string pad(3*indent,' ');
 ostringstream oss;
 oss << pad << "<StoutLaphSmearing>"<<endl;
 oss << pad << "   <LinkIterations>" << linkIterations 
     << "</LinkIterations>"<<endl;
 oss << pad << "   <LinkStapleWeight>" << linkStapleWeight 
     << "</LinkStapleWeight>"<<endl;
 oss << pad << "   <LaphSigmaCutoff>" << laphSigma 
     << "</LaphSigmaCutoff>"<<endl;
 oss << pad << "   <NumberLaphEigvecs>" << laphNumEigvecs 
     << "</NumberLaphEigvecs>"<<endl;
 oss << pad << "</StoutLaphSmearing>"<<endl;
 return oss.str();
}

void FieldSmearingInfo::output(XMLWriter& xmlout) const
{
 push(xmlout,"StoutLaphSmearing");
 write(xmlout,"LinkIterations",linkIterations);
 write(xmlout,"LinkStapleWeight",linkStapleWeight);
 write(xmlout,"LaphSigmaCutoff",laphSigma);
 write(xmlout,"NumberLaphEigvecs",laphNumEigvecs);
 pop(xmlout);
}


// *************************************************************
  }
}
