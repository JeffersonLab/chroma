#include "dilution_scheme_info.h"
#include "xml_help.h"
using namespace std;

namespace Chroma {
  namespace LaphEnv {


// *************************************************************

   // XMLReader constructor

DilutionSchemeInfo::DilutionSchemeInfo(XMLReader& xml_in)
{
 if (xml_tag_count(xml_in,"LaphDilutionScheme")!=1){
    QDPIO::cerr << "Bad XML input to DilutionSchemeInfo"<<endl;
    QDPIO::cerr << "Expected one <LaphDilutionScheme> tag"<<endl;
    QDP_abort(1);}
 XMLReader xmlr(xml_in, "./descendant-or-self::LaphDilutionScheme");
 assign_from_reader(xmlr);
}

void DilutionSchemeInfo::assign_from_reader(XMLReader& xml_in)
{
 int time_dil_type, spin_dil_type, eigvec_dil_type;
 try{
    dil_in(xml_in,"./TimeDilution", time_dil_type );
    dil_in(xml_in,"./SpinDilution", spin_dil_type );
    dil_in(xml_in,"./EigvecDilution", eigvec_dil_type );
    }
 catch(const string& err){
    QDPIO::cerr << "could not initialize DilutionSchemeInfo from XML input"<<endl;
    QDP_abort(1);}
 assign(spin_dil_type, eigvec_dil_type, time_dil_type);
}

 // *************************************************************

   // This version of the constructor assumes that header information
   // from a quark_source_sink file, for example, is passed in.

DilutionSchemeInfo::DilutionSchemeInfo(const string& header)
{
 string dilution_header;
 extract_xml_element(header,"LaphDilutionScheme",dilution_header,
                     "DilutionSchemeInfo");
 stringstream tmp;
 tmp << dilution_header;
 XMLReader xmlr0(tmp);
 XMLReader xmlr(xmlr0,"/LaphDilutionScheme");  
 assign_from_reader(xmlr);
}


  // ************************************************************

    // copy constructor

DilutionSchemeInfo::DilutionSchemeInfo(const DilutionSchemeInfo& in)
                   : spinDilutionType(in.spinDilutionType),
                     eigvecDilutionType(in.eigvecDilutionType) {}


DilutionSchemeInfo& DilutionSchemeInfo::operator=(const DilutionSchemeInfo& in)
{
 spinDilutionType=in.spinDilutionType;
 eigvecDilutionType=in.eigvecDilutionType;
 return *this;
}


// ***************************************************************


void DilutionSchemeInfo::assign(int spin_dil_type, int eigvec_dil_type, 
                                int time_dil_type)
{
 try{
    if (time_dil_type!=1) throw string("only full time dilution allowed");
    if ((spin_dil_type==-1)||(eigvec_dil_type==-1)) 
       throw string("dilution types cannot have value -1");
    if ((spin_dil_type>2)||(spin_dil_type<-2))
       throw string("spin dilution type must have value -2, 0, 1, 2");
    }
 catch(const string& errmsg){
    QDPIO::cerr << "Invalid DilutionSchemeInfo assigment:"<<endl;
    QDPIO::cerr << "   ..."<<errmsg<<endl;
    QDP_abort(1);}

 spinDilutionType=spin_dil_type;
 eigvecDilutionType=eigvec_dil_type;
}


void DilutionSchemeInfo::checkEqual(const DilutionSchemeInfo& in) const
{
 if ((spinDilutionType!=in.spinDilutionType)
   ||(eigvecDilutionType!=in.eigvecDilutionType))
    throw string("DilutionSchemeInfo does not checkEqual...abort");
}

bool DilutionSchemeInfo::operator==(const DilutionSchemeInfo& in) const
{
 return ((spinDilutionType==in.spinDilutionType)
       &&(eigvecDilutionType==in.eigvecDilutionType));
}



void DilutionSchemeInfo::getProjectors(const FieldSmearingInfo& S,
                                       vector<Projector>& dilProjs) const
{
 setProjectorMasks(spinProjs, spinDilutionType, 4);
 int nEig=S.getNumberOfLaplacianEigenvectors();
 setProjectorMasks(eigvecProjs, eigvecDilutionType, nEig);

 int nspinproj=spinProjs.size();
 int neigproj=eigvecProjs.size();
 dilProjs.resize(nspinproj*neigproj);
 int count=0;
 for (int v=0;v<neigproj;v++)
    for (int s=0;s<nspinproj;s++){
       dilProjs[count].spin_components_on=&(spinProjs[s]);
       dilProjs[count].eigvecs_on=&(eigvecProjs[v]);
       count++;}
}


void DilutionSchemeInfo::setProjectorMasks(vector<list<int> >& projs, 
                                           int dil_type, int nBasis) const
{
 projs.clear();
 int nproj;
 vector<int> projind(nBasis);
 if (dil_type==0){         // no dilution
    nproj=1;
    for (int k=0;k<nBasis;k++) projind[k]=0;
    }
 else if (dil_type==1){    // full dilution
    nproj=nBasis;
    for (int k=0;k<nBasis;k++) projind[k]=k;
    }
 else if (dil_type>1){     // block dilution
    nproj=dil_type;
    if (nproj>(nBasis/2)){
       QDPIO::cerr << "number of blocked dilution projectors too large"<<endl;
       QDP_abort(1);}
    int bksize=nBasis/nproj;
    vector<int> blocksize(nproj,bksize);
    int tb=nBasis-bksize*nproj;
    for (int sb=0;sb<tb;sb++)
       blocksize[sb]++;
    int jb=0, bc=0;
    for (int k=0;k<nBasis;k++){
       projind[k]=jb; bc++;
       if (bc==blocksize[jb]){ jb++; bc=0;}}
    }
 else if (dil_type<-1){    // interlace dilution
    nproj=-dil_type;
    if (nproj>(nBasis/2)){
       QDPIO::cerr << "number of interlaced dilution projectors too large"<<endl;
       QDP_abort(1);}
    int jp=0;
    for (int k=0;k<nBasis;k++){
       projind[k]=jp++;
       if (jp==nproj) jp=0;}
    }

 projs.resize(nproj);
 for (int k=0;k<nBasis;k++)
    projs[projind[k]].push_back(k);
}
 

int DilutionSchemeInfo::findProjectorNumber(int dil_type, int nBasis) const
{
 if (dil_type==0) return 1;                 // no dilution
 else if (dil_type==1) return nBasis;       // full dilution
 else if (dil_type>1) return dil_type;      // block dilution
 else if (dil_type<-1) return -dil_type;    // interlace dilution
}
 
int DilutionSchemeInfo::getNumberOfProjectors(const FieldSmearingInfo& S) const
{
 int nspinproj=findProjectorNumber(spinDilutionType, 4);
 int nEig=S.getNumberOfLaplacianEigenvectors();
 int neigproj=findProjectorNumber(eigvecDilutionType, nEig);
 return nspinproj*neigproj;
}


string DilutionSchemeInfo::output(int indent) const
{
 string pad(3*indent,' ');
 ostringstream oss;
 oss << pad << "<LaphDilutionScheme>"<<endl;
 oss << pad <<"   <TimeDilution> " <<endl
     << dil_out(indent,1)
     << pad << "   </TimeDilution>"<<endl;
 oss << pad << "   <SpinDilution> " <<endl
     << dil_out(indent,spinDilutionType,false)
     << pad << "   </SpinDilution>"<<endl;
 oss << pad << "   <EigvecDilution> " <<endl
     << dil_out(indent,eigvecDilutionType,true) 
     << pad << "   </EigvecDilution>"<<endl;
 oss << pad << "</LaphDilutionScheme>"<<endl;
 return oss.str();
}

void DilutionSchemeInfo::output(XMLWriter& xmlout) const
{
 push(xmlout,"LaphDilutionScheme");
 push(xmlout,"TimeDilution");
 dil_out(xmlout,1);
 pop(xmlout);
 push(xmlout,"SpinDilution");
 dil_out(xmlout,spinDilutionType,false);
 pop(xmlout);
 push(xmlout,"EigvecDilution");
 dil_out(xmlout,eigvecDilutionType,true); 
 pop(xmlout);
 pop(xmlout);
}


void DilutionSchemeInfo::dil_in(XMLReader& xml_in, const std::string& path, 
                               int& DilType)
{
 string tmp;
 DilType=0;
 try{
    read(xml_in,path+"/DilutionType",tmp);
    tmp=tidyString(tmp);
    if (tmp=="none"){
       DilType=0;
       }
    else if (tmp=="full"){
       DilType=1;
       }
    else if (tmp=="block"){
       DilType=2;
       int nproj;
       if (xml_in.count(path+"/NumberProjectors")==1){
          read(xml_in,path+"/NumberProjectors",nproj);
          if (nproj>1) DilType=nproj;
          else throw string("invalid number of block projectors");}
       }
    else if (tmp=="interlace"){
       DilType=-2;
       int nproj;
       if (xml_in.count(path+"/NumberProjectors")==1){
          read(xml_in,path+"/NumberProjectors",nproj);
          if (nproj>1) DilType=-nproj;
          else throw string("invalid number of interlace projectors");}
       }
    else{
       throw string("invalid read");}
    }
 catch(const string& errstr){
       QDPIO::cerr << "invalid DilutionSchemeInfo read from XML"<<endl;
       QDPIO::cerr << errstr << endl;
       QDP_abort(1);}
}


string DilutionSchemeInfo::dil_out(int indent, int DilType, 
                                   bool out_nproj) const
{
 string pad(3*indent,' ');
 string dtype;
 if (DilType==0) dtype="none";
 else if (DilType==1) dtype="full";
 else if (DilType>1) dtype="block";
 else if (DilType<1) dtype="interlace";
 ostringstream oss;
 oss << pad <<"      <DilutionType> "<<dtype<<" </DilutionType>"<<endl;
 if (out_nproj){
    if (DilType>1)
       oss << pad << "      <NumberProjectors> "<<DilType<<" </NumberProjectors>"<<endl;
    else if (DilType<1)
       oss << pad << "      <NumberProjectors> "<<-DilType<<" </NumberProjectors>"<<endl;
    }
 return oss.str();
}

void DilutionSchemeInfo::dil_out(XMLWriter& xmlout,
                                 int DilType, bool out_nproj) const
{
 string dtype;
 if (DilType==0) dtype="none";
 else if (DilType==1) dtype="full";
 else if (DilType>1) dtype="block";
 else if (DilType<1) dtype="interlace";
 write(xmlout,"DilutionType",dtype);
 if (out_nproj){
    if (DilType>1)
       write(xmlout,"NumberProjectors",DilType);
    else if (DilType<1)
       write(xmlout,"NumberProjectors",-DilType);
    }
}


// *************************************************************
  }
}
