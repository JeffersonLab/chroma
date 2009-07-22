#include "dilution_scheme_info.h"
#include "xml_help.h"
using namespace std;

namespace Chroma {
  namespace LaphEnv {



   // constructor gets input from XMLReader

DilutionSchemeInfo::DilutionSchemeInfo(XMLReader& xml_rdr)
{
 int time_dil_type, spin_dil_type, eigvec_dil_type;
 try{
    XMLReader xml_in(xml_rdr, "//laph_dilution_scheme");
    dil_in(xml_in,"./time_dilution", time_dil_type );
    dil_in(xml_in,"./spin_dilution", spin_dil_type );
    dil_in(xml_in,"./eigvec_dilution", eigvec_dil_type );
    }
 catch(const string& err){
    QDPIO::cerr << "could not initialize DilutionSchemeInfo from XML input"<<endl;
    QDP_abort(1);}
 assign(spin_dil_type, eigvec_dil_type, time_dil_type);
}


DilutionSchemeInfo::DilutionSchemeInfo(const DilutionSchemeInfo& in)
                   : spinDilutionType(in.spinDilutionType),
                     eigvecDilutionType(in.eigvecDilutionType) {}


DilutionSchemeInfo& DilutionSchemeInfo::operator=(const DilutionSchemeInfo& in)
{
 spinDilutionType=in.spinDilutionType;
 eigvecDilutionType=in.eigvecDilutionType;
 return *this;
}


void DilutionSchemeInfo::assign(int spin_dil_type, int eigvec_dil_type, 
                                int time_dil_type)
{
 try{
    if (time_dil_type!=1) throw "only full time dilution allowed";
    if ((spin_dil_type==-1)||(eigvec_dil_type==-1)) 
       throw "dilution types cannot have value -1";
    if ((spin_dil_type>2)||(spin_dil_type<-2))
       throw "spin dilution type must have value -2, 0, 1, 2";
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
   ||(eigvecDilutionType!=in.eigvecDilutionType)){
    QDPIO::cerr << "DilutionSchemeInfo does not checkEqual...abort"<<endl;
    QDP_abort(1);}
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
    int blocksize=((nBasis-1)/nproj)+1;
    int jb=0, bc=0;
    for (int k=0;k<nBasis;k++){
       projind[k]=jb; bc++;
       if (bc==blocksize){ jb++; bc=0;}}
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
 


string DilutionSchemeInfo::output(int indent) const
{
 string pad(3*indent,' ');
 ostringstream oss;
 oss << pad << "<laph_dilution_scheme>"<<endl;
 oss << pad <<"   <time_dilution> " <<endl
     << dil_out(indent,1)
     << pad << "   </time_dilution>"<<endl;
 oss << pad << "   <spin_dilution> " <<endl
     << dil_out(indent,spinDilutionType,false)
     << pad << "   </spin_dilution>"<<endl;
 oss << pad << "   <eigvec_dilution> " <<endl
     << dil_out(indent,eigvecDilutionType,true) 
     << pad << "   </eigvec_dilution>"<<endl;
 oss << pad << "</laph_dilution_scheme>"<<endl;
 return oss.str();
}


void DilutionSchemeInfo::dil_in(XMLReader& xml_in, const std::string& path, 
                               int& DilType)
{
 string tmp;
 DilType=0;
 try{
    read(xml_in,path+"/dilution_type",tmp);
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
       if (xml_in.count(path+"/number_projectors")==1){
          read(xml_in,path+"/number_projectors",nproj);
          if (nproj>1) DilType=nproj;
          else throw "invalid number of block projectors";}
       }
    else if (tmp=="interlace"){
       DilType=-2;
       int nproj;
       if (xml_in.count(path+"/number_projectors")==1){
          read(xml_in,path+"/number_projectors",nproj);
          if (nproj>1) DilType=-nproj;
          else throw "invalid number of interlace projectors";}
       }
    else{
       throw "invalid read";}
    }
 catch(const string& errstr){
       QDPIO::cerr << "invalid DilutionSchemeInfo read from XML"<<endl;
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
 oss << pad <<"      <dilution_type> "<<dtype<<" </dilution_type>"<<endl;
 if (out_nproj){
    if (DilType>1)
       oss << pad << "      <number_projectors> "<<DilType<<" </number_projectors>"<<endl;
    else if (DilType<1)
       oss << pad << "      <number_projectors> "<<-DilType<<" </number_projectors>"<<endl;
    }
 return oss.str();
}


// *************************************************************
  }
}
