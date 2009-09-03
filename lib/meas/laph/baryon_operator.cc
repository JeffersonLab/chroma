#include "baryon_operator.h"
#include "xml_help.h"

using namespace std;

namespace Chroma {
  namespace LaphEnv {



void BaryonOperator::initialize(const string& baryonCoefsTopDirectory)
{
 coefsTopDirectory=tidyString(baryonCoefsTopDirectory);
 if (*coefsTopDirectory.rbegin()!='/') coefsTopDirectory+="/";
 string checkfile=coefsTopDirectory+"nucleon_uud/mom_ray_000/G1g_1/SS_0";

 if (fileExists(checkfile)){
    coefsFlag=true;}
 else{
    QDPIO::cerr << "Baryon coefficient top directory not found"<<endl;
    QDP_abort(1);}
}
  
   // "current" element of xml_in should be <BaryonOp>

BaryonOperator::BaryonOperator(XMLReader& xml_in)
{
 if (!coefsFlag){
    QDPIO::cerr << "Must call BaryonOperator::Initialize(...) before creating baryon objects"<<endl;
    QDP_abort(1);}

 if (!checkXMLCurrentTagName(xml_in,"BaryonOp")){
    QDPIO::cerr << "XML input to BaryonOperator constructor must be"<<endl
                << " enclosed by single <BaryonOp> ... </BaryonOp>"<<endl;
    QDP_abort(1);}

 xmlread(xml_in, "IsospinName", isospinName, "BaryonOperator");
 isospinName=tidyString(isospinName);

 xmlread(xml_in, "Flavor", flavor, "BaryonOperator");
 flavor=tidyString(flavor);

 multi1d<int> p;
 xmlread(xml_in, "Momentum", p, "BaryonOperator");
 if (p.size()!=3){
    QDPIO::cerr << "Bad baryon momentum"<<endl;
    QDP_abort(1);}
 mom.x=p[0];
 mom.y=p[1];
 mom.z=p[2];

 xmlread(xml_in, "Irrep", irrep, "BaryonOperator");
 irrep=tidyString(irrep);

 xmlread(xml_in, "SpatialType", spatialType, "BaryonOperator");
 spatialType=tidyString(spatialType);

 xmlread(xml_in, "IrrepRow", irrepRow, "BaryonOperator");
 xmlread(xml_in, "SpatialIdNum", spatialIdNum, "BaryonOperator");
 xmlread(xml_in, "DispLength", dispLength, "BaryonOperator");
 if ((dispLength<0)||(spatialIdNum<0)||(irrepRow<1)){
    QDPIO::cerr << "Bad baryon operator input xml data"<<endl;
    QDP_abort(1);}

    // check if a single-site operator, then make the displacement length zero
 if (spatialType=="SS") dispLength=0;

    // first, determine the momentum "ray"

 string momRay;
 if (!getMomentumRay(mom.x,mom.y,mom.z,momRay)){
    QDPIO::cerr << "Baryon momentum is not part of allowed ray"<<endl;
    QDP_abort(1);}

    // now look for the baryon operator in the XML coefficient file

 ostringstream oss;
 oss << coefsTopDirectory << isospinName << "_" << flavor
     << "/" << momRay << "/" << irrep << "_" << irrepRow
     << "/" << spatialType << "_" << spatialIdNum;
 string coefs_file=oss.str();

 if (!fileExists(coefs_file)){
    QDPIO::cerr << "invalid baryon operator requested"<<endl;
    QDP_abort(1);}

 TextFileReader in(coefs_file);

 int nterms;
 Double re,im;
 int s1,s2,s3,d1,d2,d3;
 in >> nterms;
 for (int k=0;k<nterms;k++){
    in >> s1 >> s2 >> s3 >> d1 >> d2 >> d3 >> re >> im;
    terms.push_back(ElementalTerm(Elemental(s1,s2,s3,d1,d2,d3),
                    cmplx(re,im)));}

 if (in.fail()){
    in.close();
    QDPIO::cerr << "error occurred during baryon operator read"<<endl;
    QDP_abort(1);}

 in.close();
}


   //   Allowed momentum rays:
   //     000  +00  0+0  00+  ++0  +-0  +0+  +0-  0++  0+-
   //     +++  ++-  +-+  +--

bool BaryonOperator::getMomentumRay(int px, int py, int pz, string& ray) const
{
 ray.clear();
 int ppx=(px<0)?-px:px;  // absolute values
 int ppy=(py<0)?-py:py;
 int ppz=(pz<0)?-pz:pz;
 int n=0,nz[3];
 if (ppx>0) nz[n++]=ppx;
 if (ppy>0) nz[n++]=ppy;
 if (ppz>0) nz[n++]=ppz;
 if ((n==2)&&(nz[0]!=nz[1])) return false;
 if ((n==3)&&((nz[0]!=nz[1])||(nz[1]!=nz[2]))) return false;

 bool sflip=(px<0)||((px==0)&&((py<0)||((py==0)&&(pz<0))));
 char dir[3]={'-','0','+'};
 ppx=(px>0)?1:((px<0)?-1:0);  // get zero or sign
 ppy=(py>0)?1:((py<0)?-1:0);
 ppz=(pz>0)?1:((pz<0)?-1:0);
 if (sflip){ ppx=-ppx; ppy=-ppy; ppz=-ppz;}
 ppx+=1; ppy+=1; ppz+=1;
 ray="mom_ray_"; ray+=dir[ppx]; ray+=dir[ppy]; ray+=dir[ppz];
 return true;
}

string BaryonOperator::output() const
{
 ostringstream oss;
 oss << "<BaryonOp>"<<endl;
 oss << "  <IsospinName> " << isospinName << " </IsospinName>"<<endl;
 oss << "  <Flavor> "<<flavor<<" </Flavor>"<<endl;
 oss << "  <Momentum> " <<mom.x<<" "<<mom.y
      <<" "<<mom.z<<" </Momentum>"<<endl;
 oss << "  <Irrep> "<<irrep<<" </Irrep>"<<endl;
 oss << "  <IrrepRow> "<<irrepRow<<" </IrrepRow>"<<endl;
 oss << "  <SpatialType> "<<spatialType<<" </SpatialType>"<<endl;
 oss << "  <SpatialIdNum> "<<spatialIdNum
      <<" </SpatialIdNum>"<<endl;
 oss << "  <DispLength> "<<dispLength
      <<" </DispLength>"<<endl;
 oss << "</BaryonOp>"<<endl;

 return oss.str();
}

void BaryonOperator::output(XMLWriter& xmlout) const
{
 push(xmlout,"BaryonOp");
 write(xmlout,"IsospinName",isospinName);
 write(xmlout,"Flavor",flavor);
 multi1d<int> p(3);  p[0]=mom.x; p[1]=mom.y; p[2]=mom.z;
 write(xmlout,"Momentum",p);
 write(xmlout,"Irrep",irrep);
 write(xmlout,"IrrepRow",irrepRow);
 write(xmlout,"SpatialType",spatialType);
 write(xmlout,"SpatialIdNum",spatialIdNum);
 write(xmlout,"DispLength",dispLength);
 pop(xmlout);
}

string BaryonOperator::fullOutput() const
{
 ostringstream oss;
 oss << "<BaryonOp>"<<endl;
 oss << "  <IsospinName> " << isospinName << " </IsospinName>"<<endl;
 oss << "  <Flavor> "<<flavor<<" </Flavor>"<<endl;
 oss << "  <Momentum> " <<mom.x<<" "<<mom.y
      <<" "<<mom.z<<" </Momentum>"<<endl;
 oss << "  <Irrep> "<<irrep<<" </Irrep>"<<endl;
 oss << "  <IrrepRow> "<<irrepRow<<" </IrrepRow>"<<endl;
 oss << "  <SpatialType> "<<spatialType<<" </SpatialType>"<<endl;
 oss << "  <SpatialIdNum> "<<spatialIdNum
      <<" </SpatialIdNum>"<<endl;
 oss << "  <DispLength> "<<dispLength
      <<" </DispLength>"<<endl;
 oss << "  <Projection>"<<endl;
 for (BaryonOperator::ElementalIterator it=terms.begin();it!=terms.end();it++){
    oss << "     <Term> "<<it->el.output()
           << " coef=cmplx"<<it->coef<<" </Term>"<<endl;
    }
 oss << "  </Projection>"<<endl;
 oss << "</BaryonOp>"<<endl;

 return oss.str();
}

    //  static data members initialization

string BaryonOperator::coefsTopDirectory;
bool BaryonOperator::coefsFlag=false;



   //  useful routine for reading multiple BaryonOperators;
   //  also reads name of coefficient file

void createBaryonOperators(XMLReader& xml_in, 
                           list<BaryonOperator>& bops)
{
 if (xml_tag_count(xml_in,"BaryonOperators")!=1){
    QDPIO::cerr << "XML input to CreateBaryonOperators must have"<<endl
                << " a single <BaryonOperators> ... </BaryonOperators>"<<endl;
    QDP_abort(1);}
 bops.clear();
 string coef_dir;
 xmlread(xml_in,"ElementalCoefficientDirectory",coef_dir,"createBaryonOperators");
 BaryonOperator::initialize(tidyString(coef_dir));

 int num_bop=xml_tag_count(xml_in,"BaryonOperators/BaryonOp");
 for (int k=1;k<=num_bop;k++){
    ostringstream path;
    path << "./descendant-or-self::BaryonOperators/BaryonOp["<<k<<"]";
    XMLReader xml_op(xml_in,path.str());
    bops.push_back(BaryonOperator(xml_op));}
}

void createBaryonOperators(XMLReader& xml_in, 
                           vector<BaryonOperator>& bops)
{
 list<BaryonOperator> blist;
 createBaryonOperators(xml_in,blist);
 bops.assign(blist.begin(),blist.end());
}


// **************************************************************
  }
}
