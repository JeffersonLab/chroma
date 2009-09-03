#ifndef LAPH_BARYON_OPERATOR_H
#define LAPH_BARYON_OPERATOR_H

#include "qdp.h"
#include "chromabase.h"
#include <vector>
#include <list>

namespace Chroma {
  namespace LaphEnv {


// *******************************************************************
// *                                                                 *
// *   Objects of class "BaryonOperator" store identifying info      *
// *   about one particular baryon operator.  Each baryon is a       *
// *   linear superposition of so-called elemental operators.        *
// *   The Chroma tasks that compute the baryon source/sink          *
// *   functions need to read these superposition coefficients.      *
// *   The XML input into the Chroma tasks must inform Chroma        *
// *   which baryon operators to compute, so a pre-defined input     *
// *   format must be followed.  The needed input format and the     *
// *   required storage format of the elementals in each baryon      *
// *   are described below.                                          *
// *                                                                 *
// *   Chroma input format for each baryon operator:                 *
// *                                                                 *
// *       <BaryonOp>                                                *
// *          <IsospinName> nucleon </IsospinName>                   *
// *          <Flavor> uud </Flavor>                                 *
// *          <Momentum>  0 1 -1  </Momentum>                        *
// *          <Irrep> Hu </Irrep>                                    *
// *          <IrrepRow> 3 </IrrepRow>                               *
// *          <SpatialType> DDL </SpatialType>                       *
// *          <SpatialIdNum> 4 </SpatialIdNum>                       *
// *          <DispLength> 3 </DispLength>                           *
// *       </BaryonOp>                                               *
// *                                                                 *
// *   The superposition coefficients have been pre-computed         *
// *   (using a Maple program) and stored in a pre-defined           *
// *   format in a file located in a particular directory            *
// *   structure.  The elementals and coefficients for any           *
// *   given baryon are stored in a single file -- one file          *
// *   for each baryon operator.  The above requested baryon         *
// *   would be found in the file                                    *
// *                                                                 *
// *     <top_directory>/nucleon_uud/mom_ray_0+-/Hu_3/DDL_4          *
// *                                                                 *
// *   The class needs access to the files which describe the        *
// *   baryon operators in terms of elementals.  Hence, the class    *
// *   needs to know the top directory above. You must call          *                                      *
// *                                                                 *
// *       BaryonOperator::initialize(coefsTopDirectory);            *
// *                                                                 *
// *   first before using this class.  There is no default           *
// *   constructor for this class.  The constructor takes the        *
// *   above XML (in the form of an XMLReader object) as input.      *
// *   Example usage:                                                *
// *                                                                 *
// *       BaryonOperator::initialize("coefsTopDirectory");          *
// *       XMLReader xml_in;   // assigned somehow                   *
// *       BaryonOperator B(xml_in);                                 *
// *                                                                 *
// *   The constructor does some checks on the input, then attempts  *
// *   to read the coefficients in terms of elementals for that      *
// *   operator.  If found, the object is assigned and the           *
// *   expansion details and coefficients in terms of elementals     *
// *   is stored.  If not found, an error message is output and      *
// *   an abort is issued.                                           *
// *                                                                 *
// *   READING multiple baryon operators:                            *
// *                                                                 *
// *   Chroma input format:                                          *
// *                                                                 *
// *      <BaryonOperators>                                          *
// *         <ElementalCoefficientDirectory>                         *
// *             coefsTopDirectory                                   *
// *         </ElementalCoefficientDirectory>                        *
// *           <BaryonOp>                                            *
// *            ...                                                  *
// *           </BaryonOp>                                           *
// *           <BaryonOp>                                            *
// *            ...                                                  *
// *           </BaryonOp>                                           *
// *         ...                                                     *
// *      </BaryonOperators>                                         *
// *                                                                 *
// *   Given the above XML input, one just needs to issue the        *
// *   following:                                                    *
// *                                                                 *
// *      XMLReader xml_in;  // assigned somehow                     *
// *      list<BaryonOperator> blist;                                *
// *      CreateBaryonOperators(xml_in,blist);                       *
// *   or                                                            *
// *      vector<BaryonOperator> bvec;                               *
// *      CreateBaryonOperators(xml_in,bvec);                        *
// *                                                                 *
// *   The above code initializes the class, then reads all of       *
// *   baryon operators inside the <BaryonOperators> tags in         *
// *   the input XML file, storing the results in either a list or   *
// *   a multi1d.  The above functions abort execution if the read   *
// *   of any one baryon operator fails.                             *
// *                                                                 *
// *******************************************************************

class BaryonOperator
{

 public:  

  class Elemental
   {
    int spin1,spin2,spin3;
    int disp1,disp2,disp3;
    int sortIndex;

   public:

    Elemental(int s1, int s2, int s3, int d1, int d2, int d3)
     {if ((s1<1)||(s1>4)||(s2<1)||(s2>4)||(s3<1)||(s3>4)
         ||(d1<-3)||(d1>3)||(d2<-3)||(d2>3)||(d3<-3)||(d3>3)){
         QDPIO::cerr << "invalid BaryonOperator::Elemental"<<endl;
         QDP_abort(1);}
      spin1=s1; spin2=s2; spin3=s3; disp1=d1; disp2=d2; disp3=d3;
      sortIndex=((((((((((spin1-1)<<2)+spin2-1)<<2)+spin3-1)
                 <<3)+disp1+3)<<3)+disp2+3)<<3)+disp3+3;}

    Elemental(const Elemental& rhs) 
     : spin1(rhs.spin1), spin2(rhs.spin2), spin3(rhs.spin3),
       disp1(rhs.disp1), disp2(rhs.disp2), disp3(rhs.disp3),
       sortIndex(rhs.sortIndex) {}

    Elemental& operator=(const Elemental& rhs)
     {spin1=rhs.spin1; spin2=rhs.spin2; spin3=rhs.spin3;
      disp1=rhs.disp1; disp2=rhs.disp2; disp3=rhs.disp3;
      sortIndex=rhs.sortIndex; return *this;}

    int getSpin1() const {return spin1;}
    int getSpin2() const {return spin2;}
    int getSpin3() const {return spin3;}
    int getDisp1() const {return disp1;}
    int getDisp2() const {return disp2;}
    int getDisp3() const {return disp3;}

    bool operator<(const Elemental& rhs) const
     {return sortIndex<rhs.sortIndex;}

    std::string output() const
     {std::ostringstream oss;
      oss << " q1=("<<spin1<<","<<disp1<<")";
      oss << " q2=("<<spin2<<","<<disp2<<")";
      oss << " q3=("<<spin3<<","<<disp3<<") ";
      return oss.str();}
 
   };

  struct Momentum
   {
    int x,y,z;

    Momentum(){}
    Momentum(int px, int py, int pz) : x(px), y(py), z(pz) {}
    Momentum(const Momentum& rhs) : x(rhs.x), y(rhs.y), z(rhs.z) {}
    Momentum& operator=(const Momentum& rhs)
     {x=rhs.x; y=rhs.y; z=rhs.z; return *this;}
    bool operator<(const Momentum& rhs) const
     {return (x<rhs.x)||((x==rhs.x)&&((y<rhs.y)||((y==rhs.y)&&(z<rhs.z))));}
   };

  struct ElementalTerm
   {
    Elemental el;
    DComplex coef;

    ElementalTerm(const Elemental& inel, const DComplex& incf)
      : el(inel), coef(incf) {}
    ElementalTerm(const ElementalTerm& in)
      : el(in.el), coef(in.coef) {}
    ElementalTerm& operator=(const ElementalTerm& in)
     {el=in.el; coef=in.coef; return *this;}
   };

  struct ColorContractSpatialSum
   {
    Elemental el;
    int dispLength;

    ColorContractSpatialSum(const Elemental& inel, int in_dispLength)
      : el(inel), dispLength(in_dispLength) {}
    ColorContractSpatialSum(const ColorContractSpatialSum& in)
      : el(in.el), dispLength(in.dispLength) {}
    ColorContractSpatialSum& operator=(const ColorContractSpatialSum& in)
      {el=in.el; dispLength=in.dispLength; return *this;}
    bool operator<(const ColorContractSpatialSum& in) const
     {return (dispLength<in.dispLength)
         ||((dispLength==in.dispLength)&&(el<in.el));}
   };

  struct IndexCoef
   {
    int index;
    DComplex coef;

    IndexCoef(int ind, const DComplex& incf) : index(ind), coef(incf) {}
    IndexCoef(const IndexCoef& in) : index(in.index), coef(in.coef) {}
    IndexCoef& operator=(const IndexCoef& in)
     {index=in.index; coef=in.coef; return *this;}
   };

 private:

  static std::string coefsTopDirectory;
  static bool coefsFlag;

  std::string isospinName;
  std::string flavor;
  Momentum mom;
  std::string irrep;
  int irrepRow;
  std::string spatialType;
  int spatialIdNum;
  int dispLength;

  list<ElementalTerm> terms;

 public:

  static void initialize(const std::string& baryon_coefsTopDirectory);

  BaryonOperator(XMLReader& xml_in);

  BaryonOperator(const BaryonOperator& B) 
     : isospinName(B.isospinName),flavor(B.flavor),
       mom(B.mom),irrep(B.irrep),irrepRow(B.irrepRow),
       spatialType(B.spatialType),
       spatialIdNum(B.spatialIdNum),
       dispLength(B.dispLength),
       terms(B.terms) {}

  BaryonOperator& operator=(const BaryonOperator& B)
   {isospinName=B.isospinName; flavor=B.flavor; 
    mom=B.mom; irrep=B.irrep; irrepRow=B.irrepRow; 
    spatialType=B.spatialType; 
    spatialIdNum=B.spatialIdNum; 
    dispLength=B.dispLength;
    terms=B.terms;
    return *this;}

    // output functions

  std::string getIsospinName() const { return isospinName; }

  std::string getFlavor() const { return flavor; }

  Momentum getMomentumValue() const { return mom; }

  multi1d<int> getMomentumMulti1d() const
  {multi1d<int> p(3); p[0]=mom.x; p[1]=mom.y; p[2]=mom.z;
   return p;}

  int getXMomentum() const { return mom.x; }

  int getYMomentum() const { return mom.y; }

  int getZMomentum() const { return mom.z; }

  std::string getIrrep() const { return irrep; }

  int getIrrepRow() const { return irrepRow; }

  std::string getSpatialType() const { return spatialType; }

  int getSpatialIdNumber() const { return spatialIdNum; }

  int getDisplacementLength() const { return dispLength; }

  int getNumberOfElementals() const { return terms.size(); }

  typedef list<ElementalTerm>::const_iterator ElementalIterator;

  ElementalIterator begin() const { return terms.begin(); }

  ElementalIterator end() const { return terms.end(); }

  std::string output() const;  // XML output of baryon tags

  void output(XMLWriter& xmlout) const;  // XML output of baryon tags

  std::string fullOutput() const;  // XML output including elementals

    // returns the flavor std::string, but replaces all 'u' and
    // 'd' characters by 'l'

  std::string getFlavorLS() const
   {std::string tmp(flavor);
    for (int i=0;i<tmp.length();i++)
       if ((tmp[i]=='u')||(tmp[i]=='d')) tmp[i]='l';
    return tmp;}

  bool operator==(const BaryonOperator& rhs)
  {return  (spatialIdNum==rhs.spatialIdNum)
         &&(irrep==rhs.irrep)&&(irrepRow==rhs.irrepRow)
         &&(spatialType==rhs.spatialType)
         &&(isospinName==rhs.isospinName)&&(flavor==rhs.flavor)
         &&(mom.x==rhs.mom.x)&&(mom.y==mom.y)&&(mom.z==mom.z)
         &&(dispLength==rhs.dispLength);}

 private:

  bool getMomentumRay(int px, int py, int pz, std::string& ray) const;

};



// **************************************************


   //  useful routine for reading multiple BaryonOperators

void createBaryonOperators(XMLReader& xml_in, 
                           list<BaryonOperator>& bops);
void createBaryonOperators(XMLReader& xml_in, 
                           vector<BaryonOperator>& bops);


// **************************************************
  }
}
#endif  
