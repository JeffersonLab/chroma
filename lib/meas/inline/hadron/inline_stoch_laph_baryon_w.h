// -*- C++ -*-
// $Id: inline_stoch_laph_baryon_w.h,v 3.4 2009-07-15 02:52:08 jbulava Exp $
/*! \file
 * \brief Inline measurement of stochastic source and sink functions 
 * for baryons
 */

#ifndef __inline_stoch_laph_baryon_h__
#define __inline_stoch_laph_baryon_h__

#include "meas/inline/abs_inline_measurement.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{ 
	/*! \ingroup inlinehadron */
	namespace InlineStochLapHBaryonEnv 
	{
		extern const std::string name;
		bool registerAll();


		//! Parameter structure
		/*! \ingroup inlinehadron */
		struct Params
		{
			Params();
			Params(XMLReader& xml_in, const std::string& path);
			void writeXML(XMLWriter& xml_out, const std::string& path);

			unsigned long      frequency;

			//Noise ids
			//T0's 
			//cfg info
			//dilution scheme
			//smearing info (link and quark)

			struct Param_t
			{

				multi1d<BaryonOperator> bops;               /*!< Baryon operator 
																											information 
																											(see baryon_operator.h)*/

				GroupXML_t          link_smearing;         /*!< link smearing xml */

				int nnoises;         //The number of noise source to be used (must be 
														 //>= 3 
				
				std::string ud_mass_label; //The key for the light quark mass
				std::string s_mass_label; //The key for the strange quark mass
			
				std::string coeff_top_dir; //The location of the operator coefficients

			};

			struct NamedObject_t
			{
		
				//Names of files that were output of STOCH_LAPH_QUARK
		
				multi1d<std::string> quark_files; //The db files where the diluted
																					//sources and sinks are stored
																					
				std::string baryon_file; //The output db file that stores the 
																 //Operator source and sink functions.
				
				std::string gauge_id;
			
			};

			Param_t        param;      /*!< Parameters */    
			NamedObject_t  named_obj;  /*!< Named objects */
		
		};


		//! Inline measurement of stochastic group baryon operators
		/*! \ingroup inlinehadron */
		class InlineMeas : public AbsInlineMeasurement 
		{
			public:
				~InlineMeas() {}
				InlineMeas(const Params& p) : params(p) {}
				InlineMeas(const InlineMeas& p) : params(p.params) {}

				unsigned long getFrequency(void) const {return params.frequency;}

				//! Do the measurement
				void operator()(const unsigned long update_no,
						XMLWriter& xml_out); 

			protected:
				//! Do the measurement
				void func(const unsigned long update_no,
						XMLWriter& xml_out); 

			private:
				Params params;
		};


		struct Sink_qqq_t
		{
			
			Sink_qqq_t(int dil_size, int nt) : dilutions.resize(dil_size, 
					dil_size, dil_size) 
			{
				

				for (int d1 = 0 ; d1 < dil_size ; ++d1)
					for (int d2 = 0 ; d2 < dil_size ; ++d2)
						for (int d3 = 0 ; d3 < dil_size ; ++d3)
						{
							dilutions(d1, d2, d3).time.resize(nt);

							for(int t = 0 ; t < nt ; ++t)
								dilutions(d1, d2, d3).time[t] = Complex(0.0);
				
						}

			}	
			
			
			struct DilutionComponent_t
			{
				multi1d<DComplex> time;   //Will be length Lt	
			
				DilutionComponent_t& operator*(const DComplex& coeff)
				{
					DilutionComponent_t result;

					result.time = time * coeff;
				
					return result;
				}
			
			};
		
			multi3d<DilutionComponent_t> dilutions;
		
		};
		/*Accessed as follows 
		Sink_qqq_t snk;
		snk.dilutions(d1,d2,d3).time(t);
		*/

		struct Source_qqq_t
		{
			Source_qqq_t(int dil_size) : dilutions.resize(dil_size, dil_size, 
					dil_size) 
			{
				for (int d1 = 0 ; d1 < dil_size ; ++d1)
					for (int d2 = 0 ; d2 < dil_size ; ++d2)
						for (int d3 = 0 ; d3 < dil_size ; ++d3)
						{
							dilutions(d1, d2, d3) = Complex(0.0);
				
						}
			}

			multi3d<DComplex> dilutions;
		
		
		};
		/*Acessed as follows:
		Source_t src;
		src.dilutions(d1,d2,d3);
		*/


		//vector of momenta:   mutli1d< multi1d<int> > momenta( 
		//array of results: multi1d<BaronOpSourceSink_t> bresults( 
		//BaryonOpSourceSink_t(dil_size, nt), momenta.size() ) 
		 

		struct BaryonOpSourceSink_t
		{
			
			BaryonOpSourceSink_t(int dil_size, int nt) : src(dil_size), 
			snk(dil_size, nt) {}
			
			Source_qqq_t src;

			Sink_qqq_t snk;
		
			BaryonOpSourceSink_t& operator*(const DComplex& coeff)
			{

				BaryonOpSourceSink_t result;

				//Do src first
				result.src.dilutions = src.dilutions * coeff;

				//Sink
				result.snk.dilutions = snk.dilutions * coeff;

				return result;
			}

					
			BaryonOpSourceSink_t& operator+=(const BaryonOpSourceSink_t& rhs)
			{

				src.dilutions += rhs.src.dilutions;

				snk.dilutions += rhs.snk.dilutions;
				 
				return this;
			}
		
		};

	
		struct KeyBaryonOpSourceSink_t
		{
			int op_num;
			
			multi1d<int> noises; 
			int t0; 

		};

			bool operator<(const KeyBaryonOpSourceSink_t& a, 
				const KeyBaryonOpSourceSink_t& b);


			void write(BinaryWriter& bin, const KeyBaryonOpSourceSink_t& param);
		
			void read(BinaryReader& bin, KeyBaryonOpSourceSink_t& param);

			void read(XMLReader& xml, const std::string& path, 
				KeyBaryonOpSourceSink_t& param);

		void write(XMLWriter& xml, const std::string& path, 
			const KeyBaryonOpSourceSink_t& param);
	


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
// *   Chroma input format:                                          *
// *                                                                 *
// *      <baryon_operators>                                         *
// *        <elem>                                                   *
// *          <isospin_name> nucleon </isospin_name>                 *
// *          <flavor> uud </flavor>                                 *
// *          <momentum>  0 1 -1  </momentum>                        *
// *          <irrep> Hu </irrep>                                    *
// *          <irrep_row> 3 </irrep_row>                             *
// *          <spatial_type> DDL </spatial_type>                     *
// *          <spatial_id_num> 4 </spatial_id_num>                   *
// *          <disp_length> 3 </disp_length>                         *
// *        </elem>                                                  *
// *        <elem>                                                   *
// *           ....                                                  *
// *        </elem>                                                  *
// *      </baryon_operators>                                        *
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
// *       BaryonOperator::Initialize(coefs_top_directory);          *
// *                                                                 *
// *   first before using this class. The expected way of            *
// *   creating these objects is shown below:                        *
// *                                                                 *
// *       BaryonOperator::Initialize("coefs_top_directory");        *
// *       XMLReader xml_in("input.xml");                            *
// *       multi1d<BaryonOperator> B;                                *
// *       read(xml_in, "/baryon_operators", B);                     *
// *                                                                 *
// *   The above code initializes the class, then reads all of       *
// *   baryon operators inside the <baryon_operators> tags in        *
// *   the input XML file, storing the results in the vector B.      *
// *   Note that each element of B contains the information          *
// *   about its operator and it reads its projection coefficients.  *
// *                                                                 *
// *   The constructor does some checks on the input, then attempts  *
// *   to read the coefficients in terms of elementals for that      *
// *   operator.  If found, the object is assigned and the           *
// *   expansion details and coefficients in terms of elementals     *
// *   is stored.  If not found, an error message is output and      *
// *   the object remains unassigned.  The method "IsAssigned()"     *
// *   can be used to test for the success or failure of the         *
// *   construction.                                                 *
// *                                                                 *
// *   The "read" function above aborts execution if the read of     *
// *   any one baryon operator fails.                                *
// *                                                                 *
// *******************************************************************

class BaryonOperator
{

 public:  

  class Elemental
   {
    int spin1,spin2,spin3;
    int disp1,disp2,disp3;
    int sort_index;

   public:

    Elemental(int s1, int s2, int s3, int d1, int d2, int d3)
     {if ((s1<1)||(s1>4)||(s2<1)||(s2>4)||(s3<1)||(s3>4)
         ||(d1<-3)||(d1>3)||(d2<-3)||(d2>3)||(d3<-3)||(d3>3)){
         QDPIO::cerr << "invalid BaryonOperator::Elemental"<<endl; QDP_abort(1);}
      spin1=s1; spin2=s2; spin3=s3; disp1=d1; disp2=d2; disp3=d3;
      sort_index=((((((((((spin1-1)<<2)+spin2-1)<<2)+spin3-1)<<3)+disp1+3)<<3)
                 +disp2+3)<<3)+disp3+3;}

    Elemental(const Elemental& rhs) 
     : spin1(rhs.spin1), spin2(rhs.spin2), spin3(rhs.spin3),
       disp1(rhs.disp1), disp2(rhs.disp2), disp3(rhs.disp3),
       sort_index(rhs.sort_index) {}

    Elemental& operator=(const Elemental& rhs)
     {spin1=rhs.spin1; spin2=rhs.spin2; spin3=rhs.spin3;
      disp1=rhs.disp1; disp2=rhs.disp2; disp3=rhs.disp3;
      sort_index=rhs.sort_index; return *this;}

    int Spin1() const {return spin1;}
    int Spin2() const {return spin2;}
    int Spin3() const {return spin3;}
    int Disp1() const {return disp1;}
    int Disp2() const {return disp2;}
    int Disp3() const {return disp3;}

    bool operator<(const Elemental& rhs) const
     {return sort_index<rhs.sort_index;}

    string Output() const
     {ostringstream oss;
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
    int disp_length;

    ColorContractSpatialSum(const Elemental& inel, int in_disp_length)
      : el(inel), disp_length(in_disp_length) {}
    ColorContractSpatialSum(const ColorContractSpatialSum& in)
      : el(in.el), disp_length(in.disp_length) {}
    ColorContractSpatialSum& operator=(const ColorContractSpatialSum& in)
      {el=in.el; disp_length=in.disp_length; return *this;}
    bool operator<(const ColorContractSpatialSum& in) const
     {return (disp_length<in.disp_length)
         ||((disp_length==in.disp_length)&&(el<in.el));}
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

  static string coefs_top_directory;
  static bool coefs_flag;

  string isospin_name;
  string flavor;
  Momentum mom;
  string irrep;
  int irrep_row;
  string spatial_type;
  int spatial_id_num;
  int disp_length;

  bool is_assigned;

  list<ElementalTerm> terms;

 public:

  static void Initialize(const string& baryon_coefs_top_directory);

  BaryonOperator() : is_assigned(false) {}

  BaryonOperator(const BaryonOperator& B) 
     : isospin_name(B.isospin_name),flavor(B.flavor),
       mom(B.mom),irrep(B.irrep),irrep_row(B.irrep_row),
       spatial_type(B.spatial_type),
       spatial_id_num(B.spatial_id_num),
       disp_length(B.disp_length),
       is_assigned(B.is_assigned),
       terms(B.terms) {}

  BaryonOperator& operator=(const BaryonOperator& B)
   {isospin_name=B.isospin_name; flavor=B.flavor; 
    mom=B.mom; irrep=B.irrep; irrep_row=B.irrep_row; 
    spatial_type=B.spatial_type; 
    spatial_id_num=B.spatial_id_num; 
    disp_length=B.disp_length;
    is_assigned=B.is_assigned;
    terms=B.terms;
    return *this;}

    // output functions

  string IsospinName() const { return isospin_name; }

  string Flavor() const { return flavor; }

  Momentum MomentumValue() const { return mom; }

  multi1d<int> MomentumMulti1d() const
  {multi1d<int> p(3); p[0]=mom.x; p[1]=mom.y; p[2]=mom.z;
   return p;}

  int XMomentum() const { return mom.x; }

  int YMomentum() const { return mom.y; }

  int ZMomentum() const { return mom.z; }

  string Irrep() const { return irrep; }

  int IrrepRow() const { return irrep_row; }

  string SpatialType() const { return spatial_type; }

  int SpatialIdNumber() const { return spatial_id_num; }

  int DisplacementLength() const { return disp_length; }

  int NumberOfElementals() const { return terms.size(); }

  typedef list<ElementalTerm>::const_iterator ElementalIterator;

  ElementalIterator begin() const { return terms.begin(); }

  ElementalIterator end() const { return terms.end(); }

  bool IsAssigned() const { return is_assigned; }

  string Output() const;

  bool Assign(const string& isospin_tag, const string& flavor_tag,
              const multi1d<int>& momentum, const string& irrep_tag,
              int irrep_row_tag, const string& sptype_tag,
              int spid_tag, int disp_tag);

  void clear()
   {isospin_name.clear(); flavor.clear(); 
    mom.x=0; mom.y=0; mom.z=0; irrep.clear(); irrep_row=0; 
    spatial_type.clear(); spatial_id_num=0; 
    disp_length=0; terms.clear(); is_assigned=false;}

    // returns the flavor string, but replaces all 'u' and
    // 'd' characters by 'l'

  string udFlavor() const
   {string tmp(flavor);
    for (int i=0;i<tmp.length();i++)
       if ((tmp[i]=='u')||(tmp[i]=='d')) tmp[i]='l';
    return tmp;}

  bool operator==(const BaryonOperator& rhs)
  {if ((!is_assigned)||(!rhs.is_assigned)) return false;
   return  (spatial_id_num==rhs.spatial_id_num)
         &&(irrep==rhs.irrep)&&(irrep_row==rhs.irrep_row)
         &&(spatial_type==rhs.spatial_type)
         &&(isospin_name==rhs.isospin_name)&&(flavor==rhs.flavor)
         &&(mom.x==rhs.mom.x)&&(mom.y==mom.y)&&(mom.z==mom.z)
         &&(disp_length==rhs.disp_length);}


 private:

  string trim(const string& str) const;    

  bool get_momentum_ray(int px, int py, int pz, string& ray) const;

};



// **************************************************



void read(XMLReader& xml, const string& path, BaryonOperator& b);
void write(XMLWriter& xml, const string& path, const BaryonOperator& b);





	} // namespace InlineStochGroupBaryonEnv 
}

#endif
