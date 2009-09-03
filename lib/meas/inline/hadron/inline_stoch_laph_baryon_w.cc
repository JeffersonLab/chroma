// $Id: inline_stoch_laph_baryon_w.cc,v 3.4 2009-09-03 15:44:38 colin Exp $
/*! \file
 * \brief Compute the laph-diluted baryon sources and sinks. Write them 
 *  out to db files.  Uses a BaryonSourceSinkHandler.
 *
 * Baryon sources/sinks on laph diluted sources
 */


#include "inline_stoch_laph_baryon_w.h"
#include "chroma.h"

namespace Chroma {
  using namespace LaphEnv;
  namespace InlineStochLaphBaryonEnv {

    //  The crucial create measurement routine. Must be in the *.cc
    //  so that it is local to this file.  Dynamically allocates
    //  and instantiates an object of our class "StochLaphBaryonInlineMeas".

AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
                                        const std::string& path) 
{
 return new StochLaphBaryonInlineMeas(xml_in, path);
}

    //  The name of this inline measurement.   This is the name 
    //  with which the createMeasurement function is associated in the 
    //  Object Factory. You must include this name in the XML input 
    //  to Chroma through
    //     <InlineMeasurements>
    //        <elem>
    //            <Name> STOCH_LAPH_QUARK </Name>
    //             ...
    //        </elem>
    //    </InlineMeasurements>

const std::string name = "STOCH_LAPH_BARYON";

    // Registration boolean hidden in anonymous namespace.
namespace {
   bool registered = false;
}

    // Register all the factories.  This function may be called many
    // times by other measurements, so we only want to register this
    // inline measurement once.  Hence, the use of the "registered"
    // boolean above (which must be hidden in an anonymous namespace).

bool registerAll() 
{
 bool success = true; 
 if (!registered){
    success &= TheInlineMeasurementFactory::Instance().registerObject(
                      name, createMeasurement);
    registered = true;}
 return success;
}

// *********************************************************************
	
     // XML input must have form:
     //
     //   <BaryonSourceSinkInfo> ... </BaryonSourceSinkInfo>

     //   <LaphNoiseList> ... </LaphNoiseList>
     //   <SourceTimeList> ... </SourceTimeList> 

     // Inside the <LaphNoiseList> should be one or more <LaphNoise>
     // tags.
     // Inside the <SourceTimeList> should be either
     //     <Values> 2 5 8 </Values>
     // or  <All> </All>




/*
void StochLaphBaryonInlineMeas::clearSinkComputations()
{
 sinkComputations.clear();
}

void StochLaphBaryonInlineMeas::clearSourceComputations()
{
 sourceComputations.clear();
}
      
void StochLaphBaryonInlineMeas::setSinkComputations(int TimeExtent)
{
 if (!sinkComputations.empty()) sinkComputations.clear();

 if (xml_tag_count(xml_rdr,"SinkComputations")==1){
    XMLReader xmlrd(xml_rdr,"./descendant-or-self::SinkComputations");

    if (xml_tag_count(xmlrd,"NoiseList_TimeList_OneFile")==1){
       XMLReader xmlr(xmlrd,"./descendant-or-self::NoiseList_TimeList_OneFile");
       multi1d<int> source_times;
       if (xml_tag_count(xmlr,"SourceTimeList")==1){
          if (xml_tag_count(xmlr,"SourceTimeList/All")==1){
             source_times.resize(TimeExtent);
             for (int t=0;t<TimeExtent;t++) source_times[t]=t;}
          else
             xmlread(xmlr,"SourceTimeList/Values",source_times,
                     "STOCH_LAPH_QUARK");}
       int num_noises=xml_tag_count(xmlr,"LaphNoiseList/LaphNoiseInfo");
       for (int k=1;k<=num_noises;k++){
          ostringstream path;
          path << "./descendant::LaphNoiseList/LaphNoiseInfo["<<k<<"]";
          XMLReader xml_noise(xmlr,path.str());
          LaphNoiseInfo aNoise(xml_noise);
          for (int t=0;t<source_times.size();t++){
             sinkComputations.push_back(
                  SinkComputation(aNoise,source_times[t]));}}}

    if (xml_tag_count(xmlrd,"ComputationList")==1){
       XMLReader xmlr(xmlrd,"./descendant-or-self::ComputationList");
       int ncomputations=xml_tag_count(xmlr,"Computation");
       for (int k=1;k<=ncomputations;k++){
          ostringstream path;
          path << "./descendant::Computation["<<k<<"]";
          XMLReader xml_comp(xmlr,path.str());
          LaphNoiseInfo aNoise(xml_comp);
          int source_time;
          xmlread(xml_comp,"SourceTime",source_time,"STOCH_LAPH_QUARK");
          sinkComputations.push_back(
                SinkComputation(aNoise,source_time));}}

    }

 QDPIO::cout << endl << "STOCH_LAPH_QUARK sink computations:"<<endl;
 QDPIO::cout << " Number of sink computations = "<<sinkComputations.size()<<endl;
 int count=0;
 for (list<SinkComputation>::const_iterator it=sinkComputations.begin();
      it!=sinkComputations.end();count++,it++){
    QDPIO::cout <<endl<< "SinkComputation "<<count<<":"<<endl;
    QDPIO::cout << it->Noise.output();
    QDPIO::cout << "<SourceTime>"<<it->SourceTime<<"</SourceTime>"<<endl;}

}


void StochLaphBaryonInlineMeas::setSourceComputations()
{
 if (!sourceComputations.empty()) sourceComputations.clear();

 if (xml_tag_count(xml_rdr,"SourceComputations")==1){
    XMLReader xmlrd(xml_rdr,"./descendant-or-self::SourceComputations");

    if (xml_tag_count(xmlrd,"ComputationList")==1){
       XMLReader xmlr(xmlrd,"./descendant-or-self::ComputationList");
       int ncomputations=xml_tag_count(xmlr,"Computation");
       for (int k=1;k<=ncomputations;k++){
          ostringstream path;
          path << "./descendant::Computation["<<k<<"]";
          XMLReader xml_comp(xmlr,path.str());
          LaphNoiseInfo aNoise(xml_comp);
          sourceComputations.push_back(
                SourceComputation(aNoise));}}
    }

 QDPIO::cout << endl << "STOCH_LAPH_QUARK source computations:"<<endl;
 QDPIO::cout << " Number of source computations = "<<sourceComputations.size()<<endl;
 int count=0;
 for (list<SourceComputation>::const_iterator it=sourceComputations.begin();
      it!=sourceComputations.end();count++,it++){
    QDPIO::cout <<endl<< "SourceComputation "<<count<<":"<<endl;
    QDPIO::cout << it->Noise.output();}

}

*/
// *********************************************************************
	
     // Subroutine which does all of the work!!  Input parameters
     // must be as shown (specified by Chroma).  Actual input to
     // this routine is through the private data member
     //     XMLReader xlm_rdr


void StochLaphBaryonInlineMeas::operator()(unsigned long update_no,
                                           XMLWriter& xml_out) 
{

 vector<BaryonOperator> BOps;
 createBaryonOperators(xml_rdr,BOps);

 for (int k=0;k<BOps.size();k++){
    QDPIO::cout << "Baryon operator "<<k<<":"<<endl;
    QDPIO::cout << BOps[k].fullOutput()<<endl;}

 XMLBufferWriter xmlout;
 BOps[0].output(xmlout);
 QDPIO::cout << "First baryon operator from XMLWriter:"<<endl;
 QDPIO::cout << xmlout.str()<<endl;

    // create the handler and set up the info from the
    // XML <BaryonSourceSinkInfo> tag
/*
 BaryonSourceSinkHandler Q(xml_rdr);

 QDPIO::cout << endl <<"Info initialized in BaryonSourceSinkHandler"<<endl;
 {XMLBufferWriter xmlout;
 Q.getFileMap(xmlout);
 cout << xmlout.str()<<endl;}

    // set or compute the Laph eigenvectors (smears the gauge field as needed)

 Q.setLaphEigenvectors();

    // read the list of computations (noises, time sources, file indices)
    // from xml_rdr and store in the "Computations" data member

 setSourceComputations();
 setSinkComputations(Q.getGaugeConfigurationInfo().getTimeExtent());

 START_CODE();
 StopWatch outer;
 outer.start();

 int count=0;
 for (list<SourceComputation>::const_iterator it=sourceComputations.begin();
      it!=sourceComputations.end();count++,it++){
    QDPIO::cout <<endl<< "Now starting source computation "<<count<<":"<<endl;
    Q.computeSource(it->Noise);}

 count=0;
 for (list<SinkComputation>::const_iterator it=sinkComputations.begin();
      it!=sinkComputations.end();count++,it++){
    QDPIO::cout <<endl<< "Now starting sink computation "<<count<<":"<<endl;
    Q.computeSink(it->Noise,it->SourceTime);}

 outer.stop();
 QDPIO::cout << name << ": total time = " << outer.getTimeInSeconds() 
             << " secs" << endl;
 QDPIO::cout << name << ": ran successfully" << endl;

 END_CODE(); */
} 

// ******************************************************************
  }
} // namespace Chroma












/*

#include "handle.h"
#include "meas/inline/hadron/inline_stoch_laph_baryon_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/smear/displacement.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"
#include <sstream> 
#include <fstream>
#include <complex>

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
	/*!
	 * \ingroup hadron
	 *
	 * @{
	 * /
	namespace InlineStochLapHBaryonEnv 
	{ 
		//! Number of quarks to be used in this construction
		const int N_quarks = 3;


		void read(XMLReader& xml, const string& path, BaryonOperator& b)
		{
			string isospin_tag,flavor_tag,irrep_tag,sptype_tag;
			int irrep_row_tag,spid_tag,disp_tag;
			multi1d<int> p;

			XMLReader paramtop(xml, path);
			try{
				read(paramtop, "isospin_name", isospin_tag);
				read(paramtop, "flavor", flavor_tag);
				read(paramtop, "momentum", p);
				read(paramtop, "irrep", irrep_tag);
				read(paramtop, "irrep_row", irrep_row_tag);
				read(paramtop, "spatial_type", sptype_tag);
				read(paramtop, "spatial_id_num", spid_tag);
				read(paramtop, "disp_length", disp_tag);
			}
			catch(const string& error){
				QDPIO::cerr << "problem reading a baryon operator: "<<error <<endl;
				QDP_abort(1);
			}

			if (!b.Assign(isospin_tag,flavor_tag,p,irrep_tag,
						irrep_row_tag,sptype_tag,spid_tag,disp_tag)){
				QDPIO::cerr << "invalid baryon operator requested"<<endl;
				QDP_abort(1);
			}
		}

		void write(XMLWriter& xml, const string& path, const BaryonOperator& b)
		{
			push(xml, path);
			write(xml, "isospin_name", b.IsospinName());
			write(xml, "flavor", b.Flavor());
			write(xml, "momentum", b.MomentumMulti1d());
			write(xml, "irrep", b.Irrep());
			write(xml, "irrep_row", b.IrrepRow());
			write(xml, "spatial_type", b.SpatialType());
			write(xml, "spatial_id_num", b.SpatialIdNumber());
			write(xml, "disp_length", b.DisplacementLength());
			pop(xml);
		}

		void BaryonOperator::Initialize(const string& baryon_coefs_top_directory)
		{
			coefs_top_directory=baryon_coefs_top_directory;
			if (*coefs_top_directory.rbegin()!='/') coefs_top_directory+="/";
			string checkfile=coefs_top_directory
				+"nucleon_uud/mom_ray_000/G1g_1/SS_0";
			ifstream in(checkfile.c_str());
			if (in){
				coefs_flag=true;
				in.close();}
			else{
				QDPIO::cerr << "Baryon coefficient top directory not found"<<endl;
				in.close();
				QDP_abort(1);}
		}


		bool BaryonOperator::Assign(const string& isospin_tag, const string& flavor_tag,
				const multi1d<int>& p_momentum, const string& irrep_tag,
				int irrep_row_tag, const string& sptype_tag,
				int spid_tag, int disp_tag)
		{
			clear();
			if (!coefs_flag){
				QDPIO::cerr << "Must call BaryonOperator::Initialize(...) before creating baryon objects"<<endl;
				QDP_abort(1);}

				isospin_name=trim(isospin_tag);
				flavor=trim(flavor_tag);
				irrep=trim(irrep_tag);
				spatial_type=trim(sptype_tag);
				if (p_momentum.size()!=3){
					clear();
					QDPIO::cerr << "Bad baryon momentum"<<endl;
					return false;}
					mom.x=p_momentum[0];
					mom.y=p_momentum[1];
					mom.z=p_momentum[2];

					irrep_row=irrep_row_tag;
					disp_length=disp_tag;
					spatial_id_num=spid_tag;
					if ((disp_length<0)||(spatial_id_num<0)){
						QDPIO::cerr << "Bad baryon operator input xml data"<<endl;
						clear();
						return false;}

						// check if a single-site operator, then make the displacement length zero
						if (spatial_type=="SS") disp_length=0;

						// first, determine the momentum "ray"

						string mom_ray;
						if (!get_momentum_ray(mom.x,mom.y,mom.z,mom_ray)){
							QDPIO::cerr << "Baryon momentum is not part of allowed ray"<<endl;
							clear();
							return false;}

							// now look for the baryon operator in the XML coefficient file

							ostringstream oss;
							oss << coefs_top_directory << isospin_name << "_" << flavor
								<< "/" << mom_ray << "/" << irrep << "_" << irrep_row
								<< "/" << spatial_type << "_" << spatial_id_num;
							string coefs_file=oss.str();

							// cout << coefs_file << endl;

							ifstream in(coefs_file.c_str());
							if (!in){
								QDPIO::cerr << "Invalid baryon requested"<<endl;
								clear();
								in.close();
								return false;}

								int nterms;
								Double re,im;
								int s1,s2,s3,d1,d2,d3;
								in >> nterms;
								for (int k=0;k<nterms;k++){
									in >> s1 >> s2 >> s3 >> d1 >> d2 >> d3 >> re >> im;
									terms.push_back(ElementalTerm(Elemental(s1,s2,s3,d1,d2,d3),
												cmplx(re,im)));}

									if (in.fail()){
										clear();}
									else
										is_assigned=true;

									in.close();
									return true;
		}

		// remove leading and trailing blanks

		string BaryonOperator::trim(const string& str) const
		{
			int start=str.find_first_not_of(" ");
			int len=str.find_last_not_of(" ")-start+1;
			return str.substr(start,len);
		}

		//   Allowed momentum rays:
		//     000  +00  0+0  00+  ++0  +-0  +0+  +0-  0++  0+-
		//     +++  ++-  +-+  +--

		bool BaryonOperator::get_momentum_ray(int px, int py, int pz, string& ray) const
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

		string BaryonOperator::Output() const
		{
			ostringstream oss;
			if (is_assigned){
				oss << "<baryon>"<<endl;
				oss << "  <isospin_name> " << isospin_name << " </isospin_name>"<<endl;
				oss << "  <flavor> "<<flavor<<" </flavor>"<<endl;
				oss << "  <momentum> " <<mom.x<<" "<<mom.y
					<<" "<<mom.z<<" </momentum>"<<endl;
				oss << "  <irrep> "<<irrep<<" </irrep>"<<endl;
				oss << "  <irrep_row> "<<irrep_row<<" </irrep_row>"<<endl;
				oss << "  <spatial_type> "<<spatial_type<<" </spatial_type>"<<endl;
				oss << "  <spatial_id_num> "<<spatial_id_num
					<<" </spatial_id_num>"<<endl;
				oss << "  <disp_length> "<<disp_length
					<<" </disp_length>"<<endl;
				oss << "  <projection>"<<endl;
				for (BaryonOperator::ElementalIterator it=terms.begin();it!=terms.end();it++){
					oss << "     <term> "<<it->el.Output()
						<< " coef=cmplx"<<it->coef<<" </term>"<<endl;
				}
				oss << "  </projection>"<<endl;
				oss << "</baryon>"<<endl;
			}
			return oss.str();
		}



		string BaryonOperator::coefs_top_directory;
		bool BaryonOperator::coefs_flag=false;


		//
		// The spin basis matrix to goto Dirac
		//
		SpinMatrix rotate_mat(adj(DiracToDRMat()));

		// Reader for input parameters
		void read(XMLReader& xml, const string& path, InlineStochLapHBaryonEnv::Params::Param_t& param)
		{
			XMLReader paramtop(xml, path);

			int version;
			read(paramtop, "version", version);

			multi1d< multi1d<int> > temp;
			switch (version) 
			{
				case 1:

					break;

				default :

					QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
					QDP_abort(1);
			}

			read(paramtop, "BaryonOperators", param.bops);

			param.link_smearing         = readXMLGroup(paramtop, "LinkSmearing", "LinkSmearingType");

			read(paramtop, "NumNoises", param.nnoises);
			read(paramtop, "UDMass", param.ud_mass_label);
			read(paramtop, "SMass", param.s_mass_label);
			read(paramtop, "CoeffDir", param.coeff_top_dir);
			read(paramtop, "all_ords", param.all_ords);

		}


		// Writer for input parameters
		void write(XMLWriter& xml, const string& path, const InlineStochLapHBaryonEnv::Params::Param_t& param)
		{
			push(xml, path);

			int version = 1;

			write(xml, "version", version);
			write(xml, "BaryonOperators", param.bops);
			write(xml, "NumNoises", param.nnoises);
			write(xml, "UDMass", param.ud_mass_label);
			write(xml, "SMass", param.s_mass_label);
			write(xml, "CoeffDir", param.coeff_top_dir);
			write(xml, "AllOrds", param.all_ords);
			xml << param.link_smearing.xml;

			pop(xml);
		}


		//! Read named objects 
		void read(XMLReader& xml, const string& path, InlineStochLapHBaryonEnv::Params::NamedObject_t& input)
		{
			XMLReader inputtop(xml, path);

			read(inputtop, "gauge_id", input.gauge_id);
			read(inputtop, "BaryonOutfile", input.baryon_file);
			read(inputtop, "QuarkFiles", input.quark_files);
		}

		//! Write named objects
		void write(XMLWriter& xml, const string& path, const InlineStochGroupBaryonEnv::Params::NamedObject_t& input)
		{
			push(xml, path);

			write(xml, "gauge_id", input.gauge_id);
			write(xml, "BaryonOutfile", input.baryon_file);
			write(xml, "QuarkFiles", input.quark_files);

			pop(xml);
		}
	}


	namespace InlineStochLapHBaryonEnv 
	{ 
		// Anonymous namespace for registration
		namespace
		{
			AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					const std::string& path) 
			{
				return new InlineMeas(Params(xml_in, path));
			}

			//! Local registration flag
			bool registered = false;
		}

		const std::string name = "STOCH_LAPH_BARYON";

		//! Register all the factories
		bool registerAll() 
		{
			bool success = true; 
			if (! registered)
			{
				success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
				registered = true;
			}
			return success;
		}


		//----------------------------------------------------------------------------
		// Param stuff
		Params::Params()
		{ 
			frequency = 0; 
			param.mom2_max = 0;
		}

		Params::Params(XMLReader& xml_in, const std::string& path) 
		{
			try 
			{
				XMLReader paramtop(xml_in, path);

				if (paramtop.count("Frequency") == 1)
					read(paramtop, "Frequency", frequency);
				else
					frequency = 1;

				// Read program parameters
				read(paramtop, "Param", param);

				// Read in the output propagator/source configuration info
				read(paramtop, "NamedObject", named_obj);

				// Possible alternate XML file pattern
				if (paramtop.count("xml_file") != 0) 
				{
					read(paramtop, "xml_file", xml_file);
				}
			}
			catch(const std::string& e) 
			{
				QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << endl;
				QDP_abort(1);
			}
		}


		void
			Params::writeXML(XMLWriter& xml_out, const std::string& path) 
			{
				push(xml_out, path);

				// Parameters for source construction
				write(xml_out, "Param", param);

				// Write out the output propagator/source configuration info
				write(xml_out, "NamedObject", named_obj);

				pop(xml_out);
			}



		//--------------------------------------------------------------------------
		//Support for the diquarks

		void makeDiquark( multi1d<LatticeComplex> & diquark, const multi1d<LatticeComplex> & q0,
				const multi1d<LatticeComplex> & q1, const Subset & subset )
		{


			//The signs for the diquark are taken from
			//the colorContract function in qdp_primcolorvec.h
			diquark[0][subset] =  q0[0]*q1[1] - q0[1]*q1[0];
			diquark[1][subset] =  q0[1]*q1[2] - q0[2]*q1[1];
			diquark[2][subset] =  q0[2]*q1[0] - q0[0]*q1[2];


		}


		void makeColorSinglet (LatticeComplex & singlet, const multi1d<LatticeComplex> & diquark, 
				const multi1d<LatticeComplex> & q2, const Subset & subset)
		{

			singlet[subset] = diquark[0] * q2[2];
			singlet[subset] += diquark[1] * q2[0];  
			singlet[subset] += diquark[2] * q2[1];  
		}


		//--------------------------------------------------------------------------
		// Function call
		void 
			InlineMeas::operator()(unsigned long update_no,
					XMLWriter& xml_out) 
			{
				// If xml file not empty, then use alternate
				if (params.xml_file != "")
				{
					string xml_file = makeXMLFileName(params.xml_file, update_no);

					push(xml_out, "stoch_laph_baryon");
					write(xml_out, "update_no", update_no);
					write(xml_out, "xml_file", xml_file);
					pop(xml_out);

					XMLFileWriter xml(xml_file);
					func(update_no, xml);
				}
				else
				{
					func(update_no, xml_out);
				}
			}


		// Function call
		void 
			InlineMeas::func(unsigned long update_no,
					XMLWriter& xml_out) 
			{
				START_CODE();

				StopWatch snoop;
				snoop.reset();
				snoop.start();


				StopWatch swiss;

				// Test and grab a reference to the gauge field
				XMLBufferWriter gauge_xml;
				try
				{
					TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
					TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
				}
				catch( std::bad_cast ) 
				{
					QDPIO::cerr << InlineStochGroupBaryonEnv::name << ": caught dynamic cast error" 
						<< endl;
					QDP_abort(1);
				}
				catch (const string& e) 
				{
					QDPIO::cerr << InlineStochGroupBaryonEnv::name << ": map call failed: " << e 
						<< endl;
					QDP_abort(1);
				}
				const multi1d<LatticeColorMatrix>& u = 
					TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

				push(xml_out, "StochLapHBaryon");
				write(xml_out, "update_no", update_no);

				QDPIO::cout << InlineStochGroupBaryonEnv::name << 
					": Stochastic LapH-Diluted Baryon Operators" << endl;

				proginfo(xml_out);    // Print out basic program info

				// Write out the input
				params.writeXML(xml_out, "Input");

				// Write out the config info
				write(xml_out, "Config_info", gauge_xml);

				push(xml_out, "Output_version");
				write(xml_out, "out_version", 1);
				pop(xml_out);

				//First calculate some gauge invariant observables just for info.
				//This is really cheap.
				MesPlq(xml_out, "Observables", u);

				//
				// Initialize the slow Fourier transform phases
				//

				// Sanity check - if this doesn't work we have serious problems
				if (phases.numSubsets() != QDP::Layout::lattSize()[decay_dir])
				{
					QDPIO::cerr << name << ": number of time slices not equal to that in the decay direction: " 
						<< QDP::Layout::lattSize()[decay_dir]
						<< endl;
					QDP_abort(1);
				}

				//
				// Smear the gauge field if needed
				//
				multi1d<LatticeColorMatrix> u_smr = u;

				try
				{
					std::istringstream  xml_l(params.param.link_smearing.xml);
					XMLReader  linktop(xml_l);
					QDPIO::cout << "Link smearing type = " << params.param.link_smearing.id << endl;


					Handle< LinkSmearing >
						linkSmearing(TheLinkSmearingFactory::Instance().createObject(params.param.link_smearing.id,
									linktop, params.param.link_smearing.path));

					(*linkSmearing)(u_smr);
				}
				catch(const std::string& e) 
				{
					QDPIO::cerr << name << ": Caught Exception link smearing: " << e << endl;
					QDP_abort(1);
				}

				MesPlq(xml_out, "Smeared_Observables", u_smr);	

				//Initialize the quark database manager	
				QuarkDBManager  qdb_man(params.named_obj.quark_files, 
						params.param.nnoises);
				
				QDPIO::cout << "Num Orderings = " << qdp_man.getNumOrds() << endl;

				//
				// Baryon operators
				//
				BaryonOperator::Initialize(params.param.coeff_top_dir);


				for(int t0 = 0; t0 < participating_timeslices.size() ; ++t0)
				{

					//Loop over orderings

					//Loop over elementals

					//Initialize the Slow Fourier Transform Object for this elemental 
					int decay_dir = diluted_quarks[0]->getDecayDir();

					SftMom phases(params.param.moms, decay_dir);

					//Loop over moms 



					//Make the source op
					
					//Loop over source orderings 
					QDPIO::cout << "Source Ordering = " << ord << endl;

					//Loop over d1
					//Loop over d2
					{

						keySmearedDispColorVector[0].dil = i;
						keySmearedDispColorVector[1].dil = j;

						//Form the di-quark to save on recalculating 
						multi1d<LatticeComplex> diquark(Nc);

						const multi1d<LatticeComplex> &q0 = smrd_disp_srcs.getDispSource(n0, 
								keySmearedDispColorVector[0]); 

						const multi1d<LatticeComplex> &q1 = smrd_disp_srcs.getDispSource(n1, 
								keySmearedDispColorVector[1]);


						watch.reset();
						watch.start();
						//For the source, restrict this operation to a subset
						makeDiquark( diquark, q0 , q1, phases.getSet()[ participating_timeslices[t0] ] ); 
						watch.stop();

						/*QDPIO::cout<< " Made diquark : time = " << 
							watch.getTimeInSeconds() << "secs" << endl;
							* /		

						for(int k = 0 ; k < diluted_quarks[n2]->getDilSize(t0) ; ++k)	
						{

							keySmearedDispColorVector[2].dil = k;

							// Contract over color indices with antisym tensor.
							// NOTE: the creation operator only lives on a time slice, so restrict
							// the operation to that time slice

							LatticeComplex c_oper;

							const multi1d<LatticeComplex> &q2 = smrd_disp_srcs.getDispSource(n2, 
									keySmearedDispColorVector[2]);

							watch.reset();
							watch.start();

							makeColorSinglet( c_oper, diquark, q2, phases.getSet()[ 
									participating_timeslices[t0] ] );

							watch.stop();

							/*QDPIO::cout<< "Made Color singlet : time =  " <<  
								watch.getTimeInSeconds() << "secs" << endl;
								* /	

							// Slow fourier-transform
							// We can restrict what the FT routine requires to a subset.
							watch.reset();
							watch.start();



							multi2d<DComplex> c_sum;
							int num_mom;

							c_sum = phases.sft(c_oper, participating_timeslices[t0]);
							num_mom = phases.numMom();

							watch.stop();


						} // end for k
					} // end for j
				} // end for i

				//End loop over orderings

		swiss.stop();

		QDPIO::cout << "Source operator construction: operator= " << l 
			<< "  time= "
			<< swiss.getTimeInSeconds() 
			<< " secs" << endl;



		int ord = 0;
		//for(int ord = 0 ; ord < num_orderings ; ++ord)
		{
			QDPIO::cout << "Ordering = " << ord << endl;

			annih_oper.time_slices[0].orderings[ord].perm = perms[ord];

			const int n0 = perms[ord][0];
			const int n1 = perms[ord][1];
			const int n2 = perms[ord][2];

			// The operator must hold all the dilutions
			// We know that all time slices match. However, not all time slices of the
			// lattice maybe used

			// Creation operator
			BaryonOperator_t::TimeSlices_t::Orderings_t& aop = annih_oper.time_slices[0].orderings[ord];

			aop.dilutions.resize(diluted_quarks[n0]->getDilSize(t0), diluted_quarks[n1]->getDilSize(t0),
					diluted_quarks[n2]->getDilSize(t0) );

			for (int n = 0 ; n < N_quarks ; ++n)
			{
				keySmearedDispColorVector[n].t0 = t0;
			}

			for(int i = 0 ; i <  diluted_quarks[n0]->getDilSize(t0) ; ++i)
			{
				for(int j = 0 ; j < diluted_quarks[n1]->getDilSize(t0) ; ++j)	      
				{

					keySmearedDispColorVector[0].dil = i;
					keySmearedDispColorVector[1].dil = j;

					//Form the di-quark to save on recalculating 
					multi1d<LatticeComplex> diquark(Nc);

					const multi1d<LatticeComplex> &q0 = smrd_disp_snks.getDispSolution(n0, 
							keySmearedDispColorVector[0]); 

					const multi1d<LatticeComplex> &q1 = smrd_disp_snks.getDispSolution(n1, 
							keySmearedDispColorVector[1]);


					//QDPIO::cout<<"q0[0] testval= "<< peekSite(q0[0], orig)
					//	<< endl; 

					//QDPIO::cout<<"q1[0] testval= "<< peekSite(q1[0], orig)
					//	<< endl; 


					watch.reset();
					watch.start();

					makeDiquark( diquark, q0 , q1, all ); 

					watch.stop();
					/*QDPIO::cout << "Made diquark: time = " << 
						watch.getTimeInSeconds() << "secs " << endl;
						* /

					for(int k = 0 ; k < diluted_quarks[n2]->getDilSize(t0) ; ++k)	
					{

						keySmearedDispColorVector[2].dil = k;

						// Contract over color indices with antisym tensor.
						// There is a potential optimization here - the colorcontract of
						// the first two quarks could be pulled outside the innermost dilution
						// loop.
						// NOTE: the creation operator only lives on a time slice, so restrict
						// the operation to that time slice

						LatticeComplex a_oper;

						const multi1d<LatticeComplex> &q2 = smrd_disp_snks.getDispSolution(n2, 
								keySmearedDispColorVector[2]);

						//QDPIO::cout<<"q2[0] testval= "<< peekSite(q2[0], orig)
						//<< endl;

						watch.reset();
						watch.start();

						makeColorSinglet( a_oper, diquark, q2, all);

						watch.stop();

						/*
							 QDPIO::cout <<	"Made Color Singlet: time = " <<
							 watch.getTimeInSeconds() << "secs" << endl;
							 */
						/*QDPIO::cout << "testval = " << peekSite(a_oper, orig) 
							<< endl;
							* /

						watch.reset();
						watch.start();

						// Slow fourier-transform
						multi2d<DComplex> a_sum;
						int num_mom;

						a_sum = phases.sft(
								a_oper);
						num_mom = phases.numMom();

						watch.stop();
						/*
							 QDPIO::cout << "Spatial Sums completed: time " << 
							 watch.getTimeInSeconds() << "secs" << endl;
							 * /		
						// Unpack into separate momentum and correlator
						aop.dilutions(i,j,k).mom_projs.resize(num_mom);

						for(int mom_num = 0 ; mom_num < num_mom ; ++mom_num) 
						{
							aop.dilutions(i,j,k).mom_projs[mom_num].mom = params.param.moms[mom_num];

							aop.dilutions(i,j,k).mom_projs[mom_num].op = a_sum[mom_num];

						}

					} // end for k
				} // end for j
			} // end for i
		}//end ord 
		swiss.stop();


		QDPIO::cout << "Sink operator construction: operator= " << l 
			<< "  time= "
			<< swiss.getTimeInSeconds() 
			<< " secs" << endl;

		QDPIO::cout << "Sink op testval( t0 = " << 
			participating_timeslices[t0] << ") = " << 
			annih_oper.time_slices[0].orderings[0].dilutions(0,0,0).mom_projs[0].op[0] 
			<< endl;

		//Hard code the elemental op name for now 
		std::stringstream cnvrt;
		cnvrt <<  annih_oper.id  << "_t" << participating_timeslices[t0] << "_snk.lime";

		std::string filename;

		filename = cnvrt.str(); 

		// Write the meta-data and the binary for this operator
		swiss.reset();
		swiss.start();
		{
			XMLBufferWriter     src_record_xml, file_xml;
			BinaryBufferWriter  src_record_bin;

			push(file_xml, "SinkBaryonOperator");
			write(file_xml, "Params", params.param);
			write(file_xml, "Config_info", gauge_xml);
			write(file_xml, "Op_Info",qqq_oplist.ops[l]);
			push(file_xml, "QuarkSources");

			push(file_xml, "Quark_l");
			push(file_xml, "TimeSlice");
			push(file_xml, "Dilutions");
			for (int dil = 0; dil < diluted_quarks[0]->getDilSize(t0) ; ++dil)
			{
				write( file_xml, "elem", 
						diluted_quarks[0]->getSourceHeader(t0, dil) );
			}
			pop(file_xml); //dilutions 
			pop(file_xml); //TimeSlice
			pop(file_xml); //Quark_l

			push(file_xml, "Quark_m");
			push(file_xml, "TimeSlice");
			push(file_xml, "Dilutions");
			for (int dil = 0; dil < diluted_quarks[1]->getDilSize(t0) ; ++dil)
			{
				write( file_xml, "elem", 
						diluted_quarks[1]->getSourceHeader(t0, dil) );
			}
			pop(file_xml); //dilutions 
			pop(file_xml); //TimeSlice
			pop(file_xml); //Quark_m

			push(file_xml, "Quark_r");
			push(file_xml, "TimeSlice");
			push(file_xml, "Dilutions");
			for (int dil = 0; dil < diluted_quarks[2]->getDilSize(t0) ; ++dil)
			{
				write( file_xml, "elem", 
						diluted_quarks[2]->getSourceHeader(t0, dil) );
			}
			pop(file_xml); //dilutions 
			pop(file_xml); //TimeSlice
			pop(file_xml); //Quark_r

			pop(file_xml);//QuarkSources
			push(file_xml, "QuarkSinks");

			push(file_xml, "Quark_l");
			write(file_xml, "PropHeader", diluted_quarks[0]->getPropHeader(0,0) );
			pop(file_xml);

			push(file_xml, "Quark_m");
			write(file_xml, "PropHeader", diluted_quarks[1]->getPropHeader(0,0) );
			pop(file_xml);

			push(file_xml, "Quark_r");
			write(file_xml, "PropHeader", diluted_quarks[2]->getPropHeader(0,0) );
			pop(file_xml);

			pop(file_xml);//QuarkSinks 
			pop(file_xml);//SinkBaryonOperator

			QDPFileWriter qdp_file(file_xml, filename,     // are there one or two files???
					QDPIO_SINGLEFILE, QDPIO_SERIAL, QDPIO_OPEN);


			write(src_record_xml, "BaryonAnnihilationOperator", annih_oper);
			write(src_record_bin, annih_oper);

			write(qdp_file, src_record_xml, src_record_bin);

		}
		swiss.stop();

		QDPIO::cout << "Sink Operator writing: operator = " << l
			<< "  time= " << swiss.getTimeInSeconds() << " secs" << endl;

	} // end for l (operator )

} //End Make annihilation operator

} //end t0 
// Close the namelist output file XMLDAT
pop(xml_out);     // StochBaryon

snoop.stop();
	QDPIO::cout << InlineStochGroupBaryonEnv::name << ": total time = " 
<< snoop.getTimeInSeconds() 
	<< " secs" << endl;

	QDPIO::cout << InlineStochGroupBaryonEnv::name << ": ran successfully" << endl;

	END_CODE();
	} // func

} // namespace InlineStochGroupBaryonEnv

/*! @} * /  // end of group hadron

} // namespace Chroma
*/
