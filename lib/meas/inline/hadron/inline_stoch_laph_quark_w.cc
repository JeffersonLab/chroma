// $Id: inline_stoch_laph_quark_w.cc,v 3.10 2009-08-21 14:54:51 colin Exp $
/*! \file
 * \brief Compute the laph-diluted quark sources and sinks. Write them 
 *  out to db files.  Uses a QuarkSourceSinkHandler.
 *
 * Propagator calculation on laph diluted sources
 */


#include "inline_stoch_laph_quark_w.h"
#include "chroma.h"

namespace Chroma {
  using namespace LaphEnv;
  namespace InlineStochLaphQuarkEnv {

    //  The crucial create measurement routine. Must be in the *.cc
    //  so that it is local to this file.  Dynamically allocates
    //  and instantiates an object of our class "StochLaphQuarkInlineMeas".

AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
                                        const std::string& path) 
{
 return new StochLaphQuarkInlineMeas(xml_in, path);
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

const std::string name = "STOCH_LAPH_QUARK";

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
     //   <QuarkSourceSinkInfo> ... </QuarkSourceSinkInfo>

     //   <LaphNoiseList> ... </LaphNoiseList>
     //   <SourceTimeList> ... </SourceTimeList> 

     // Inside the <LaphNoiseList> should be one or more <LaphNoise>
     // tags.
     // Inside the <SourceTimeList> should be either
     //     <Values> 2 5 8 </Values>
     // or  <All> </All>





void StochLaphQuarkInlineMeas::clearSinkComputations()
{
 sinkComputations.clear();
}

void StochLaphQuarkInlineMeas::clearSourceComputations()
{
 sourceComputations.clear();
}
      
void StochLaphQuarkInlineMeas::setSinkComputations(int TimeExtent)
{
 if (!sinkComputations.empty()) sinkComputations.clear();

 if (xml_tag_count(xml_rdr,"SinkComputations")==1){
    XMLReader xmlrd(xml_rdr,"./descendant-or-self::SinkComputations");

    int default_file_index=-1;
    if (xml_tag_count(xmlrd,"DefaultFileIndex")==1)
       xmlread(xmlrd,"DefaultFileIndex",default_file_index,"STOCH_LAPH_QUARK");

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
             sinkComputations.push_back(SinkComputation(aNoise,source_times[t],
                                                        default_file_index));}}}

    if (xml_tag_count(xmlrd,"ComputationList")==1){
       XMLReader xmlr(xmlrd,"./descendant-or-self::ComputationList");
       int ncomputations=xml_tag_count(xmlr,"Computation");
       for (int k=1;k<=ncomputations;k++){
          ostringstream path;
          path << "./descendant::Computation["<<k<<"]";
          XMLReader xml_comp(xmlr,path.str());
          LaphNoiseInfo aNoise(xml_comp);
          int file_index=default_file_index;
          if (xml_tag_count(xml_comp,"FileIndex")==1)
             xmlread(xml_comp,"FileIndex",file_index,"STOCH_LAPH_QUARK");
          int source_time;
          xmlread(xml_comp,"SourceTime",source_time,"STOCH_LAPH_QUARK");
          sinkComputations.push_back(SinkComputation(aNoise,source_time,
                                                     file_index));}}

    }

 QDPIO::cout << endl << "STOCH_LAPH_QUARK sink computations:"<<endl;
 QDPIO::cout << " Number of sink computations = "<<sinkComputations.size()<<endl;
 int count=0;
 for (list<SinkComputation>::const_iterator it=sinkComputations.begin();
      it!=sinkComputations.end();count++,it++){
    QDPIO::cout <<endl<< "SinkComputation "<<count<<":"<<endl;
    QDPIO::cout << it->Noise.output();
    QDPIO::cout << "<SourceTime>"<<it->SourceTime<<"</SourceTime>"<<endl;
    QDPIO::cout << "<FileIndex>"<<it->FileIndex<<"</FileIndex>"<<endl;}

}


void StochLaphQuarkInlineMeas::setSourceComputations()
{
 if (!sourceComputations.empty()) sourceComputations.clear();

 if (xml_tag_count(xml_rdr,"SourceComputations")==1){
    XMLReader xmlrd(xml_rdr,"./descendant-or-self::SourceComputations");

    int default_file_index=-1;
    if (xml_tag_count(xmlrd,"DefaultFileIndex")==1)
       xmlread(xmlrd,"DefaultFileIndex",default_file_index,"STOCH_LAPH_QUARK");

    if (xml_tag_count(xmlrd,"ComputationList")==1){
       XMLReader xmlr(xmlrd,"./descendant-or-self::ComputationList");
       int ncomputations=xml_tag_count(xmlr,"Computation");
       for (int k=1;k<=ncomputations;k++){
          ostringstream path;
          path << "./descendant::Computation["<<k<<"]";
          XMLReader xml_comp(xmlr,path.str());
          LaphNoiseInfo aNoise(xml_comp);
          int file_index=default_file_index;
          if (xml_tag_count(xml_comp,"FileIndex")==1)
             xmlread(xml_comp,"FileIndex",file_index,"STOCH_LAPH_QUARK");
          sourceComputations.push_back(SourceComputation(aNoise,
                                                         file_index));}}

    }

 QDPIO::cout << endl << "STOCH_LAPH_QUARK source computations:"<<endl;
 QDPIO::cout << " Number of source computations = "<<sourceComputations.size()<<endl;
 int count=0;
 for (list<SourceComputation>::const_iterator it=sourceComputations.begin();
      it!=sourceComputations.end();count++,it++){
    QDPIO::cout <<endl<< "SourceComputation "<<count<<":"<<endl;
    QDPIO::cout << it->Noise.output();
    QDPIO::cout << "<FileIndex>"<<it->FileIndex<<"</FileIndex>"<<endl;}

}


// *********************************************************************
	
     // Subroutine which does all of the work!!  Input parameters
     // must be as shown (specified by Chroma).  Actual input to
     // this routine is through the private data member
     //     XMLReader xlm_rdr


void StochLaphQuarkInlineMeas::operator()(unsigned long update_no,
                                          XMLWriter& xml_out) 
{

    // create the handler and set up the info from the
    // XML <QuarkSourceSinkInfo> tag

 QuarkSourceSinkHandler Q(xml_rdr);  

    // read the list of computations (noises, time sources, file indices)
    // from xml_rdr and store in the "Computations" data member

 setSourceComputations();
 setSinkComputations(Q.getGaugeConfigurationInfo().getTimeExtent());

 int count=0;
 for (list<SourceComputation>::const_iterator it=sourceComputations.begin();
      it!=sourceComputations.end();count++,it++){
    QDPIO::cout << "Now starting source computation "<<count<<":"<<endl;
    Q.computeSource(it->Noise,it->FileIndex);}

 count=0;
 for (list<SinkComputation>::const_iterator it=sinkComputations.begin();
      it!=sinkComputations.end();count++,it++){
    QDPIO::cout << "Now starting sink computation "<<count<<":"<<endl;
    Q.computeSink(it->Noise,it->SourceTime,it->FileIndex);}

/*
 GaugeConfigurationInfo G(xml_rdr);

 GaugeConfigurationInfo G2(G);

 cout << "G2 time extent = "<<G2.getTimeExtent()<<endl;
 cout << "G2 time dir = "<<G2.getTimeDir()<<endl;
 cout << "G2 num dir = "<<G2.getNumberOfDirections()<<endl;
 cout << "G2 num traj = "<<G2.getTrajNum()<<endl;
 cout << "G2 gauge_id = "<<G2.getGaugeId()<<endl;

 cout << "G2 output = "<<endl;
 cout << G2.output()<<endl;
 cout << "G2 done"<<endl;

 cout << "are equal? "; assertEqual(G2,G,"inlineMeas");
 cout << "okay"<<endl;
 if (G2==G) cout << "equal true"<<endl;
 else cout << "equal false"<< endl;

// string header("<GaugeConfigurationInfo><onetag>here</onetag></GaugeConfigurationInfo>");
 //GaugeConfigurationInfo G3(header);


 string header=G.getGaugeConfigHeader();
 GaugeConfigurationInfo G4(header);

 try{
    G.checkEqual(G4);}
 catch(const string& err){
    cout << "Did not match"<<endl;}
 cout << "match okay!!"<<endl;

cout << G.getFullRecordXML()<<endl;

// assertEqual(G,G3,"inlineMeas");


 QuarkInfo SQ(xml_rdr,G);

 cout << " mass is "<<SQ.getMass()<<endl;

 string qhead = SQ.getQuarkHeader();

 cout << "qhead = "<<qhead<<endl;

 QuarkInfo SQ2(qhead);

 cout << "SQ2 = "<<SQ2.getQuarkHeader()<<endl;
 cout << "SQ2 "<<SQ2.getMassName()<<" = "<<SQ2.getMass()<<endl;
 cout << "SQ2 action id = "<<SQ2.getActionId()<<endl;

 assertEqual(SQ,SQ2,"inlineMeas");
 cout << "they match"<<endl;
*/
/*
 cout << "  ****************************************"<<endl<<endl;

 FieldSmearingInfo fs(xml_rdr);
 cout << "fs header is "<<fs.getHeader()<<endl;

 string fshead=fs.getHeader();

 FieldSmearingInfo fs2(fshead);
 cout << "fs2 header is "<<fs2.getHeader()<<endl;

 try{
    fs.checkEqual(fs2);
    cout << "fs and fs2 match"<<endl;}
 catch(const string& err2){
    cout << "fs and fs2 did not match"<<endl;}


 LaphNoiseInfo noise1(xml_rdr);
 cout << "noise1 header is "<<noise1.getHeader()<<endl;

 string nhead = noise1.getHeader();
 LaphNoiseInfo noise2(nhead);
 cout << "noise2 header is "<<noise2.getHeader()<<endl;

 try{
    noise1.checkEqual(noise2);
    cout <<"noise1 and noise2 are equal"<<endl;}
 catch(const string& err3){
    cout <<"noise1 and noise2 do not match"<<endl;}



 DilutionSchemeInfo dilscheme1(xml_rdr);
 cout << "dilscheme1 header is "<<dilscheme1.getHeader()<<endl;

 string dhead = dilscheme1.getHeader();
 DilutionSchemeInfo dilscheme2(dhead);
 cout << "dilscheme2 header is "<<dilscheme2.getHeader()<<endl;

 try{
    dilscheme1.checkEqual(dilscheme2);
    cout <<"dilscheme1 and dilscheme2 are equal"<<endl;}
 catch(const string& err3){
    cout <<"dilscheme1 and dilscheme2 do not match"<<endl;}
*/

/*
 GaugeConfigurationHandler uHandler(xml_rdr);
 cout << "uHandler header:" << uHandler.getGaugeConfigHeader()<<endl;
 uHandler.setData();
 const multi1d<LatticeColorMatrix>& U=uHandler.getData();

 {XMLBufferWriter out_xml;
 MesPlq(out_xml, "UnsmearedWLoops", U);
 QDPIO::cout << out_xml.str() << endl;}

 uHandler.setInfo(xml_rdr);
 cout << uHandler.getInfo().getTimeExtent()<<endl;
*/
/*
 stringstream oss;
 oss << " <laph_dilution_scheme>" << endl;                           
 oss << "    <time_dilution>" << endl;                               
 oss << "       <dilution_type> full </dilution_type>" << endl;      
 oss << "    </time_dilution>" << endl;                              
 oss << "    <spin_dilution>" << endl;                               
 oss << "       <dilution_type> none </dilution_type>" << endl;      
 oss << "    </spin_dilution>" << endl;                              
 oss << "    <eigvec_dilution>" << endl;                             
 oss << "       <dilution_type> block </dilution_type>" << endl;     
 oss << "       <number_projectors> 4 </number_projectors>" << endl; 
 oss << "    </eigvec_dilution>" << endl;                            
 oss << " </laph_dilution_scheme>" << endl;                          
 XMLReader xml2(oss);
 XMLReader xmlB(xml2,"/");

 DilutionSchemeInfo dilB(xmlB);

 cout << dilB.output() << endl;
*/

				

// cout << Q.getHeader()<<endl;
 
/*
   FieldSmearingHandler FSH(xml_rdr);

   QDPIO::cout << "Constructed FieldSmearing Handler" << endl;
				
                              
 FSH.computeSmearedGaugeField();	 
				
QDPIO::cout << "Computed Smeared Gauge Field" << endl;


{XMLBufferWriter out_xml;
 MesPlq(out_xml, "SmearedWLoops", FSH.getSmearedGaugeField());
 QDPIO::cout << out_xml.str() << endl;}
					 
FSH.computeLaphEigenvectors();
                                
const multi1d<LatticeColorVector>& ev=FSH.getLaphEigenvectors();
    
  */  


START_CODE();

				/*
					 StopWatch snoop;
					 snoop.reset();
					 snoop.start();

					 SpinMatrix sp_one = 1.0;

					 SpinMatrix gamma_4 = Gamma(8) * sp_one;

					 SpinMatrix rotate_mat(adj(DiracToDRMat()));


					 push(xml_out, "LaphDilutedProps");
					 write(xml_out, "update_no", update_no);

					 QDPIO::cout << name << ": LapH diluted propagator calculation" << endl;

					 proginfo(xml_out);    // Print out basic program info

				// Write out the input
				write(xml_out, "Input", params);

				// Write out the config header
				write(xml_out, "Config_info", gauge_xml);

				push(xml_out, "Output_version");
				write(xml_out, "out_version", 1);
				pop(xml_out);

				// Calculate some gauge invariant observables just for info.
				MesPlq(xml_out, "Observables", u);


				//Sanity checks on noise sources
				int n_noises = params.param.ran_seeds.size();

				for (int r1 = 0 ; r1 < n_noises ; ++r1)
				for (int r2 = 0 ; r2 < r1 ; ++r2)
				{
				if (toBool(params.param.ran_seeds[r1] == params.param.ran_seeds[r2]) )
				{
				QDPIO::cerr << "ERROR: Seed " << r1 << " = " << r2 << endl;
				QDP_abort(1);
				}
				}

				//Sanity checks on the eigenvector dilutions

				//First check that the dilution components sum to unity
				int n_ev_dil = params.param.eig_vec_dils.size();

				std::vector<int> total_dils;

				for (int d = 0 ;  d < n_ev_dil ; ++d)
				{
				int curr_size = params.param.eig_vec_dils[d].size();

				for (int e = 0 ; e < curr_size ; ++e)
				total_dils.push_back(params.param.eig_vec_dils[d][e]);
				}

				sort( total_dils.begin(), total_dils.end() );

				for (int v = 0 ; v < total_dils.size() ; ++v)
				{
				if (total_dils[v] != v)
				{
				QDPIO::cerr << 
				"ERROR: the eigenvector dilution projectors do not sum to unity"
				<< endl;

				QDP_abort(1);
			}
			}


			//
			// DB storage
			//
			BinaryStoreDB< SerialDBKey<KeyLaphDilutedProp_t>, SerialDBData<LatticeFermion> > qdp_db_props;

			// Open the file, and write the meta-data and the binary for this operator
			{
				XMLBufferWriter file_xml;

				push(file_xml, "DBMetaData");
				write(file_xml, "id", string("dilutedLaphProp"));
				write(file_xml, "lattSize", QDP::Layout::lattSize());
				write(file_xml, "decay_dir", params.param.decay_dir);
				proginfo(file_xml);    // Print out basic program info
				write(file_xml, "Params", params.param);
				write(file_xml, "Config_info", gauge_xml);
				pop(file_xml);

				std::string file_str(file_xml.str());
				qdp_db_props.setMaxUserInfoLen(file_str.size());

				qdp_db_props.open(params.named_obj.prop_file, O_RDWR | O_CREAT, 0664);

				qdp_db_props.insertUserdata(file_str);
			}

			// Total number of iterations
			int ncg_had = 0;
			//
			// Try the factories
			//
			try
			{
				StopWatch swatch;
				swatch.reset();
				QDPIO::cout << "Try the various factories" << endl;

				// Typedefs to save typing
				typedef LatticeFermion               T;
				typedef multi1d<LatticeColorMatrix>  P;
				typedef multi1d<LatticeColorMatrix>  Q;

				//
				// Initialize fermion action
				//
				std::istringstream  xml_s(params.param.prop.fermact.xml);
				XMLReader  fermacttop(xml_s);
				QDPIO::cout << "FermAct = " << params.param.prop.fermact.id << endl;

				// Generic Wilson-Type stuff
				Handle< FermionAction<T,P,Q> >
					S_f(TheFermionActionFactory::Instance().createObject(params.param.prop.fermact.id,
								fermacttop,
								params.param.prop.fermact.path));

				Handle< FermState<T,P,Q> > state(S_f->createState(u));

				Handle< SystemSolver<LatticeFermion> > PP = S_f->qprop(state,
						params.param.prop.invParam);

				QDPIO::cout << "Suitable factory found: compute all the quark props" << endl;
				swatch.start();



				//Initialize quark smearing operator for the sinks
				Handle< QuarkSmearing<LatticeFermion> > quarkSmearing;

				try
				{
					QDPIO::cout << "Create quark smearing object" << endl;

					XMLBufferWriter smr_buf;
					push(smr_buf, "QuarkSmearing");
					write(smr_buf, "wvf_kind", "VECTOR_SMEAR");
					write(smr_buf, "sigma", 0.0);
					write(smr_buf, "subset_vecs_id", params.named_obj.eigvec_id);
					write(smr_buf, "no_smear_dir", params.param.decay_dir);
					pop(smr_buf);
					// Create the quark smearing object
					XMLReader  smeartop(smr_buf);

					quarkSmearing =
						TheFermSmearingFactory::Instance().createObject("VECTOR_SMEAR",
								smeartop, "/QuarkSmearing");

				}
				catch(const std::string& e) 
				{
					QDPIO::cerr << ": Caught Exception creating quark smearing object: " << e << endl;
					QDP_abort(1);
				}
				catch(...)
				{
					QDPIO::cerr << ": Caught generic exception creating smearing object" << endl;
					QDP_abort(1);
				}


				//
				// Loop over the source color and spin, creating the source
				// and calling the relevant propagator routines.
				//
				const int decay_dir           = params.param.decay_dir;
				const multi1d<int>& t_sources = params.param.t_sources;

				SftMom phases(0, true, Nd-1);

				// Binary output
				push(xml_out, "LaphDilutedProps");

				// Loop over each operator 
				QDPIO::cout << "Spin_dil = " << params.param.spin_dil << endl;


				for (int r = 0 ; r < n_noises ; ++r)
				{
					const Seed& curr_seed = params.param.ran_seeds[r];

					int nspins;

					if (params.param.spin_dil)
						nspins = Ns;
					else
						nspins = 1;

					for(int spin_source = 0 ; spin_source < nspins ; ++spin_source)
					{
						QDPIO::cout << "spin_dil = " << spin_source << endl; 

						for (int vec_dil = 0 ; vec_dil < params.param.eig_vec_dils.size() ; 
								++vec_dil)
						{
							const multi1d<int>& curr_vecs = params.param.eig_vec_dils[vec_dil];


							LatticeFermion source_out = zero;

							for(int tt=0; tt < t_sources.size(); ++tt)
							{
								int t_source = t_sources[tt];
								QDPIO::cout << "t_source = " << t_source << endl; 

								multi1d<int> t_arr(1);
								t_arr[0] = t_source;
								//Make diluted noisy source xml	
								XMLBufferWriter src_buf;
								push(src_buf, "SourceParams");
								write(src_buf, "SourceType", "RAND_DILUTE_EIGVEC_ZN_SOURCE");
								write(src_buf, "version", 1);
								write(src_buf, "ran_seed", curr_seed);
								write(src_buf, "N", 4);
								write(src_buf, "j_decay", params.param.decay_dir);
								write(src_buf, "t_sources", t_arr);
								write(src_buf, "eigen_vec_id", params.named_obj.eigvec_id);
								write(src_buf, "eigen_vectors", curr_vecs);

								multi1d<int> spin_msk;

								if (params.param.spin_dil)
								{
									spin_msk.resize(1);
									spin_msk[0] = spin_source;
								}
								else
								{
									spin_msk.resize(Ns);
									for (int sp = 0 ; sp < Ns ; ++sp)
										spin_msk[sp] = sp;
								}

								write(src_buf, "spin_mask", spin_msk);
								pop(src_buf);

								XMLReader src_rdr(src_buf);

								//Create noisy source
								Handle< QuarkSourceConstruction<LatticeFermion> >
									sourceConstruction(TheFermSourceConstructionFactory::Instance().createObject(
												"RAND_DILUTE_EIGVEC_ZN_SOURCE",
												src_rdr,
												"/SourceParams"));

								//This source lives only on a single timeslice
								LatticeFermion curr_source = (*sourceConstruction)(u);


								//Multiply this source by gamma_4
								curr_source = gamma_4 * curr_source;	

								//Perform the inversion
								LatticeFermion curr_solution = zero;

								SystemSolverResults_t res = (*PP)(curr_solution, curr_source);

								//Rotate the source and solution from the DR to the DP basis
								curr_source = rotate_mat * curr_source;
								curr_solution = rotate_mat * curr_solution;

								//Keep all the t0's in a single LatticeFermion
								source_out[phases.getSet()[t_source]] += curr_source;

								//Smear Solution with the vector smearing operator
								(*quarkSmearing)(curr_solution,u);

								//Insert the solution into the map
								KeyLaphDilutedProp_t snk_key; 

								snk_key.src_or_snk = "SNK";
								snk_key.t0 = t_source;
								snk_key.spin_dil = spin_source;
								snk_key.evec_dil = vec_dil;
								snk_key.noise_src = r;


								SerialDBKey<KeyLaphDilutedProp_t> snk_db_key;
								snk_db_key.key() = snk_key;

								SerialDBData<LatticeFermion> snk_db_data;
								snk_db_data.data() = curr_solution;

								qdp_db_props.insert(snk_db_key, snk_db_data);

							} // for t_source

							//Insert source into map 
							KeyLaphDilutedProp_t src_key;
							src_key.src_or_snk = "SRC";
							src_key.t0 = 0; 
							src_key.spin_dil = spin_source;
							src_key.evec_dil = vec_dil;
							src_key.noise_src = r;

							SerialDBKey<KeyLaphDilutedProp_t> src_db_key;
							src_db_key.key() = src_key;

							SerialDBData<LatticeFermion> src_db_data;
							src_db_data.data() = source_out;

							qdp_db_props.insert(src_db_key, src_db_data);

						}//for vec dil

					} // for spin_source

				}// for r
				pop(xml_out);

				swatch.stop();
				QDPIO::cout << "Propagators computed: time= " 
					<< swatch.getTimeInSeconds() 
					<< " secs" << endl;
			}
			catch (const std::string& e) 
			{
				QDPIO::cout << name << ": caught exception around qprop: " << e << endl;
				QDP_abort(1);
			}

			push(xml_out,"Relaxation_Iterations");
			write(xml_out, "ncg_had", ncg_had);
			pop(xml_out);

			pop(xml_out);  // prop_matelem_colorvec

			snoop.stop();
			QDPIO::cout << name << ": total time = "
				<< snoop.getTimeInSeconds() 
				<< " secs" << endl;

			*/
				QDPIO::cout << name << ": ran successfully" << endl;

			END_CODE();
			} 

	}

} // namespace Chroma
