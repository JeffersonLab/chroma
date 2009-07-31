// $Id: inline_stoch_laph_quark_w.cc,v 3.9 2009-07-31 19:15:03 colin Exp $
/*! \file
 * \brief Compute the matrix element of  LatticeColorVector*M^-1*LatticeColorVector
 *
 * Propagator calculation on a colorvector
 */

#include "fermact.h"
#include "meas/inline/hadron/inline_stoch_laph_quark_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ferm/subset_vectors.h"
#include "util/ferm/map_obj.h"
#include "util/ferm/key_val_db.h"
#include "util/ferm/key_laph_dil_prop.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/smear/quark_smearing_factory.h"
#include "meas/smear/quark_smearing_aggregate.h"
#include "meas/sources/source_smearing_aggregate.h"
#include "meas/sources/source_smearing_factory.h"
#include "meas/sinks/sink_smearing_aggregate.h"
#include "meas/sinks/sink_smearing_factory.h"
#include "meas/laph/laph.h"

#include "meas/inline/io/named_objmap.h"

#include "meas/inline/hadron/inline_make_source_ferm_w.h"
#include "meas/sources/source_const_factory.h"
#include "meas/sources/source_const_aggregate.h"
#include "util/ferm/diractodr.h"

namespace Chroma {
  using namespace LaphEnv;
  namespace InlineStochLaphQuarkEnv {
    namespace {

AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, const std::string& path) 
 {return new StochLaphQuarkInlineMeas(xml_in, path);}

   //! Local registration flag
bool registered = false;
}

const std::string name = "STOCH_LAPH_QUARK";

  //! Register all the factories
bool registerAll() 
{
 bool success = true; 
 if (!registered){
    success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
    registered = true;}
 return success;
}

// *********************************************************************
	
     // Subroutine which does all of the work!!

void StochLaphQuarkInlineMeas::operator()(unsigned long update_no,
                                          XMLWriter& xml_out) 
{
				
 QuarkSourceSinkHandler Q(xml_rdr);
  
 QDPIO::cout << "Constructed Gauge Handler" << endl;
 
   FieldSmearingHandler FSH(xml_rdr);

   QDPIO::cout << "Constructed FieldSmearing Handler" << endl;
				
                              
				FSH.computeSmearedGaugeField();	 
				
				QDPIO::cout << "Computed Smeared Gauge Field" << endl;
				MesPlq(xml_out, "Observables", 
						FSH.getSmearedGaugeField() );
					 
                                //FSH.computeLaphEigenvectors();
                                
   // const multi1d<LatticeColorVector>& ev=FSH.getLaphEigenvectors();
    
    


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
