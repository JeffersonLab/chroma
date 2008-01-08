// $Id: inline_stoch_group_baryon_w.cc,v 1.8 2008-01-08 18:59:34 jbulava Exp $
/*! \file
 * \brief Inline measurement of stochastic group baryon operator
 *
 */

#include "handle.h"
#include "meas/inline/hadron/inline_stoch_group_baryon_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/smear/quark_smearing_factory.h"
#include "meas/smear/quark_smearing_aggregate.h"
#include "meas/sources/source_smearing_aggregate.h"
#include "meas/sources/source_smearing_factory.h"
#include "meas/sinks/sink_smearing_aggregate.h"
#include "meas/sinks/sink_smearing_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/smear/displacement.h"
#include "util/ferm/diractodr.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"
#include <sstream> 

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  /*!
   * \ingroup hadron
   *
   * @{
   */
  namespace InlineStochGroupBaryonEnv 
  { 
    // Reader for input parameters
    void read(XMLReader& xml, const string& path, InlineStochGroupBaryonEnv::Params::Param_t& param)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
	/**************************************************************************/
	break;

      default :
	/**************************************************************************/

	QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
	QDP_abort(1);
      }

      read(paramtop, "mom2_max", param.mom2_max);
      read(paramtop, "displacement_length", param.displacement_length);

      param.source_quark_smearing = readXMLGroup(paramtop, "CreationOperatorSmearing", "wvf_kind");
      param.sink_quark_smearing   = readXMLGroup(paramtop, "AnnihilationOperatorSmearing", "wvf_kind");
      param.link_smearing         = readXMLGroup(paramtop, "LinkSmearing", "LinkSmearingType");
    }


    // Writer for input parameters
    void write(XMLWriter& xml, const string& path, const InlineStochGroupBaryonEnv::Params::Param_t& param)
    {
      push(xml, path);

      int version = 1;

      write(xml, "version", version);
      write(xml, "mom2_max", param.mom2_max);
      write(xml, "displacement_length", param.displacement_length);
      xml << param.source_quark_smearing.xml;
      xml << param.sink_quark_smearing.xml;
      xml << param.link_smearing.xml;

      pop(xml);
    }


    //Reader for 3-quark operator file
    void read(XMLReader& xml, const string& path, InlineStochGroupBaryonEnv::Params::NamedObject_t::ThreeQuarkOpsFile_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "ops_file", input.ops_file);
      read(inputtop, "id", input.id);
    }


    // Writer for 3-quark operator file
    void write(XMLWriter& xml, const string& path, const InlineStochGroupBaryonEnv::Params::NamedObject_t::ThreeQuarkOpsFile_t& input)
    {
      push(xml, path);
      write(xml, "ops_file", input.ops_file);
      write(xml, "id", input.id);
      pop(xml);
    }


    //Read spin dilution files   
    void read(XMLReader& xml, const string& path, InlineStochGroupBaryonEnv::Params::NamedObject_t::QuarkFiles_t::TimeDilutions_t::SpinDilutions_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "dilution_files", input.dilution_files);
    }


    //Write spin dilution files 
    void write(XMLWriter& xml, const string& path, const InlineStochGroupBaryonEnv::Params::NamedObject_t::QuarkFiles_t::TimeDilutions_t::SpinDilutions_t& input)
    {
      push(xml, path);
      write(xml, "dilution_files", input.dilution_files);
      pop(xml);
    }

    //Read time dilution files 
    void read(XMLReader& xml, const string& path, InlineStochGroupBaryonEnv::Params::NamedObject_t::QuarkFiles_t::TimeDilutions_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "spin_files", input.spin_files);
		}


    //Write time dilution files 
    void write(XMLWriter& xml, const string& path, const InlineStochGroupBaryonEnv::Params::NamedObject_t::QuarkFiles_t::TimeDilutions_t& input)
    {
      push(xml, path);
      write(xml, "spin_files", input.spin_files);
      pop(xml);
    }



    //Read Quark dilution files 
    void read(XMLReader& xml, const string& path, InlineStochGroupBaryonEnv::Params::NamedObject_t::QuarkFiles_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "time_files", input.time_files);
    }


    //Write Quark dilution files 
    void write(XMLWriter& xml, const string& path, const InlineStochGroupBaryonEnv::Params::NamedObject_t::QuarkFiles_t& input)
    {
      push(xml, path);
      write(xml, "time_files", input.time_files);
      pop(xml);
    }


    //! Read named objects 
    void read(XMLReader& xml, const string& path, InlineStochGroupBaryonEnv::Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "operators_file", input.operators_file);
      read(inputtop, "Quarks", input.quarks);
      read(inputtop, "Quark_ids", input.quark_ids);
    }

    //! Write named objects
    void write(XMLWriter& xml, const string& path, const InlineStochGroupBaryonEnv::Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "operators_file", input.operators_file);
      write(xml, "Quarks", input.quarks);

      pop(xml);
    }
  }


  namespace InlineStochGroupBaryonEnv 
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

    const std::string name = "STOCH_GROUP_BARYON";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= QuarkSourceSmearingEnv::registerAll();
	success &= QuarkSinkSmearingEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }



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



    //--------------------------------------------------------------

    //! Structure holding a source and its solutions
    struct QuarkSourceSolutions_t
    {
      //! Structure holding solutions
      struct TimeDilutions_t
      {
				struct SpinDilutions_t
				{
					struct Dilutions_t
					{
	  
						
						LatticeFermion     source;
	  				LatticeFermion     soln;
	  				PropSourceConst_t  source_header;
	  				ChromaProp_t       prop_header;
					};

	  			multi1d<int>       s0;   //Spins included in this dilution
					
					multi1d<Dilutions_t>  dilutions; //dilutions per timeslice per spin
				
				};
				
	  		multi1d<int>       t0;  //Times included in this dilution
				
				multi1d<SpinDilutions_t> spin_dilutions;
			
				int findSpin (const int spin) const //return the spin dilution element to which 'spin' belongs
				{
					int ret_s0 = 0;

					bool stop = false; 

					for (int s0 = 0 ; s0 < spin_dilutions.size() ; ++s0)
					{
						for (int spin0 = 0 ; spin0 < spin_dilutions[s0].s0.size() ; ++spin0)
						{
							if (spin == spin_dilutions[s0].s0[spin0])
							{
								ret_s0 = s0;
								stop = true; 
								break;
							}
						}
						
						if (stop)
						{
							break; 
						}
							
					}

					return ret_s0; 
				
				}
			
			};

      
			int   decay_dir;
      Seed  seed;
     	multi1d<TimeDilutions_t>  time_dilutions;
    
			int findTime( const int time) const   //return the time dilution element to which a particular time 
			{																						//belongs
				
				int ret_t0 = 0;
				bool stop = false; 
				
				for ( int t0 = 0 ; t0 < time_dilutions.size() ;  ++t0 )
				{
					for ( int time0 = 0 ; time0 < time_dilutions[t0].t0.size() ; ++time0 )
					{
						if ( time == time_dilutions[t0].t0[ time0 ] ) 
						{	
							ret_t0 = t0;	
							stop = true;
							break;
						}
				
					}
					
					if (stop)
					{
						break;
					}
				
				}
			
				return ret_t0;
			}			
		
		};


    //! 3-quark operator structure
	struct ThreeQuarkOps_t
	{
		
		struct ThreeQuarkOp_t
		{
	  
			struct QuarkInfo_t
	  	{
	    	int  displacement;    /*!< Orig plus/minus 1-based directional displacements */
	    	int  spin;            /*!< 0-based spin index */
	    	int  disp_dir;        /*!< 0-based direction */
	    	int  disp_len;        /*!< 0-based length */
	  	};

	  	multi1d<QuarkInfo_t>  quark;    /*!< Displacement and spin for each quark */
	  	std::string name;                 /*!< Name of the 3-quark operator */
		};
	
		multi1d<ThreeQuarkOp_t> ops; /*!< 3-quark ops within a file */
	};


    //! Structure holding structures

    //! The key for smeared and displaced color vectors
    struct KeySmearedDispColorVector_t
    {
      int  displacement;    /*!< Orig plus/minus 1-based directional displacements */
      int  spin;            /*!< 0-based spin index */
    };


    //! Support for the keys of smeared and displaced color vectors
    bool operator<(const KeySmearedDispColorVector_t& a, const KeySmearedDispColorVector_t& b)
    {
      multi1d<int> lga(2);
      lga[0] = a.displacement;
      lga[1] = a.spin;

      multi1d<int> lgb(2);
      lgb[0] = b.displacement;
      lgb[1] = b.spin;

      return (lga < lgb);
    }


    //! The smeared and displaced color vectors
    struct SmearedDispColorVector_t
    {
      struct Quarks_t
      {
	struct TimeDilutions_t
	{
		struct SpinDilutions_t
		{
			
	  	struct Dilutions_t
	  	{
	    	LatticeColorVector  source;
	    	LatticeColorVector  soln;
	  	};

		
			multi1d<int> s0;	//Array of participating spins for this dilution element 
	 		multi1d<Dilutions_t> dilutions; //Additional dilutions   
		
		};

		multi1d<int> t0;  //Array of participating timeslices for this dilution element
		SpinDilutions_t spin_dilutions;    //Only one spin dilution element needs to be calculated
	};	
	
	multi1d<TimeDilutions_t> time_dilutions;
      };

      multi1d<Quarks_t>  quarks;
    };

//    map< KeySmearedDispColorVector_t, SmearedDispColorVector_t > quarks;


    //! Baryon operator
    struct BaryonOperator_t
    {
      //! Quark orderings within a baryon operator
      struct Orderings_t
      {
	//! Baryon operator time slices
	struct TimeDilutions_t
	{
	  struct SpinDilutions_t               //A single set of 3 spin dilution indicies 
		{
			//! Baryon operator dilutions
	  	struct Dilutions_t
	  	{
	    	//! Momentum projected correlator
	    	struct Mom_t
	    	{
	      	multi1d<int>       mom;       /*!< D-1 momentum of this operator */
	      	multi1d<DComplex>  op;        /*!< Momentum projected operator */
	    	};

	  		multi1d<Mom_t> mom_projs;       /*!< Holds momentum projections of the operator */
	  
			};
	  	
			multi1d<int> s0_l;  // participating spins for left quark
			multi1d<int> s0_m;  // participating spins for middle quark
			multi1d<int> s0_r;  // participating spins for right quark
			
			multi3d<Dilutions_t> dilutions;   /*!< Additional dilution indices  */
		};

	  multi1d<int>  t0_l;          /*!< Left quark participating time_slices for this dilution element */
	  multi1d<int>  t0_m;          /*!< Middle quark participating time_slices for this dilution element */
	  multi1d<int>  t0_r;          /*!< Right quark participating time_slices for this dilution element */
	
		SpinDilutions_t spin;        /*!< Only one combination of the 3 spin dilution elements needs to be calculated */
		
	};
	  
	multi1d<int> perm;                  /*!< This particular permutation of quark orderings */
	
	multi1d<TimeDilutions_t> time_dilutions;  /*!< combination of time dilution elements for each quark */
      };

      multi1d< multi1d<int> > perms;   /*!< Permutations of quark enumeration */

      GroupXML_t    smearing;          /*!< String holding quark smearing xml */

      Seed          seed_l;            /*!< Id of left quark */
      Seed          seed_m;            /*!< Id of middle quark */
      Seed          seed_r;            /*!< Id of right quark */

      int           operator_num;      /*!< Operator number within file */
      std::string   id;                /*!< Tag/ID used in analysis codes */

      int           mom2_max;          /*!< |\vec{p}|^2 */
      int           decay_dir;         /*!< Direction of decay */
      multi1d<Orderings_t> orderings;  /*!< Array is over quark orderings */
    };


    //! BaryonOperator header writer
    void write(XMLWriter& xml, const string& path, const BaryonOperator_t& param)
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);
      write(xml, "id", param.id);
      write(xml, "operator_num", param.operator_num);
      write(xml, "mom2_max", param.mom2_max);
      write(xml, "decay_dir", param.decay_dir);
      write(xml, "seed_l", param.seed_l);
      write(xml, "seed_m", param.seed_m);
      write(xml, "seed_r", param.seed_r);
      write(xml, "perms", param.perms);
      xml <<  param.smearing.xml;

      pop(xml);
    }


    //! BaryonOperator binary writer
    void write(BinaryWriter& bin, const BaryonOperator_t::Orderings_t::TimeDilutions_t::SpinDilutions_t::Dilutions_t::Mom_t& param)
    {
      write(bin, param.mom);
      write(bin, param.op);
    }

    //! BaryonOperator binary writer
    void write(BinaryWriter& bin, const BaryonOperator_t::Orderings_t::TimeDilutions_t::SpinDilutions_t::Dilutions_t& param)
    {
      write(bin, param.mom_projs);
    }

    //! BaryonOperator binary writer
    void write(BinaryWriter& bin, const BaryonOperator_t::Orderings_t::TimeDilutions_t::SpinDilutions_t& param)
    {
      write(bin, param.dilutions);
    	
			write(bin, param.s0_l);
    	write(bin, param.s0_m);
    	write(bin, param.s0_r);
		}
    
		//! BaryonOperator binary writer
    void write(BinaryWriter& bin, const BaryonOperator_t::Orderings_t::TimeDilutions_t& param)
    {
      write(bin, param.spin);
    	
			write(bin, param.t0_l);
    	write(bin, param.t0_m);
    	write(bin, param.t0_r);
		}

    //! BaryonOperator binary writer
    void write(BinaryWriter& bin, const BaryonOperator_t::Orderings_t& param)
    {
      write(bin, param.perm);
      write(bin, param.time_dilutions);
    }

    //! BaryonOperator binary writer
    void write(BinaryWriter& bin, const BaryonOperator_t& param)
    {
      write(bin, param.seed_l);
      write(bin, param.seed_m);
      write(bin, param.seed_r);
      write(bin, param.mom2_max);
      write(bin, param.decay_dir);
      write(bin, param.perms);
      write(bin, param.orderings);
    }




    //! Read 3-quark operators file
    void readOps(ThreeQuarkOps_t& oplist, 
		    const std::string& ops_file,
		    int displacement_length)
    {
      START_CODE();

      TextFileReader reader(ops_file);

      int num_ops;
      reader >> num_ops;
      oplist.ops.resize(num_ops);
			
			//Loop over ops within a file
      for(int n=0; n < oplist.ops.size(); ++n)
      {
				std::string name; 
				reader >> name;
				oplist.ops[n].name = name;  
	
	  ThreeQuarkOps_t::ThreeQuarkOp_t& qqq = oplist.ops[n];
	  qqq.quark.resize(3);

	  // Make spin index 0 based
	  {
	    multi1d<int> spin(3);
	    reader >> spin[0] >> spin[1] >> spin[2];

	    for(int i=0; i < spin.size(); ++i)
	      qqq.quark[i].spin = spin[i] - 1;
	  }

	  // Convert displacements to  disp_dir, disp_len
	  {
	    multi1d<int> displacement(3);
	    reader >> displacement[0] >> displacement[1] >> displacement[2];

	    for(int i=0; i < displacement.size(); ++i)
	    {
	      qqq.quark[i].displacement = displacement[i];

	      if (displacement[i] == 0)
	      {
		qqq.quark[i].disp_dir = 0;
		qqq.quark[i].disp_len = 0;
	      }
	      else if (displacement[i] > 0)
	      {
		qqq.quark[i].disp_dir = qqq.quark[i].displacement - 1;
		qqq.quark[i].disp_len = displacement_length;
	      }
	      else
	      {
		qqq.quark[i].disp_dir = -qqq.quark[i].displacement - 1;
		qqq.quark[i].disp_len = -displacement_length;
	      }
	    } 
	  }

      } //n
	
      reader.close();

      END_CODE();
    } //void


    //! Construct array of maps of displacements
    void displaceQuarks(map<KeySmearedDispColorVector_t, SmearedDispColorVector_t>& disp_quarks,
			const multi1d<LatticeColorMatrix>& u_smr,
			const ThreeQuarkOps_t&  oplist,
			const multi1d<QuarkSourceSolutions_t> quarks)
    {
      START_CODE();

			StopWatch snoop;
			
      cout << __func__ << ": entering" << endl;

	  // Loop over the ops
	  for(int l=0; l < oplist.ops.size(); ++l)
	  {
	    const ThreeQuarkOps_t::ThreeQuarkOp_t& qqq_op = oplist.ops[l];
	    
	    //Loop over quark terms within the op
			for(int i = 0; i < qqq_op.quark.size(); ++i)
	    {
	      const ThreeQuarkOps_t::ThreeQuarkOp_t::QuarkInfo_t& q_info = qqq_op.quark[i];

	      // Check if this displacement and spin are in the map,
	      // if not, then *ALL* quarks and dilutions must be shifted
	      KeySmearedDispColorVector_t key;
	      key.displacement = q_info.displacement;
	      key.spin         = q_info.spin;

	      // If no entry, then create a displaced version of the quark
	      if (disp_quarks.find(key) == disp_quarks.end())
	      {
//	      cout << __func__ 
//		   << ": n=" << n
//		   << " l=" << l
//		   << " i=" << i 
//		   << " disp=" << term.quark[i].displacement
//		   << " len=" << term.quark[i].disp_len
//		   << " dir=" << term.quark[i].disp_dir
//		   << endl;

		// Insert an empty entry and then modify it. This saves on
		// copying the data around
		{
		  SmearedDispColorVector_t disp_empty;
		  disp_quarks.insert(std::make_pair(key, disp_empty));

		  // Sanity check - the entry better be there
		  if (disp_quarks.find(key) == disp_quarks.end())
		  {
		    QDPIO::cerr << __func__ 
				<< ": internal error - could not insert empty key in map"
				<< endl;
		    QDP_abort(1);
		  }		      
		}

		// Modify the previous empty entry
		SmearedDispColorVector_t& disp_q = disp_quarks.find(key)->second;
		disp_q.quarks.resize(quarks.size());
		
		snoop.reset();
		snoop.start();
	
		//Loop over all quarks 
		for(int n=0; n < disp_q.quarks.size(); ++n)
		{
		  disp_q.quarks[n].time_dilutions.resize(quarks[n].time_dilutions.size());

		  for(int t0 = 0; t0 < disp_q.quarks[n].time_dilutions.size(); ++t0)
		  {
			
					  //Denotes which dilution element contains the desired spin 
						const int spin = quarks[n].time_dilutions[t0].findSpin( key.spin );
						
						disp_q.quarks[n].time_dilutions[t0].spin_dilutions.dilutions.resize(quarks[n].time_dilutions[t0].spin_dilutions[ spin ].dilutions.size());
						
						for(int m=0; m < disp_q.quarks[n].time_dilutions[t0].spin_dilutions.dilutions.size(); ++m)
		    		{
		      	
							// Short-hand
		      		const QuarkSourceSolutions_t::TimeDilutions_t::SpinDilutions_t::Dilutions_t& qq = 
			quarks[n].time_dilutions[t0].spin_dilutions[ spin ].dilutions[ m ];

		      		SmearedDispColorVector_t::Quarks_t::TimeDilutions_t::SpinDilutions_t::Dilutions_t& dil = 
			disp_q.quarks[n].time_dilutions[t0].spin_dilutions.dilutions[ m ];
					
					

		      		// Pull out the appropriate spin component, then displace it
							
		      		dil.source = peekSpin(qq.source, q_info.spin);
		      		dil.soln   = peekSpin(qq.soln,   q_info.spin);

					
		     			displacement(u_smr, dil.source, q_info.disp_len, q_info.disp_dir);
		     			displacement(u_smr, dil.soln,   q_info.disp_len, q_info.disp_dir);

						} // for m
					} // for t0
		} // for n
		
		snoop.stop();

		QDPIO::cout << "Displaced Quarks: Spin = "<<key.spin<<" Disp = "<<
			key.displacement <<" Time = "<<snoop.getTimeInSeconds() <<" sec"<<endl;
		// Insert
		disp_quarks.insert(std::make_pair(key, disp_q));

	      } // if find in map
	    } // for i
	  } // for l
	
      cout << __func__ << ": exiting" << endl;

      END_CODE();
    }
    


    //-------------------------------------------------------------------------------
    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "stoch_baryon");
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

      push(xml_out, "stoch_baryon");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << InlineStochGroupBaryonEnv::name << ": Stochastic Baryon Operator" << endl;

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

      // Save current seed
      Seed ran_seed;
      QDP::RNG::savern(ran_seed);

      //
      // Read the source and solutions
      //
      StopWatch swatch, swiss;
      swatch.reset();
      swatch.start();

      multi1d<QuarkSourceSolutions_t>  quarks(params.named_obj.quarks.size());
      QDPIO::cout << "Number of quarks = " << params.named_obj.quarks.size() << endl;

      // Grab the decay direction
      int j_decay;

      try
      {
	QDPIO::cout << "quarks.size = " << quarks.size() << endl;
	

	
	for(int n = 0; n < quarks.size(); ++n)
	{
		bool initq = false;

	  QDPIO::cout << "Attempt to read solutions for source number " << n << endl;
		
		quarks[n].time_dilutions.resize(params.named_obj.quarks[n].time_files.size());

		//Ensure that each time dilution element contains the same number of spin dilutions
		
		int n_spin_dil = params.named_obj.quarks[n].time_files[0].spin_files.size();
		
		
	  QDPIO::cout << "time_dilutions.size = " << quarks[n].time_dilutions.size() << endl;
	  
		for(int t0 = 0 ; t0 < quarks[n].time_dilutions.size() ; ++t0)
	  {
	    bool init_t = false;
			
			if ( n_spin_dil != params.named_obj.quarks[n].time_files[t0].spin_files.size() )
			{
				QDPIO::cerr << "Unequal number of spin dilutions for each time dilution : t0 = " 
					<< t0 << endl; 
				
				QDP_abort(1);

			}

			//Ensure that each additional dilution per spin contains the same number of elements
			int n_dil = params.named_obj.quarks[n].time_files[t0].spin_files[0].dilution_files.size();
		
			
			quarks[n].time_dilutions[t0].spin_dilutions.resize(params.named_obj.quarks[n].time_files[t0].spin_files.size());
			
	    QDPIO::cout << "spin_dilutions.size = " << quarks[n].time_dilutions[t0].spin_dilutions.size() << endl;
	    
			for( int s0 = 0 ; s0 < quarks[n].time_dilutions[t0].spin_dilutions.size() ; ++s0 )
			{
				bool init_s = false;
				
	    	if (n_dil != params.named_obj.quarks[n].time_files[t0].spin_files[s0].dilution_files.size())
				{
					QDPIO::cerr << "Unequal number of additional dilutions for each spin dilution : s0 = "
						<< s0 << endl; 
					
					QDP_abort(1);
				}
				
				quarks[n].time_dilutions[t0].spin_dilutions[s0].dilutions.resize(params.named_obj.quarks[n].time_files[t0].spin_files[s0].dilution_files.size());
			
	    	QDPIO::cout << "dilutions.size = " << quarks[n].time_dilutions[t0].spin_dilutions[s0].dilutions.size() << endl;
			
				for(int i = 0 ; i < quarks[n].time_dilutions[t0].spin_dilutions[s0].dilutions.size() ; ++i)
	    	{
	      	// Short-hand
	      	const std::string& dilution_file = params.named_obj.quarks[n].time_files[t0].spin_files[s0].dilution_files[i];

	      	QuarkSourceSolutions_t::TimeDilutions_t::SpinDilutions_t::Dilutions_t& qq = 
						quarks[n].time_dilutions[t0].spin_dilutions[s0].dilutions[i];

	      	XMLReader file_xml, record_xml;

	      	QDPIO::cout << "reading file= " << dilution_file << endl;
	      	QDPFileReader from(file_xml, dilution_file, QDPIO_SERIAL);

					read(from, record_xml, qq.soln);
					close(from);
	
					read(record_xml, "/Propagator/PropSource", qq.source_header);
	      	read(record_xml, "/Propagator/ForwardProp", qq.prop_header);

					//Get source from the named object map 
					std::stringstream srcstrm;
					srcstrm << 	"zN_source_q" << n + 1 << "_t" << 
					qq.source_header.t_source;

					std::string source_name = srcstrm.str();

					qq.source = TheNamedObjMap::Instance().getData< LatticeFermion >(source_name);
	    
			
					if (!initq)
					{
						read(record_xml, "/Propagator/PropSource/Source/ran_seed",
						quarks[n].seed);
					
	      		j_decay = qq.source_header.j_decay;
						
						initq = true;
					}

					if (!init_t)
					{
						//Must be fixed, assumes full time dilution
						quarks[n].time_dilutions[t0].t0.resize(1);
						quarks[n].time_dilutions[t0].t0[0] = qq.source_header.t_source;
	   			
						init_t = true;
						
					}	


					if (!init_s)
					{
						read(record_xml, "/Propagator/PropSource/Source/spin_mask" , 
								quarks[n].time_dilutions[t0].spin_dilutions[s0].s0 ) ;	
					
					
						init_s = true;
					}
				} //i
	  	} //s0
		} //t0
 	}//n
      } //try
      catch (const string& e) 
      {
	QDPIO::cerr << "Error extracting headers: " << e << endl;
	QDP_abort(1);
      }
      swatch.stop();

      QDPIO::cout << "Sources and solutions successfully read: time= "
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;


      //
      // Initialize the slow Fourier transform phases
      //
      SftMom phases(params.param.mom2_max, false, j_decay);
    
      // Sanity check - if this doesn't work we have serious problems
      if (phases.numSubsets() != QDP::Layout::lattSize()[j_decay])
      {
	QDPIO::cerr << name << ": number of time slices not equal to that in the decay direction: " 
		    << QDP::Layout::lattSize()[j_decay]
		    << endl;
	QDP_abort(1);
      }

      //
      // The first sanity check - the time slices for all the dilutions of each
      // quark must be the same. The dilutions do not have to match, though. 
      // We will keep a list of the participating time slices. This should be
      // the entire time length of the lattice, but this condition may be
      // relaxed in the future. However, we will require that the dilutions
      // for a quark saturate each time slice.
      //
      // To prime the work, grab a first chunk of time slices for the
      // rest of the checks
      //
/*
			multi1d<int> participating_time_slices(quarks[0].time_slices.size());
      for(int t=0; t < quarks[0].time_slices.size(); ++t)
      {
	participating_time_slices[t] = quarks[0].time_slices[t].dilutions[0].t0;
      }
*/
			
//------------------------------------------------------      
/*
#if defined(STOCH_USE_ALL_TIME_SLICES)
      // Sanity check - this may be relaxed later
      // The time slices should be the entire length of the time axis
      if (participating_time_slices.size() != QDP::Layout::lattSize()[j_decay])
      {
	QDPIO::cerr << name << ": number of time slices not equal to that in the decay direction: " 
		    << QDP::Layout::lattSize()[j_decay]
		    << endl;
	QDP_abort(1);
      }
      else
      {
	for(int t=0; t < participating_time_slices.size(); ++t)
	{
	  if (participating_time_slices[t] != t)
	  {
	    QDPIO::cerr << name << ": number of time slices not equal to that in the decay direction: " 
			<< QDP::Layout::lattSize()[j_decay]
			<< endl;
	    QDP_abort(1);
	  }
	}
      }
#endif
*/
//------------------------------------------------------      



      //
      // Check for each quark source that the solutions have their diluted
      // on every site only once
      //
     // swatch.start();

     /* try
      {
	push(xml_out, "Norms"); */
	
			
	/*		
			
	for(int n=0; n < quarks.size(); ++n)
	{
	  
		bool first = true;
	  int  N;
	  LatticeFermion quark_noise;      // noisy source on entire lattice
		
	  // Sanity check - the number of time slices better match
	  
		if (participating_time_slices.size() != quarks[n].time_slices.size())
	  {
	    QDPIO::cerr << name << ": incompatible number of time slices for quark_num= " << n
			<< endl;
	    QDP_abort(1);
	  }
		
	  for(int t=0; t < quarks[n].time_slices.size(); ++t)
	  {
	    for(int i=0; i < quarks[n].time_slices[t].dilutions.size(); ++i)
	    {
	      // Short-hand
	      QuarkSourceSolutions_t::TimeSlices_t::Dilutions_t& qq = 
		quarks[n].time_slices[t].dilutions[i];

	      QDPIO::cout << "Quark_num= " << n << "  time_slice= " << t << "  dilution_num= " << i << endl;

	      // Check time slices
	      if (participating_time_slices[t] != qq.t0)
	      {
		QDPIO::cerr << name << ": time slices incompatible: expected t0="
			    << participating_time_slices[t] 
			    << "  found t=" << qq.t0
			    << endl;
		QDP_abort(1);
	      }

	      // Build source construction
	      QDPIO::cout << "Source_id = " << qq.source_header.source.id << endl;
	      QDPIO::cout << "Source = XX" << qq.source_header.source.xml << "XX" << endl;

	      std::istringstream  xml_s(qq.source_header.source.xml);
	      XMLReader  sourcetop(xml_s);

	      if (qq.source_header.source.id != DiluteZNQuarkSourceConstEnv::name)
	      {
		QDPIO::cerr << "Expected source_type = " << DiluteZNQuarkSourceConstEnv::name << endl;
		QDP_abort(1);
	      }

	      // Manually create the params so I can peek into them and use the source constructor
	      DiluteZNQuarkSourceConstEnv::Params  srcParams(sourcetop, 
							     qq.source_header.source.path);
	      DiluteZNQuarkSourceConstEnv::SourceConst<LatticeFermion>  srcConst(srcParams);
      
	      if (first) 
	      {
		first = false;

		quarks[0].decay_dir = srcParams.j_decay;

		// Grab N
		N = srcParams.N;

		// Set the seed to desired value
		quarks[n].seed = srcParams.ran_seed;
		QDP::RNG::setrn(quarks[n].seed);

		// Create the noisy quark source on the entire lattice
		zN_src(quark_noise, N);
	      }

	      // The seeds must always agree - here the seed is the unique id of the source
	      if ( toBool(srcParams.ran_seed != quarks[n].seed) )
	      {
		QDPIO::cerr << "quark source=" << n << "  dilution=" << i << " seed does not match" << endl;
		QDP_abort(1);
	      }

	      // The N's must always agree
	      if ( toBool(srcParams.N != N) )
	      {
		QDPIO::cerr << "quark source=" << n << "  dilution=" << i << " N does not match" << endl;
		QDP_abort(1);
	      }

	      // Use a trick here, create the source and subtract it from the global noisy
	      // Check at the end that the global noisy is zero everywhere.
	      // NOTE: the seed will be set every call
	      qq.source = srcConst(u);
	      quark_noise -= qq.source;

#if 0
	      // Diagnostic
	      {
		// Keep a copy of the phases with NO momenta
		SftMom phases_nomom(0, true, qq.source_header.j_decay);

		multi1d<Double> source_corr = sumMulti(localNorm2(qq.source), 
						       phases_nomom.getSet());

		multi1d<Double> soln_corr = sumMulti(localNorm2(qq.soln), 
						     phases_nomom.getSet());

		push(xml_out, "elem");
		write(xml_out, "n", n);
		write(xml_out, "i", i);
		write(xml_out, "source_corr", source_corr);
		write(xml_out, "soln_corr", soln_corr);
		pop(xml_out);
	      }
#endif 

			} // end for t
	  } // end for i

	  Double dcnt = norm2(quark_noise);
	  if (toDouble(dcnt) != 0.0)  // problematic - seems to work with unnormalized sources 
	  {
	    QDPIO::cerr << "Noise not saturated by all potential solutions: dcnt=" << dcnt << endl;
	    QDP_abort(1);
	  }

	} // end for n

	pop(xml_out);  // norms
      } // end try
      catch(const std::string& e) 
      {
	QDPIO::cerr << ": Caught Exception creating source: " << e << endl;
	QDP_abort(1);
      }

      swatch.stop();

      QDPIO::cout << "Sources saturated: time= "
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;
*/

		
      //
      // Another sanity check. The seeds of all the quarks must be different
      //
			//Fix this; what if quarks 1 and 2 are the same?
      for(int n=1; n < quarks.size(); ++n)
      {
		 
	if ( toBool(quarks[n].seed == quarks[0].seed) ) 
	{
	  QDPIO::cerr << name << ": error, baryon op seeds are the same" << endl;
	  QDP_abort(1);
	}
      }

      //
      // Start the file writer
      //
      // There is a file XML
      // Then, there are separate records for each creation operator
      // and each annihilation operator
      //
      XMLBufferWriter file_xml;
      push(file_xml, "BaryonOperator");
      write(file_xml, "Params", params.param);
      write(file_xml, "Config_info", gauge_xml);
      pop(file_xml);


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
								       linktop,
								       params.param.link_smearing.path));
	(*linkSmearing)(u_smr);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception link smearing: " << e << endl;
	QDP_abort(1);
      }

			
      MesPlq(xml_out, "Smeared_Observables", u_smr);
      //
      // The spin basis matrix to goto Dirac
      //
      SpinMatrix rotate_mat(adj(DiracToDRMat()));


      //
      // Read operator coefficients
      //
      QDPIO::cout << "Reading 3-quark operators" << endl;
      ThreeQuarkOps_t qqq_oplist; 

	
			const std::string& file = params.named_obj.operators_file.ops_file;

	
			readOps(qqq_oplist, file, params.param.displacement_length);


      //
      // Smear all the quark sources up front, overwrite the originals to save space,
      // and rotate them to the Dirac spin basis
      //
      try
      {
	QDPIO::cout << "Smear all the quark sources up front" << endl;

	// Create the source quark smearing object
	std::istringstream  xml_s(params.param.source_quark_smearing.xml);
	XMLReader  smeartop(xml_s);
	
	Handle< QuarkSmearing<LatticeFermion> > sourceQuarkSmearing =
	  TheFermSmearingFactory::Instance().createObject(params.param.source_quark_smearing.id,
							  smeartop,
							  params.param.source_quark_smearing.path);
	
	// Source smear all the sources up front, multiply by gamma_4 and 
	// go to the Dirac spin basis 
	// NOTE: the creation operator is non-zero on only 1 time slice; however,
	// the smearing functions work on the entire lattice. Oh well.
	for(int n=0; n < quarks.size(); ++n)
	{
	  QDPIO::cout << "Smearing sources for quark[" << n << "]  over time_dilutions = " 
		      << quarks[n].time_dilutions.size() << endl;

	  for(int t0 = 0; t0 < quarks[n].time_dilutions.size(); ++t0)
	  {
			for ( int s0 = 0 ; s0 < quarks[n].time_dilutions[t0].spin_dilutions.size() ; ++s0)
			{
	    	for(int i = 0; i < quarks[n].time_dilutions[t0].spin_dilutions[s0].dilutions.size(); ++i)
	    	{
	      	LatticeFermion src(quarks[n].time_dilutions[t0].spin_dilutions[s0].dilutions[i].source);
	      	(*sourceQuarkSmearing)(src, u_smr);

				//multiply by gamma_4 as well here
				SpinMatrix mat = rotate_mat * Gamma(8);


	      quarks[n].time_dilutions[t0].spin_dilutions[s0].dilutions[i].source = mat * src;
	    	}
	  	}
		}
	}
	QDPIO::cout << "Source smearings done" << endl;
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << ": Caught Exception smearing quark sources: " << e << endl;
	QDP_abort(1);
      }
      catch(...)
      {
	QDPIO::cerr << ": Caught generic exception smearing quark sources" << endl;
	QDP_abort(1);
      }



      //
      // Smear all the quark solutions up front, overwrite the originals to save space,
      // and rotate them to the Dirac spin basis
      //
      try
      {
	QDPIO::cout << "Smear all the quark solutions up front" << endl;

	std::istringstream  xml_s(params.param.sink_quark_smearing.xml);
	XMLReader  smeartop(xml_s);
	
	Handle< QuarkSmearing<LatticeFermion> > sinkQuarkSmearing;
	sinkQuarkSmearing = 
	  TheFermSmearingFactory::Instance().createObject(params.param.sink_quark_smearing.id,
							  smeartop,
							  params.param.sink_quark_smearing.path);

	// Smear the solutions up front
	for(int n=0; n < quarks.size(); ++n)
	{
	  QDPIO::cout << "Smearing solutions for quark[" << n << "]  over time_dilutions = " 
		      << quarks[n].time_dilutions.size() << endl;

	  for(int t0 = 0 ; t0 < quarks[n].time_dilutions.size() ; ++t0)
	  {
			for (int s0 = 0 ; s0 < quarks[n].time_dilutions[t0].spin_dilutions.size() ; ++s0)
			{
				for(int i = 0 ; i < quarks[n].time_dilutions[t0].spin_dilutions[s0].dilutions.size() ; ++i)
	    	{
	      	LatticeFermion soln(quarks[n].time_dilutions[t0].spin_dilutions[s0].dilutions[i].soln);
	      	(*sinkQuarkSmearing)(soln, u_smr);

	      	quarks[n].time_dilutions[t0].spin_dilutions[s0].dilutions[i].soln = rotate_mat * soln;
	    	}
	  	}
		}
	}

	QDPIO::cout << "Solution smearings done" << endl;
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << ": Caught Exception smearing quark solutions: " << e << endl;
	QDP_abort(1);
      }
      catch(...)
      {
	QDPIO::cerr << ": Caught generic exception smearing quark solutions" << endl;
	QDP_abort(1);
      }


      //
      // We want to avoid excessive flops and shifts. So, we want to extract
      // spin components and displace them up front and avoid doing this in
      // a nested loop which will shift the innermost terms many times.
      // However, because we will permute both the source and sink of all quarks,
      // each displacement of a selected spin component will occur for *ALL*
      // quarks, not just say the 1st one. So, we will loop over all the coeff
      // entries and for each quark entry will spin extract and shift *ALL*
      // the quark dilutions since we know they will all eventually be used.
      //
      // This process can be memory intensive. We are intentionally making a 
      // space versus time trade-off.
      //

      // The map of all the displaced spin components
      map<KeySmearedDispColorVector_t, SmearedDispColorVector_t> disp_quarks;

      // The big displacement
      displaceQuarks(disp_quarks, u_smr, qqq_oplist, quarks);



      //
      // Baryon operators
      //
      if (quarks.size() != 3)
      {
	QDPIO::cerr << "expecting 3 quarks but have num quarks= " << quarks.size() << endl;
	QDP_abort(1);
      }

      //
      // Permutations of quarks within an operator
      //
      
			

			const std::string& pstring = params.named_obj.quark_ids; 
			int num_orderings = 1; 
			
			//Should think of a cleverer algorithm for n quarks 
			if  (pstring.size() != 3)
			{
				QDPIO::cerr << "Invalid size for 'quark_ids'. Must be 3 but is " << pstring.size() << endl;
				QDP_abort(1);
			}

			//If 2 identical quarks, different one must be in the third position
			if ( ( pstring[0] == pstring[2] ) && ( pstring[1] != pstring[2] ) )
			{
				QDPIO::cerr << "Invalid format for 'quark_ids'. Identical q's must be last 2 entries."
					<< endl;
				QDP_abort(1);
			}
			
			if ( pstring[1] == pstring[2] )
			{
				num_orderings = 2;
			}
			
			if  (pstring[0] == pstring[2]) 
			{
				num_orderings = 6;
			}
					
      multi1d< multi1d<int> >  perms(num_orderings);
      {
	multi1d<int> p(3);

	if (num_orderings >= 1)
	{
	  p[0] = 0; p[1] = 1; p[2] = 2;
	  perms[0] = p;
	}

	if (num_orderings >= 2)
	{
	  p[0] = 0; p[1] = 2; p[2] = 1;
	  perms[1] = p;
	}

	if (num_orderings >= 3)
	{
	  p[0] = 1; p[1] = 0; p[2] = 2;
	  perms[2] = p;
	}

	if (num_orderings >= 4)
	{
	  p[0] = 1; p[1] = 2; p[2] = 0;
	  perms[3] = p;
	}

	if (num_orderings >= 5)
	{
	  p[0] = 2; p[1] = 1; p[2] = 0;
	  perms[4] = p;
	}
 
	if (num_orderings >= 6)
	{
	  p[0] = 2; p[1] = 0; p[2] = 1;
	  perms[5] = p;
	}
      }

      // Creation operator
      BaryonOperator_t  creat_oper;
      creat_oper.mom2_max    = params.param.mom2_max;
      creat_oper.decay_dir   = j_decay;
      creat_oper.seed_l      = quarks[0].seed;
      creat_oper.seed_m      = quarks[1].seed;
      creat_oper.seed_r      = quarks[2].seed;
      creat_oper.smearing    = params.param.source_quark_smearing;
      creat_oper.perms       = perms;
      creat_oper.orderings.resize(num_orderings);

			 // Annihilation operator
      BaryonOperator_t  annih_oper;
      annih_oper.mom2_max    = params.param.mom2_max;
      annih_oper.decay_dir   = j_decay;
      annih_oper.seed_l      = quarks[0].seed;
      annih_oper.seed_m      = quarks[1].seed;
      annih_oper.seed_r      = quarks[2].seed;
      annih_oper.smearing    = params.param.sink_quark_smearing;
      annih_oper.perms       = perms;
      annih_oper.orderings.resize(num_orderings);


      // Construct creation operator
      QDPIO::cout << "Build operators" << endl;



	//Who's the finest quark of them all?
		int finest_num = 0, time_grid_size = Layout::lattSize()[j_decay]; 

		for (int n = 0 ; n < quarks.size() ; ++n)
		{
			//assume each time dilution element contains the same # of timeslices
			if (quarks[n].time_dilutions[0].t0.size() < time_grid_size)
			{
				time_grid_size = quarks[n].time_dilutions[0].t0.size();
				finest_num = n; 
			}

		}

		const QuarkSourceSolutions_t& finest_quark = quarks[finest_num];

		
			
	// Loop over each operator 
	for(int l=0; l < qqq_oplist.ops.size(); ++l)
	{
	  //QDPIO::cout << "Creation operator: op = " << l << endl;

    push(xml_out, "BaryonOperator");

		creat_oper.id = qqq_oplist.ops[l].name;
		annih_oper.id = qqq_oplist.ops[l].name;
	

		write(xml_out, "Name", creat_oper.id);

	  // The operator number. This is just the operator number within the 
	  // coefficient file
	  creat_oper.operator_num = l;
	  annih_oper.operator_num = l;

	  // Loop over all orderings and build the operator
	  swiss.reset();
	  swiss.start();

	  for(int ord = 0; ord < creat_oper.orderings.size(); ++ord)
	  {
	    QDPIO::cout << "Ordering = " << ord << endl;

	    creat_oper.orderings[ord].perm = perms[ord];
	    creat_oper.orderings[ord].time_dilutions.resize(finest_quark.time_dilutions.size());
		
	    annih_oper.orderings[ord].perm = perms[ord];
	    annih_oper.orderings[ord].time_dilutions.resize(finest_quark.time_dilutions.size());
     
	    const int n0 = perms[ord][0];
	    const int n1 = perms[ord][1];
	    const int n2 = perms[ord][2];

	    
			//Fetch the quarks

		  const ThreeQuarkOps_t::ThreeQuarkOp_t& qqq_op = qqq_oplist.ops[l];
		     
			KeySmearedDispColorVector_t key0;
		  key0.displacement = qqq_op.quark[0].displacement;
		  key0.spin         = qqq_op.quark[0].spin;
	      
		  KeySmearedDispColorVector_t key1;
		  key1.displacement = qqq_op.quark[1].displacement;
		  key1.spin         = qqq_op.quark[1].spin;
	      
		  KeySmearedDispColorVector_t key2;
		  key2.displacement = qqq_op.quark[2].displacement;
		  key2.spin         = qqq_op.quark[2].spin;

		  // Sanity check - this entry better be in the map or blow-up
		  if (disp_quarks.find(key0) == disp_quarks.end())
		  {
				QDPIO::cerr << name 
				    << ": internal error - map of displacements missing an entry" 
				    << endl;
				QDP_abort(1);
		  }
		      
		  // The location of the displaced quark
		  const SmearedDispColorVector_t& disp_q0 = disp_quarks.find(key0)->second;
		  const SmearedDispColorVector_t& disp_q1 = disp_quarks.find(key1)->second;
		  const SmearedDispColorVector_t& disp_q2 = disp_quarks.find(key2)->second;
		
			
			//Loop over the time dilutions included by every quark		
	    for(int t0 = 0; t0 < finest_quark.time_dilutions.size(); ++t0)
	    {
				QDPIO::cout<<" evaluation of time dilution element : "<< t0 <<endl; 
				
				BaryonOperator_t::Orderings_t::TimeDilutions_t& cop = 
					creat_oper.orderings[ord].time_dilutions[ t0 ];

	      BaryonOperator_t::Orderings_t::TimeDilutions_t& aop = 
					annih_oper.orderings[ord].time_dilutions[ t0 ];
				

				//time slices included in this dilution element for the finest quark
				const multi1d<int>& times = finest_quark.time_dilutions[t0].t0;
				
				//Loop over the timeslices included in this dilution element
				for (int time = 0 ; time < times.size() ; ++time)
				{
					QDPIO::cout << "time = " <<time << endl;
					//Determine in which dilution element time0 occurs
					const int timeDil0 = quarks[n0].findTime( times[ time ] );
					const int timeDil1 = quarks[n1].findTime( times[ time ] );
					const int timeDil2 = quarks[n2].findTime( times[ time ] );
				

					//Determine in which dilution elements the spins occur
					const int spinDil0 = quarks[n0].time_dilutions[timeDil0].findSpin( key0.spin );
					const int spinDil1 = quarks[n1].time_dilutions[timeDil1].findSpin( key1.spin );
					const int spinDil2 = quarks[n2].time_dilutions[timeDil2].findSpin( key2.spin );
				
				
					//Store which time_slices participate in this dilution element for each quark
	      	cop.t0_l = quarks[n0].time_dilutions[ timeDil0 ].t0;
	      	cop.t0_m = quarks[n1].time_dilutions[ timeDil1 ].t0;
	      	cop.t0_r = quarks[n2].time_dilutions[ timeDil2 ].t0;
				
	      	aop.t0_l = quarks[n0].time_dilutions[ timeDil0 ].t0;
	      	aop.t0_m = quarks[n1].time_dilutions[ timeDil1 ].t0;
	      	aop.t0_r = quarks[n2].time_dilutions[ timeDil2 ].t0;
				
					//Store which spins participate in this dilution element for each quark
			  	cop.spin.s0_l =	quarks[n0].time_dilutions[ timeDil0 ].spin_dilutions[ spinDil0 ].s0;
				 	cop.spin.s0_m =	quarks[n1].time_dilutions[ timeDil1 ].spin_dilutions[ spinDil1 ].s0;
			    cop.spin.s0_r = quarks[n2].time_dilutions[ timeDil2 ].spin_dilutions[ spinDil2 ].s0;

			  	aop.spin.s0_l = quarks[n0].time_dilutions[ timeDil0 ].spin_dilutions[ spinDil0 ].s0;
				 	aop.spin.s0_m = quarks[n1].time_dilutions[ timeDil1 ].spin_dilutions[ spinDil1 ].s0;
			    aop.spin.s0_r = quarks[n2].time_dilutions[ timeDil2 ].spin_dilutions[ spinDil2 ].s0;

				
				
					cop.spin.dilutions.resize(quarks[n0].time_dilutions[ timeDil0 ].spin_dilutions[ spinDil0 ].dilutions.size(), 
				   	quarks[n1].time_dilutions[ timeDil1 ].spin_dilutions[ spinDil1 ].dilutions.size(), 
				   	quarks[n2].time_dilutions[ timeDil2 ].spin_dilutions[ spinDil2 ].dilutions.size() );

					aop.spin.dilutions.resize(quarks[n0].time_dilutions[ timeDil0 ].spin_dilutions[ spinDil0 ].dilutions.size(), 
				   	quarks[n1].time_dilutions[ timeDil1 ].spin_dilutions[ spinDil1 ].dilutions.size(), 
				   	quarks[n2].time_dilutions[ timeDil2 ].spin_dilutions[ spinDil2 ].dilutions.size() );

				
	      	for( int i = 0 ; i < cop.spin.dilutions.size1() ; ++i )
	      	{
			for( int j = 0 ; j < cop.spin.dilutions.size2() ; ++j )
			{
		  	for( int k = 0 ; k < cop.spin.dilutions.size3() ; ++k)
		  	{

		      // Contract over color indices with antisym tensor.
		      // There is a potential optimization here - the colorcontract of
		      // the first two quarks could be pulled outside the innermost dilution
		      // loop.
		      // NOTE: the creation operator only lives on a few time slice, so restrict
		      // the operation to that time slice
				
					LatticeComplex c_oper;
		      c_oper[phases.getSet()[ times[time] ] ] =
			colorContract(disp_q0.quarks[n0].time_dilutions[ timeDil0 ].spin_dilutions.dilutions[i].source,
				      disp_q1.quarks[n1].time_dilutions[ timeDil1 ].spin_dilutions.dilutions[j].source,
				      disp_q2.quarks[n2].time_dilutions[ timeDil2 ].spin_dilutions.dilutions[k].source);
		   
					//This operation will be performed too much for non-maximal time dilution 
					LatticeComplex a_oper = 
						colorContract(disp_q0.quarks[n0].time_dilutions[ timeDil0 ].spin_dilutions.dilutions[i].soln,
				      disp_q1.quarks[n1].time_dilutions[ timeDil1 ].spin_dilutions.dilutions[j].soln,
				      disp_q2.quarks[n2].time_dilutions[ timeDil2 ].spin_dilutions.dilutions[k].soln);
					
		    	// Slow fourier-transform
		    	// We can restrict what the FT routine requires to a subset.
		    
					multi2d<DComplex> c_sum( phases.sft(c_oper, times[time]) );
					multi2d<DComplex> a_sum( phases.sft(a_oper) );

		    	// Unpack into separate momentum and correlator
		    	cop.spin.dilutions(i,j,k).mom_projs.resize(phases.numMom());
		    	aop.spin.dilutions(i,j,k).mom_projs.resize(phases.numMom());
		    
					for(int mom_num = 0 ; mom_num < phases.numMom() ; ++mom_num) 
		    	{
		      	cop.spin.dilutions(i,j,k).mom_projs[mom_num].mom = phases.numToMom(mom_num);
		      	aop.spin.dilutions(i,j,k).mom_projs[mom_num].mom = phases.numToMom(mom_num);
		      
						cop.spin.dilutions(i,j,k).mom_projs[mom_num].op[ time ] = c_sum[mom_num][ times[time] ];
						aop.spin.dilutions(i,j,k).mom_projs[mom_num].op = a_sum[mom_num];
		    	}
		  	} // end for k
			} // end for j
	      	} // end for i
	    
				} //time
			} // end for t0
	  } // end for ord

		

	  swiss.stop();

		
	  QDPIO::cout << "Operator construction: operator= " << l 
		      << "  time= "
		      << swiss.getTimeInSeconds() 
		      << " secs" << endl;
/*
for (int pr = 0 ; pr < 6 ; ++pr)
{
		QDPIO::cout<<"Source testval(p = " << pr << " ) = "<<creat_oper.orderings[pr].time_slices[0].dilutions(0,0,0).mom_projs[0].op[0]<<endl;
	//	QDPIO::cout<<"Source SD testval = "<<baryon_ops_Source.ops[1].orderings[0].time_slices[6].dilutions(0,0,0).mom_projs[0].op[0]<<endl;
}
*/
					//Hard code the elemental op name for now 
			
			std::stringstream cnvrt;

			cnvrt <<  creat_oper.id  << ".lime";

		
			std::string filename;

			filename = cnvrt.str(); 

			QDPFileWriter qdp_file(file_xml, filename,     // are there one or two files???
			     QDPIO_SINGLEFILE, QDPIO_SERIAL, QDPIO_OPEN);

	  // Write the meta-data and the binary for this operator
	  swiss.reset();
		swiss.start();
	  {
	    XMLBufferWriter     src_record_xml, snk_record_xml, file_xml;
	    BinaryBufferWriter  src_record_bin, snk_record_bin;

      push(file_xml, "BaryonOperator");
      write(file_xml, "Params", params.param);
      write(file_xml, "Config_info", gauge_xml);
      pop(file_xml);

		  write(src_record_xml, "BaryonCreationOperator", creat_oper);
	    write(src_record_bin, creat_oper);

	    write(qdp_file, src_record_xml, src_record_bin);
	  
	  

	    write(snk_record_xml, "BaryonAnnihilationOperator", annih_oper);
	    write(snk_record_bin, annih_oper);

	    write(qdp_file, snk_record_xml, snk_record_bin);
	  }
	  swiss.stop();

	  QDPIO::cout << "Operator writing: operator = " << l 
		      << "  time= "
		      << swiss.getTimeInSeconds() 
		      << " secs" << endl;
	
		
		 pop(xml_out); // Operator 

		
  } // end for l (operator )



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

  /*! @} */  // end of group hadron

} // namespace Chroma
