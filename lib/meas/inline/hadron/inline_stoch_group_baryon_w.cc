// $Id: inline_stoch_group_baryon_w.cc,v 1.7 2007-12-06 18:31:58 jbulava Exp $
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
#include "meas/sources/dilutezN_source_const.h"
#include "meas/sources/zN_src.h"
#include "meas/glue/mesplq.h"
#include "meas/smear/displacement.h"
#include "util/ferm/diractodr.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"
#include <sstream> 

#include "meas/inline/io/named_objmap.h"

//#define STOCH_USE_ALL_TIME_SLICES 1

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


    // Reader for input parameters
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


    // Coefficient files
    void read(XMLReader& xml, const string& path, InlineStochGroupBaryonEnv::Params::NamedObject_t::CoeffFiles_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "coeff_file", input.coeff_file);
      read(inputtop, "id", input.id);
    }


    // Coefficient files
    void write(XMLWriter& xml, const string& path, const InlineStochGroupBaryonEnv::Params::NamedObject_t::CoeffFiles_t& input)
    {
      push(xml, path);
      write(xml, "coeff_file", input.coeff_file);
      write(xml, "id", input.id);
      pop(xml);
    }


    // Dilutions for each time slice
    void read(XMLReader& xml, const string& path, InlineStochGroupBaryonEnv::Params::NamedObject_t::Operator_t::TimeSlices_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "dilution_files", input.dilution_files);
    }


    // Dilutions for each time slice
    void write(XMLWriter& xml, const string& path, const InlineStochGroupBaryonEnv::Params::NamedObject_t::Operator_t::TimeSlices_t& input)
    {
      push(xml, path);
      write(xml, "dilution_files", input.dilution_files);
      pop(xml);
    }


    // Operator parameters
    void read(XMLReader& xml, const string& path, InlineStochGroupBaryonEnv::Params::NamedObject_t::Operator_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "soln_files", input.soln_files);
    }


    // Operator parameters
    void write(XMLWriter& xml, const string& path, const InlineStochGroupBaryonEnv::Params::NamedObject_t::Operator_t& input)
    {
      push(xml, path);
      write(xml, "soln_files", input.soln_files);
      pop(xml);
    }


    //! Propagator parameters
    void read(XMLReader& xml, const string& path, InlineStochGroupBaryonEnv::Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "operator_file", input.operator_file);
      read(inputtop, "operator_coeff_files", input.operator_coeff_files);
      read(inputtop, "Quarks", input.quarks);
    }

    //! Propagator parameters
    void write(XMLWriter& xml, const string& path, const InlineStochGroupBaryonEnv::Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "operator_file", input.operator_file);
      write(xml, "operator_coeff_files", input.operator_coeff_files);
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
      struct TimeSlices_t
      {
	struct Dilutions_t
	{
	  int                t0;
	  LatticeFermion     source;
	  LatticeFermion     soln;
	  PropSourceConst_t  source_header;
	  ChromaProp_t       prop_header;
	};

	multi1d<Dilutions_t>  dilutions;
      };

      int   decay_dir;
      Seed  seed;
      multi1d<TimeSlices_t>  time_slices;
    };


    //! Operator coefficient structure
    struct OperCoeffs_t
    {
      struct CoeffTerms_t
      {
	struct CoeffTerm_t
	{
	  struct QuarkTerm_t
	  {
	    int  displacement;    /*!< Orig plus/minus 1-based directional displacements */

	    int  spin;            /*!< 0-based spin index */
	    int  disp_dir;        /*!< 0-based direction */
	    int  disp_len;        /*!< 0-based length */
	  };

	  multi1d<QuarkTerm_t>  quark;    /*!< Displacement and spin for each quark */
	  Complex               coeff;    /*!< Weight on color contraction */
	};

	multi1d<CoeffTerm_t> op;          /*!< Terms within a single operator */
	std::string name;                 /*!< Name of the operator */
			};

      multi1d<CoeffTerms_t>  ops;         /*!< All the operators within a file */
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
	struct TimeSlices_t
	{
	  struct Dilutions_t
	  {
	    int                 t0;
	    LatticeColorVector  source;
	    LatticeColorVector  soln;
	  };

	  multi1d<Dilutions_t> dilutions;
	};

	multi1d<TimeSlices_t> time_slices;
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
	struct TimeSlices_t
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

	  int                  t0;          /*!< Source time location */
	  multi3d<Dilutions_t> dilutions;   /*!< Hybrid list indices */
	};
	  
	multi1d<int> perm;                  /*!< This particular permutation of quark orderings */
	multi1d<TimeSlices_t> time_slices;  /*!< Time slices of the lattice that are used */
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
    void write(BinaryWriter& bin, const BaryonOperator_t::Orderings_t::TimeSlices_t::Dilutions_t::Mom_t& param)
    {
      write(bin, param.mom);
      write(bin, param.op);
    }

    //! BaryonOperator binary writer
    void write(BinaryWriter& bin, const BaryonOperator_t::Orderings_t::TimeSlices_t::Dilutions_t& param)
    {
      write(bin, param.mom_projs);
    }

    //! BaryonOperator binary writer
    void write(BinaryWriter& bin, const BaryonOperator_t::Orderings_t::TimeSlices_t& param)
    {
      write(bin, param.dilutions);
    }

    //! BaryonOperator binary writer
    void write(BinaryWriter& bin, const BaryonOperator_t::Orderings_t& param)
    {
      write(bin, param.perm);
      write(bin, param.time_slices);
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




    //! Read operator coefficient files
    void readCoeffs(OperCoeffs_t& coeffs, 
		    const std::string& operator_coeff_file,
		    int displacement_length)
    {
      START_CODE();

      TextFileReader reader(operator_coeff_file);

      int num_ops;
      reader >> num_ops;
      coeffs.ops.resize(num_ops);

      for(int n=0; n < coeffs.ops.size(); ++n)
      {
	int num_terms;
	std::string name; 
	reader >> num_terms >> name;
	coeffs.ops[n].op.resize(num_terms);
	coeffs.ops[n].name = name;  
	
	for(int l=0; l < coeffs.ops[n].op.size(); ++l)
	{
	  OperCoeffs_t::CoeffTerms_t::CoeffTerm_t& term = coeffs.ops[n].op[l];
	  term.quark.resize(3);

	  // Make spin index 0 based
	  {
	    multi1d<int> spin(3);
	    reader >> spin[0] >> spin[1] >> spin[2];

	    for(int i=0; i < spin.size(); ++i)
	      term.quark[i].spin = spin[i] - 1;
	  }

	  // Convert displacements to  disp_dir, disp_len
	  {
	    multi1d<int> displacement(3);
	    reader >> displacement[0] >> displacement[1] >> displacement[2];

	    for(int i=0; i < displacement.size(); ++i)
	    {
	      term.quark[i].displacement = displacement[i];

	      if (displacement[i] == 0)
	      {
		term.quark[i].disp_dir = 0;
		term.quark[i].disp_len = 0;
	      }
	      else if (displacement[i] > 0)
	      {
		term.quark[i].disp_dir = term.quark[i].displacement - 1;
		term.quark[i].disp_len = displacement_length;
	      }
	      else
	      {
		term.quark[i].disp_dir = -term.quark[i].displacement - 1;
		term.quark[i].disp_len = -displacement_length;
	      }
	    }
	  }

	  // Read the garbage around a complex
	  {
	    Real re, im;
	    char lparen, comma, rparen;

	    reader >> lparen >> re >> comma >> im >> rparen;
	    term.coeff = cmplx(re,im);
	  }
	}
      }
	
      reader.close();

      END_CODE();
    }


    //! Construct array of maps of displacements
    void displaceQuarks(map<KeySmearedDispColorVector_t, SmearedDispColorVector_t>& disp_quarks,
			const multi1d<LatticeColorMatrix>& u_smr,
			const multi1d<OperCoeffs_t>&  coeffs,
			const multi1d<QuarkSourceSolutions_t> quarks)
    {
      START_CODE();

			StopWatch snoop;
			
      cout << __func__ << ": entering" << endl;

      // Loop over all files of operators
      for(int f=0; f < coeffs.size(); ++f)
      {
	// Loop over each operator within a file
	for(int c=0; c < coeffs[f].ops.size(); ++c)
	{
	  // Loop over the rows/terms within an operator
	  for(int l=0; l < coeffs[f].ops[c].op.size(); ++l)
	  {
	    const OperCoeffs_t::CoeffTerms_t::CoeffTerm_t& term = coeffs[f].ops[c].op[l];
	    
	    for(int i=0; i < term.quark.size(); ++i)
	    {
	      const OperCoeffs_t::CoeffTerms_t::CoeffTerm_t::QuarkTerm_t& term_q = term.quark[i];

	      // Check if this displacement and spin are in the map,
	      // if not, then *ALL* quarks and dilutions must be shifted
	      KeySmearedDispColorVector_t key;
	      key.displacement = term_q.displacement;
	      key.spin         = term_q.spin;

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
		
		for(int n=0; n < disp_q.quarks.size(); ++n)
		{
		  disp_q.quarks[n].time_slices.resize(quarks[n].time_slices.size());

		  for(int t=0; t < disp_q.quarks[n].time_slices.size(); ++t)
		  {
		    disp_q.quarks[n].time_slices[t].dilutions.resize(quarks[n].time_slices[t].dilutions.size());

		    for(int m=0; m < disp_q.quarks[n].time_slices[t].dilutions.size(); ++m)
		    {
		      // Short-hand
		      const QuarkSourceSolutions_t::TimeSlices_t::Dilutions_t& qq = 
			quarks[n].time_slices[t].dilutions[m];

		      SmearedDispColorVector_t::Quarks_t::TimeSlices_t::Dilutions_t& dil = 
			disp_q.quarks[n].time_slices[t].dilutions[m];
					

		      // Pull out the appropriate spin component, then displace it
		      dil.t0     = qq.t0;
		      dil.source = peekSpin(qq.source, term_q.spin);
		      dil.soln   = peekSpin(qq.soln,   term_q.spin);

					
					//if (dil.source != zero)
					{
		      	displacement(u_smr, dil.source, term_q.disp_len, term_q.disp_dir);
		      	displacement(u_smr, dil.soln,   term_q.disp_len, term_q.disp_dir);
		    	}

					} // for m
		  } // for t
		} // for n
		
		snoop.stop();

		QDPIO::cout << "Displaced Quarks: Spin = "<<key.spin<<" Disp = "<<
			key.displacement <<" Time = "<<snoop.getTimeInSeconds() <<" sec"<<endl;
		// Insert
		disp_quarks.insert(std::make_pair(key, disp_q));

	      } // if find in map
	    } // for i
	  } // for l
	}  // for n
      }  // for c
	
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
      QDPIO::cout << "Number of quarks= " << params.named_obj.quarks.size() << endl;

      // Grab the decay direction
      int j_decay;

      try
      {
	QDPIO::cout << "quarks.size= " << quarks.size() << endl;
	for(int n=0; n < quarks.size(); ++n)
	{
		bool initq = false;

	  QDPIO::cout << "Attempt to read solutions for source number= " << n << endl;
	  quarks[n].time_slices.resize(params.named_obj.quarks[n].soln_files.size());

	  QDPIO::cout << "time_slices.size= " << quarks[n].time_slices.size() << endl;
	  for(int t=0; t < quarks[n].time_slices.size(); ++t)
	  {
	    quarks[n].time_slices[t].dilutions.resize(params.named_obj.quarks[n].soln_files[t].dilution_files.size());
	    QDPIO::cout << "dilutions.size= " << quarks[n].time_slices[t].dilutions.size() << endl;
	    for(int i=0; i < quarks[n].time_slices[t].dilutions.size(); ++i)
	    {
	      // Short-hand
	      const std::string& dilution_file = params.named_obj.quarks[n].soln_files[t].dilution_files[i];

	      QuarkSourceSolutions_t::TimeSlices_t::Dilutions_t& qq = 
		quarks[n].time_slices[t].dilutions[i];

	      XMLReader file_xml, record_xml, record_xml_source;

	      QDPIO::cout << "reading file= " << dilution_file << endl;
	      QDPFileReader from(file_xml, dilution_file, QDPIO_SERIAL);

				//For now, read both source and solution
				read(from, record_xml, qq.soln);
	    //  read(from, record_xml_source, qq.source);
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
					
					initq = true;
				}

	      qq.t0 = qq.source_header.t_source;
	      j_decay = qq.source_header.j_decay;
	   
			}
	  }
	}
      }
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
      multi1d<int> participating_time_slices(quarks[0].time_slices.size());
      for(int t=0; t < quarks[0].time_slices.size(); ++t)
      {
	participating_time_slices[t] = quarks[0].time_slices[t].dilutions[0].t0;
      }

//------------------------------------------------------      
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
//------------------------------------------------------      



      //
      // Check for each quark source that the solutions have their diluted
      // on every site only once
      //
     // swatch.start();

     /* try
      {
	push(xml_out, "Norms"); */
	for(int n=0; n < quarks.size(); ++n)
	{
	 /* bool first = true;
	  int  N;
	  LatticeFermion quark_noise;      // noisy source on entire lattice
*/
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
/*
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
*/
			} // end for t
	  } // end for i
/*
	  Double dcnt = norm2(quark_noise);
	  if (toDouble(dcnt) != 0.0)  // problematic - seems to work with unnormalized sources 
	  {
	    QDPIO::cerr << "Noise not saturated by all potential solutions: dcnt=" << dcnt << endl;
	    QDP_abort(1);
	  }
*/
	} // end for n
/*
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

//      QDPFileWriter qdp_file(file_xml, params.named_obj.operator_file,     // are there one or two files???
	//		     QDPIO_SINGLEFILE, QDPIO_SERIAL, QDPIO_OPEN);

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
      QDPIO::cout << "Reading operator coefficient files" << endl;
      multi1d<OperCoeffs_t>  coeffs(params.named_obj.operator_coeff_files.size());

      for(int i=0; i < coeffs.size(); ++i)
      {
	const std::string& file = params.named_obj.operator_coeff_files[i].coeff_file;
	QDPIO::cout << "Reading operator coefficient file = " << file << endl;

	readCoeffs(coeffs[i], file, params.param.displacement_length);
      }


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
	
	// Source smear all the sources up front
	// NOTE: the creation operator is non-zero on only 1 time slice; however,
	// the smearing functions work on the entire lattice. Oh well.
	for(int n=0; n < quarks.size(); ++n)
	{
	  QDPIO::cout << "Smearing sources for quark[" << n << "]  over dilutions = " 
		      << quarks[n].time_slices.size() << endl;

	  for(int t=0; t < quarks[n].time_slices.size(); ++t)
	  {
	    for(int i=0; i < quarks[n].time_slices[t].dilutions.size(); ++i)
	    {
	      LatticeFermion src(quarks[n].time_slices[t].dilutions[i].source);
	      (*sourceQuarkSmearing)(src, u_smr);

				//multiply by gamma_4 as well here
	      quarks[n].time_slices[t].dilutions[i].source = rotate_mat * Gamma(1) * src;
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
	  QDPIO::cout << "Smearing solutions for quark[" << n << "]  over dilutions = " 
		      << quarks[n].time_slices.size() << endl;

	  for(int t=0; t < quarks[n].time_slices.size(); ++t)
	  {
	    for(int i=0; i < quarks[n].time_slices[t].dilutions.size(); ++i)
	    {
	      LatticeFermion soln(quarks[n].time_slices[t].dilutions[i].soln);
	      (*sinkQuarkSmearing)(soln, u_smr);

	      quarks[n].time_slices[t].dilutions[i].soln = rotate_mat * soln;
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
      displaceQuarks(disp_quarks, u_smr, coeffs, quarks);



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
      int num_orderings = 6;   // number of permutations of the numbers  0,1,2
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
     // swatch.start();
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
     // swatch.start();
      BaryonOperator_t  annih_oper;
      annih_oper.mom2_max    = params.param.mom2_max;
      annih_oper.decay_dir   = j_decay;
      annih_oper.seed_l      = quarks[0].seed;
      annih_oper.seed_m      = quarks[1].seed;
      annih_oper.seed_r      = quarks[2].seed;
      annih_oper.smearing    = params.param.sink_quark_smearing;
      annih_oper.perms       = perms;
      annih_oper.orderings.resize(num_orderings);

     // push(xml_out, "BaryonCreationOperator");

      // Construct creation operator
      QDPIO::cout << "Build operators" << endl;

      // Loop over all files of operators
      for(int f=0; f < coeffs.size(); ++f)
      {
	// Set the id to be used in the analysis codes
//	creat_oper.id = params.named_obj.operator_coeff_files[f].id; 

	// Loop over each operator within a file
	for(int c=0; c < coeffs[f].ops.size(); ++c)
	{
	  QDPIO::cout << "Creation operator: f=" << f << "  op= " << c << endl;

    push(xml_out, "BaryonCreationOperator");

		creat_oper.id = coeffs[f].ops[c].name;

		write(xml_out, "Name", creat_oper.id);

	  // The operator number. This is just the operator number within the 
	  // coefficient file
	  creat_oper.operator_num = c;

	  // Loop over all orderings and build the operator
	  swiss.reset();
	  swiss.start();

	  for(int ord=0; ord < creat_oper.orderings.size(); ++ord)
	  {
	    QDPIO::cout << "Creation operator: ordering = " << ord << endl;

	    creat_oper.orderings[ord].perm = perms[ord];
	    creat_oper.orderings[ord].time_slices.resize(participating_time_slices.size());
     
	    const int n0 = perms[ord][0];
	    const int n1 = perms[ord][1];
	    const int n2 = perms[ord][2];

	    // The operator must hold all the dilutions
	    // We know that all time slices match. However, not all time slices of the
	    // lattice maybe used
	    for(int t=0; t < participating_time_slices.size(); ++t)
	    {
	      BaryonOperator_t::Orderings_t::TimeSlices_t& cop = creat_oper.orderings[ord].time_slices[t];

	      cop.t0 = participating_time_slices[t];
	      cop.dilutions.resize(quarks[n0].time_slices[t].dilutions.size(), 
				   quarks[n1].time_slices[t].dilutions.size(), 
				   quarks[n2].time_slices[t].dilutions.size());

	      for(int i=0; i < quarks[n0].time_slices[t].dilutions.size(); ++i)
	      {
		for(int j=0; j < quarks[n1].time_slices[t].dilutions.size(); ++j)
		{
		  for(int k=0; k < quarks[n2].time_slices[t].dilutions.size(); ++k)
		  {
		    // The correlator
		    LatticeComplex bar = zero;

		    // Loop over the rows/terms within an operator
		    for(int l=0; l < coeffs[f].ops[c].op.size(); ++l)
		    {
		      const OperCoeffs_t::CoeffTerms_t::CoeffTerm_t& term_q = coeffs[f].ops[c].op[l];
		     
					KeySmearedDispColorVector_t key0;
		      key0.displacement = term_q.quark[0].displacement;
		      key0.spin         = term_q.quark[0].spin;
	      
		      KeySmearedDispColorVector_t key1;
		      key1.displacement = term_q.quark[1].displacement;
		      key1.spin         = term_q.quark[1].spin;
	      
		      KeySmearedDispColorVector_t key2;
		      key2.displacement = term_q.quark[2].displacement;
		      key2.spin         = term_q.quark[2].spin;


					/* 
		      KeySmearedDispColorVector_t key0;
		      key0.displacement = term_q.quark[n0].displacement;
		      key0.spin         = term_q.quark[n0].spin;
	      
		      KeySmearedDispColorVector_t key1;
		      key1.displacement = term_q.quark[n1].displacement;
		      key1.spin         = term_q.quark[n1].spin;
	      
		      KeySmearedDispColorVector_t key2;
		      key2.displacement = term_q.quark[n2].displacement;
		      key2.spin         = term_q.quark[n2].spin;
	      */
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

		      // Contract over color indices with antisym tensor.
		      // There is a potential optimization here - the colorcontract of
		      // the first two quarks could be pulled outside the innermost dilution
		      // loop.
		      // NOTE: the creation operator only lives on 1 time slice, so restrict
		      // the operation to that time slice
		      LatticeComplex b_oper;
		      b_oper[phases.getSet()[cop.t0]] =
			colorContract(disp_q0.quarks[n0].time_slices[t].dilutions[i].source,
				      disp_q1.quarks[n1].time_slices[t].dilutions[j].source,
				      disp_q2.quarks[n2].time_slices[t].dilutions[k].source);
		    
		      bar[phases.getSet()[cop.t0]] += term_q.coeff * b_oper;
		    } // end for l

		    // Slow fourier-transform
		    // We can restrict the FT routine requires to the t0 subset.
		    multi2d<DComplex> hsum(phases.sft(bar, cop.t0));

		    // Unpack into separate momentum and correlator
		    // NOTE: because creation operator is only non-zero on 1 time slice, we only
		    // keep that part
		    cop.dilutions(i,j,k).mom_projs.resize(phases.numMom());
		    for(int sink_mom_num=0; sink_mom_num < phases.numMom(); ++sink_mom_num) 
		    {
		      cop.dilutions(i,j,k).mom_projs[sink_mom_num].mom   = phases.numToMom(sink_mom_num);
		      cop.dilutions(i,j,k).mom_projs[sink_mom_num].op.resize(1);
		      cop.dilutions(i,j,k).mom_projs[sink_mom_num].op[0] = hsum[sink_mom_num][cop.t0];
		    }
		  } // end for k
		} // end for j
	      } // end for i
	    } // end for t
	  } // end for ord

		

	  swiss.stop();

	  QDPIO::cout << "Creation operator construction: file= " << f 
		      << "  operator= " << c 
		      << "  time= "
		      << swiss.getTimeInSeconds() 
		      << " secs" << endl;

for (int pr = 0 ; pr < 6 ; ++pr)
{
		QDPIO::cout<<"Source testval(p = " << pr << " ) = "<<creat_oper.orderings[pr].time_slices[0].dilutions(0,0,0).mom_projs[0].op[0]<<endl;
	//	QDPIO::cout<<"Source SD testval = "<<baryon_ops_Source.ops[1].orderings[0].time_slices[6].dilutions(0,0,0).mom_projs[0].op[0]<<endl;
}
	//Hard code the elemental op name for now 
			
			std::stringstream cnvrt;

			cnvrt << creat_oper.id << "_" << 
				creat_oper.orderings[0].time_slices[0].t0 << ".lime";

		
			std::string filename;

			filename = cnvrt.str(); 

			QDPFileWriter qdp_file(file_xml, filename,     // are there one or two files???
			     QDPIO_SINGLEFILE, QDPIO_SERIAL, QDPIO_OPEN);

	  // Write the meta-data and the binary for this operator
	  swiss.reset();
		swiss.start();
	  {
	    XMLBufferWriter     record_xml, file_xml;
	    BinaryBufferWriter  record_bin;

      push(file_xml, "BaryonOperator");
      write(file_xml, "Params", params.param);
      write(file_xml, "Config_info", gauge_xml);
      pop(file_xml);

		  write(record_xml, "BaryonCreationOperator", creat_oper);
	    write(record_bin, creat_oper);

	    write(qdp_file, record_xml, record_bin);
	  }
	  swiss.stop();

	  QDPIO::cout << "Creation operator writing: file= " << f 
		      << "  operator= " << c 
		      << "  time= "
		      << swiss.getTimeInSeconds() 
		      << " secs" << endl;
	
	
		pop(xml_out); // BaryonCreationOperator

	
	
	
	
	
	
//	} // end for c (operator within a coeff file)
  //    } // end for f (coeff file)

    //  pop(xml_out); // BaryonCreationOperator

     /* swatch.stop();

      QDPIO::cout << "Creation Operator computed: time= "
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;
*/

      // Annihilation operator
     /* swatch.start();
      BaryonOperator_t  annih_oper;
      annih_oper.mom2_max    = params.param.mom2_max;
      annih_oper.decay_dir   = j_decay;
      annih_oper.seed_l      = quarks[0].seed;
      annih_oper.seed_m      = quarks[1].seed;
      annih_oper.seed_r      = quarks[2].seed;
      annih_oper.smearing    = params.param.sink_quark_smearing;
      annih_oper.perms       = perms;
      annih_oper.orderings.resize(num_orderings);
*/
      push(xml_out, "BaryonAnnihilationOperator");

      // Construct annihilation operator
  //    QDPIO::cout << "Build annihilation operator" << endl;

      // Loop over all files of operators
   //   for(int f=0; f < coeffs.size(); ++f)
     // {
	// Set the id to be used in the analysis codes
//	annih_oper.id = params.named_obj.operator_coeff_files[f].id; 

	// Loop over each operator within a file
//	for(int c=0; c < coeffs[f].ops.size(); ++c)
	//{
	  QDPIO::cout << "Annihilation operator: f=" << f << "  op= " << c << endl;
		 
		annih_oper.id = coeffs[f].ops[c].name;

		write(xml_out, "Name", annih_oper.id);
	  // The operator number. This is just the operator number within the 
	  // coefficient file
	  annih_oper.operator_num = c;

	  // Loop over all orderings and build the operator
	  swiss.reset();
	  swiss.start();

	  for(int ord=0; ord < annih_oper.orderings.size(); ++ord)
	  {
	    QDPIO::cout << "Annihilation operator: ordering = " << ord << endl;
	  
	    annih_oper.orderings[ord].perm = perms[ord];
	    annih_oper.orderings[ord].time_slices.resize(participating_time_slices.size());
     
	    const int n0 = perms[ord][0];
	    const int n1 = perms[ord][1];
	    const int n2 = perms[ord][2];

	    // The operator must hold all the dilutions
	    // We know that all time slices match. However, not all time slices of the
	    // lattice maybe used
	    for(int t=0; t < participating_time_slices.size(); ++t)
	    {
	      BaryonOperator_t::Orderings_t::TimeSlices_t& aop = annih_oper.orderings[ord].time_slices[t];

	      aop.t0 = participating_time_slices[t];
	      aop.dilutions.resize(quarks[n0].time_slices[t].dilutions.size(), 
				   quarks[n1].time_slices[t].dilutions.size(), 
				   quarks[n2].time_slices[t].dilutions.size());

	      for(int i=0; i < quarks[n0].time_slices[t].dilutions.size(); ++i)
	      {
		for(int j=0; j < quarks[n1].time_slices[t].dilutions.size(); ++j)
		{
		  for(int k=0; k < quarks[n2].time_slices[t].dilutions.size(); ++k)
		  {
		    // The correlator
		    LatticeComplex bar = zero;

		    // Loop over the rows/terms within an operator
		    for(int l=0; l < coeffs[f].ops[c].op.size(); ++l)
		    {
		      const OperCoeffs_t::CoeffTerms_t::CoeffTerm_t& term_q = coeffs[f].ops[c].op[l];
	     
					
					KeySmearedDispColorVector_t key0;
		      key0.displacement = term_q.quark[0].displacement;
		      key0.spin         = term_q.quark[0].spin;
	      
		      KeySmearedDispColorVector_t key1;
		      key1.displacement = term_q.quark[1].displacement;
		      key1.spin         = term_q.quark[1].spin;
	      
		      KeySmearedDispColorVector_t key2;
		      key2.displacement = term_q.quark[2].displacement;
		      key2.spin         = term_q.quark[2].spin;
					
					
					/* 
		      KeySmearedDispColorVector_t key0;
		      key0.displacement = term_q.quark[n0].displacement;
		      key0.spin         = term_q.quark[n0].spin;
	      
		      KeySmearedDispColorVector_t key1;
		      key1.displacement = term_q.quark[n1].displacement;
		      key1.spin         = term_q.quark[n1].spin;
	      
		      KeySmearedDispColorVector_t key2;
		      key2.displacement = term_q.quark[n2].displacement;
		      key2.spin         = term_q.quark[n2].spin;
	      */
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

		      // Contract over color indices with antisym tensor
		      LatticeComplex b_oper =
			colorContract(disp_q0.quarks[n0].time_slices[t].dilutions[i].soln,
				      disp_q1.quarks[n1].time_slices[t].dilutions[j].soln,
				      disp_q2.quarks[n2].time_slices[t].dilutions[k].soln);
		      
		      bar += term_q.coeff * b_oper;
		    } // end for l

		    // Slow fourier-transform
		    multi2d<DComplex> hsum(phases.sft(bar));

		    // Unpack into separate momentum and correlator
		    aop.dilutions(i,j,k).mom_projs.resize(phases.numMom());
		    for(int sink_mom_num=0; sink_mom_num < phases.numMom(); ++sink_mom_num) 
		    {
		      aop.dilutions(i,j,k).mom_projs[sink_mom_num].mom  = phases.numToMom(sink_mom_num);
		      aop.dilutions(i,j,k).mom_projs[sink_mom_num].op   = hsum[sink_mom_num];
		    }
		  } // end for k
		} // end for j
	      } // end for i
	    } // end for t
	  } // end for ord

	  swiss.stop();

	  QDPIO::cout << "Annihilation operator construction: file= " << f 
		      << "  operator= " << c 
		      << "  time= "
		      << swiss.getTimeInSeconds() 
		      << " secs" << endl;

					for (int pr = 0 ; pr < 6 ; ++pr )
					{
		QDPIO::cout<<"Sink testval = "<<annih_oper.orderings[pr].time_slices[0].dilutions(0,0,0).mom_projs[0].op[0]<<endl;
					}
	  // Write the meta-data and the binary for this operator
	  swiss.reset();
	  swiss.start();
	  {
	    XMLBufferWriter     record_xml;
	    BinaryBufferWriter  record_bin;

	    write(record_xml, "BaryonAnnihilationOperator", annih_oper);
	    write(record_bin, annih_oper);

	    write(qdp_file, record_xml, record_bin);
	  }
	  swiss.stop();

	  QDPIO::cout << "Annihilation operator writing: file= " << f 
		      << "  operator= " << c 
		      << "  time= "
		      << swiss.getTimeInSeconds() 
		      << " secs" << endl;
	
		
		 pop(xml_out); // OperatorB


		
		
		} // end for c (operator within a coeff file)
      } // end for f (coeff file)

          // swatch.stop();

      /*QDPIO::cout << "Annihilation operator computed: time= "
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;
*/

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
