// $Id: inline_stoch_group_baryon_w.cc,v 1.1 2007-06-18 19:40:03 edwards Exp $
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

      param.source_quark_smearing = readXMLGroup(paramtop, "SourceQuarkSmearing", "wvf_kind");
      param.sink_quark_smearing   = readXMLGroup(paramtop, "SinkQuarkSmearing", "wvf_kind");
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
      read(inputtop, "operator_coeff_files", input.operator_coeff_files);
      read(inputtop, "operator_file", input.operator_file);
      read(inputtop, "Operator", input.op);
    }

    //! Propagator parameters
    void write(XMLWriter& xml, const string& path, const InlineStochGroupBaryonEnv::Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "operator_coeff_files", input.operator_coeff_files);
      write(xml, "operator_file", input.operator_file);
      write(xml, "Operator", input.op);

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
      struct QuarkSolution_t
      {
	LatticeFermion     source;
	LatticeFermion     soln;
	PropSourceConst_t  source_header;
	ChromaProp_t       prop_header;
      };

      int   decay_dir;
      Seed  seed;
      multi1d<QuarkSolution_t>  dilutions;
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
	struct Dilutions_t
	{
	  LatticeColorVector  source;
	  LatticeColorVector  soln;
	};

	multi1d<Dilutions_t> dilutions;
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
	//! Baryon operator dilution
	struct Dilutions_t
	{
	  //! Momentum projected correlator
	  struct Mom_t
	  {
	    multi1d<int>       mom;    /*!< D-1 momentum of this correlator*/
	    multi1d<DComplex>  corr;   /*!< Momentum projected correlator */
	  };

	  multi1d<Mom_t> corrs;        /*!< Holds momentum projected correlators */
	};
	  
	multi1d<int> perm;                     /*!< This particular permutation of quark orderings */
	multi3d<Dilutions_t> dilutions;        /*!< Hybrid list indices */
      };

      multi1d< multi1d<int> > perms;   /*!< Permutations of quark enumeration */

      GroupXML_t    smearing;          /*!< String holding quark smearing xml */

      Seed          seed_l;            /*!< Id of left quark */
      Seed          seed_m;            /*!< Id of middle quark */
      Seed          seed_r;            /*!< Id of right quark */

      int           mom2_max;
      int           decay_dir;
      multi1d<Orderings_t> orderings;  /*!< Array is over quark orderings */
    };


    //! BaryonOperator header writer
    void write(XMLWriter& xml, const string& path, const BaryonOperator_t& param)
    {
      if( path != "." )
	push(xml, path);

      int version = 1;
      write(xml, "version", version);
      write(xml, "mom2_max", param.mom2_max);
      write(xml, "decay_dir", param.decay_dir);
      write(xml, "seed_l", param.seed_l);
      write(xml, "seed_m", param.seed_m);
      write(xml, "seed_r", param.seed_r);
      xml <<  param.smearing.xml;
      write(xml, "perms", param.perms);

      if( path != "." )
	pop(xml);
    }


    //! BaryonOperator binary writer
    void write(BinaryWriter& bin, const BaryonOperator_t::Orderings_t::Dilutions_t::Mom_t& param)
    {
      write(bin, param.mom);
      write(bin, param.corr);
    }

    //! BaryonOperator binary writer
    void write(BinaryWriter& bin, const BaryonOperator_t::Orderings_t::Dilutions_t& param)
    {
      write(bin, param.corrs);
    }

    //! BaryonOperator binary writer
    void write(BinaryWriter& bin, const BaryonOperator_t::Orderings_t& param)
    {
      write(bin, param.perm);
      write(bin, param.dilutions);
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
	reader >> num_terms;
	coeffs.ops[n].op.resize(num_terms);
	
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

		for(int n=0; n < disp_q.quarks.size(); ++n)
		{
		  disp_q.quarks[n].dilutions.resize(quarks[n].dilutions.size());

		  for(int m=0; m < disp_q.quarks[n].dilutions.size(); ++m)
		  {
		    // Short-hand
		    SmearedDispColorVector_t::Quarks_t::Dilutions_t& dil = 
		      disp_q.quarks[n].dilutions[m];

		    // Pull out the appropriate spin component, then displace it
		    dil.source = peekSpin(quarks[n].dilutions[m].source, term_q.spin);
		    dil.soln   = peekSpin(quarks[n].dilutions[m].soln,   term_q.spin);

		    displacement(u_smr, dil.source, term_q.disp_len, term_q.disp_dir);
		    displacement(u_smr, dil.soln,   term_q.disp_len, term_q.disp_dir);
		  } // for m
		} // for n
		
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

      // First calculate some gauge invariant observables just for info.
      // This is really cheap.
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

      multi1d<QuarkSourceSolutions_t>  quarks(params.named_obj.op.size());
      QDPIO::cout << "num_quarks= " << params.named_obj.op.size() << endl;

      try
      {
	QDPIO::cout << "quarks.size= " << quarks.size() << endl;
	for(int n=0; n < quarks.size(); ++n)
	{
	  QDPIO::cout << "Attempt to read solutions for source number=" << n << endl;
	  quarks[n].dilutions.resize(params.named_obj.op[n].soln_files.size());

	  QDPIO::cout << "dilutions.size= " << quarks[n].dilutions.size() << endl;
	  for(int i=0; i < quarks[n].dilutions.size(); ++i)
	  {
	    XMLReader file_xml, record_xml;

	    QDPIO::cout << "reading file= " << params.named_obj.op[n].soln_files[i] << endl;
	    QDPFileReader from(file_xml, params.named_obj.op[n].soln_files[i], QDPIO_SERIAL);
	    read(from, record_xml, quarks[n].dilutions[i].soln);
	    close(from);
	
	    read(record_xml, "/Propagator/PropSource", quarks[n].dilutions[i].source_header);
	    read(record_xml, "/Propagator/ForwardProp", quarks[n].dilutions[i].prop_header);
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
      // Check for each quark source that the solutions have their diluted
      // on every site only once
      //
      swatch.start();

      try
      {
	push(xml_out, "Norms");
	for(int n=0; n < quarks.size(); ++n)
	{
	  bool first = true;
	  int  N;
	  LatticeFermion quark_noise;      // noisy source on entire lattice

	  for(int i=0; i < quarks[n].dilutions.size(); ++i)
	  {
	    std::istringstream  xml_s(quarks[n].dilutions[i].source_header.source.xml);
	    XMLReader  sourcetop(xml_s);
//	QDPIO::cout << "Source = " << quarks[n].dilutions[i].source_header.source.id << endl;

	    if (quarks[n].dilutions[i].source_header.source.id != DiluteZNQuarkSourceConstEnv::name)
	    {
	      QDPIO::cerr << "Expected source_type = " << DiluteZNQuarkSourceConstEnv::name << endl;
	      QDP_abort(1);
	    }

	    QDPIO::cout << "Quark num= " << n << "  dilution num= " << i << endl;

	    // Manually create the params so I can peek into them and use the source constructor
	    DiluteZNQuarkSourceConstEnv::Params  srcParams(sourcetop, 
							   quarks[n].dilutions[i].source_header.source.path);
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
	    quarks[n].dilutions[i].source = srcConst(u);
	    quark_noise -= quarks[n].dilutions[i].source;

#if 0
	    // Diagnostic
	    {
	      // Keep a copy of the phases with NO momenta
	      SftMom phases_nomom(0, true, quarks[n].dilutions[i].source_header.j_decay);

	      multi1d<Double> source_corr = sumMulti(localNorm2(quarks[n].dilutions[i].source), 
						     phases_nomom.getSet());

	      multi1d<Double> soln_corr = sumMulti(localNorm2(quarks[n].dilutions[i].soln), 
						   phases_nomom.getSet());

	      push(xml_out, "elem");
	      write(xml_out, "n", n);
	      write(xml_out, "i", i);
	      write(xml_out, "source_corr", source_corr);
	      write(xml_out, "soln_corr", soln_corr);
	      pop(xml_out);
	    }
#endif
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

      QDPFileWriter qdp_file(file_xml, params.named_obj.operator_file,     // are there one or two files???
			     QDPIO_SINGLEFILE, QDPIO_SERIAL, QDPIO_OPEN);

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
	const std::string& file = params.named_obj.operator_coeff_files[i];
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
	for(int n=0; n < quarks.size(); ++n)
	{
	  QDPIO::cout << "Smearing sources for quark[" << n << "]  over dilutions = " 
		      << quarks[n].dilutions.size() << endl;

	  for(int i=0; i < quarks[n].dilutions.size(); ++i)
	  {
	    LatticeFermion src(quarks[n].dilutions[i].source);
	    (*sourceQuarkSmearing)(src, u);

	    quarks[n].dilutions[i].source = rotate_mat * src;
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
		      << quarks[n].dilutions.size() << endl;

	  for(int i=0; i < quarks[n].dilutions.size(); ++i)
	  {
	    LatticeFermion soln(quarks[n].dilutions[i].soln);
	    (*sinkQuarkSmearing)(soln, u);

	    quarks[n].dilutions[i].soln = rotate_mat * soln;
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
      int j_decay = quarks[0].decay_dir;

      // Initialize the slow Fourier transform phases
      SftMom phases(params.param.mom2_max, false, j_decay);
    
      // Length of lattice in decay direction
      int length = phases.numSubsets();

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
      swatch.start();
      BaryonOperator_t  baryon_opA;
      baryon_opA.mom2_max    = params.param.mom2_max;
      baryon_opA.decay_dir   = j_decay;
      baryon_opA.seed_l      = quarks[0].seed;
      baryon_opA.seed_m      = quarks[1].seed;
      baryon_opA.seed_r      = quarks[2].seed;
      baryon_opA.smearing    = params.param.source_quark_smearing;
      baryon_opA.orderings.resize(num_orderings);
      baryon_opA.perms.resize(num_orderings);

      push(xml_out, "OperatorA");

      // Sanity check
      if ( toBool(baryon_opA.seed_l == baryon_opA.seed_m) )
      {
	QDPIO::cerr << "baryon op seeds are the same" << endl;
	QDP_abort(1);
      }

      // Sanity check
      if ( toBool(baryon_opA.seed_l == baryon_opA.seed_r) )
      {
	QDPIO::cerr << "baryon op seeds are the same" << endl;
	QDP_abort(1);
      }

      // Sanity check
      if ( toBool(baryon_opA.seed_m == baryon_opA.seed_r) )
      {
	QDPIO::cerr << "baryon op seeds are the same" << endl;
	QDP_abort(1);
      }


      // Construct operator A
      QDPIO::cout << "Build creation operator" << endl;

      // Loop over all files of operators
      for(int f=0; f < coeffs.size(); ++f)
      {
	// Loop over each operator within a file
	for(int c=0; c < coeffs[f].ops.size(); ++c)
	{
	  QDPIO::cout << "Creation operator: f=" << f << "  op= " << c << endl;

	  // Loop over all orderings and build the operator
	  swiss.reset();
	  swiss.start();

	  for(int ord=0; ord < baryon_opA.orderings.size(); ++ord)
	  {
	    QDPIO::cout << "Creation operator: ordering = " << ord << endl;
	  
	    baryon_opA.perms[ord] = perms[ord];
     
	    const int n0 = perms[ord][0];
	    const int n1 = perms[ord][1];
	    const int n2 = perms[ord][2];

	    // The operator must hold all the dilutions
	    BaryonOperator_t::Orderings_t& ordering = baryon_opA.orderings[ord];
	    ordering.dilutions.resize(quarks[n0].dilutions.size(), 
				      quarks[n1].dilutions.size(), 
				      quarks[n2].dilutions.size());

	    for(int i=0; i < quarks[n0].dilutions.size(); ++i)
	    {
	      for(int j=0; j < quarks[n1].dilutions.size(); ++j)
	      {
		for(int k=0; k < quarks[n2].dilutions.size(); ++k)
		{
		  // The correlator
		  LatticeComplex bar = zero;

		  // Loop over the rows/terms within an operator
		  for(int l=0; l < coeffs[f].ops[c].op.size(); ++l)
		  {
		    const OperCoeffs_t::CoeffTerms_t::CoeffTerm_t& term_q = coeffs[f].ops[c].op[l];
	      
		    KeySmearedDispColorVector_t key0;
		    key0.displacement = term_q.quark[n0].displacement;
		    key0.spin         = term_q.quark[n0].spin;
	      
		    KeySmearedDispColorVector_t key1;
		    key1.displacement = term_q.quark[n1].displacement;
		    key1.spin         = term_q.quark[n1].spin;
	      
		    KeySmearedDispColorVector_t key2;
		    key2.displacement = term_q.quark[n2].displacement;
		    key2.spin         = term_q.quark[n2].spin;
	      
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
		      colorContract(disp_q0.quarks[n0].dilutions[i].source,
				    disp_q1.quarks[n1].dilutions[j].source,
				    disp_q2.quarks[n2].dilutions[k].source);
		    
		    bar += term_q.coeff * b_oper;
		  } // end for l

		  // Slow fourier-transform
		  multi2d<DComplex> hsum(phases.sft(bar));

		  // Unpack into separate momentum and correlator
		  ordering.dilutions(i,j,k).corrs.resize(phases.numMom());
		  for(int sink_mom_num=0; sink_mom_num < phases.numMom(); ++sink_mom_num) 
		  {
		    ordering.dilutions(i,j,k).corrs[sink_mom_num].mom  = phases.numToMom(sink_mom_num);
		    ordering.dilutions(i,j,k).corrs[sink_mom_num].corr = hsum[sink_mom_num];
		  }
		} // end for k
	      } // end for j
	    } // end for i
	  } // end for ord

	  swiss.stop();

	  QDPIO::cout << "Creation operator construction: file= " << f 
		      << "  operator= " << c 
		      << "  time= "
		      << swiss.getTimeInSeconds() 
		      << " secs" << endl;

	  // Write the meta-data and the binary for this operator
	  swiss.start();
	  XMLBufferWriter     record_xml;
	  BinaryBufferWriter  record_bin;

	  write(record_xml, "BaryonCreationOperator", baryon_opA);
	  write(record_bin, baryon_opA);

	  write(qdp_file, record_xml, record_bin);
	  swiss.stop();

	  QDPIO::cout << "Creation operator writing: file= " << f 
		      << "  operator= " << c 
		      << "  time= "
		      << swiss.getTimeInSeconds() 
		      << " secs" << endl;
	} // end for c (operator within a coeff file)
      } // end for f (coeff file)

      pop(xml_out); // OperatorA

      swatch.stop();

      QDPIO::cout << "Creation Operator computed: time= "
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;


      // Operator B
      swatch.start();
      BaryonOperator_t  baryon_opB;
      baryon_opB.mom2_max    = params.param.mom2_max;
      baryon_opB.decay_dir   = j_decay;
      baryon_opB.seed_l      = quarks[0].seed;
      baryon_opB.seed_m      = quarks[1].seed;
      baryon_opB.seed_r      = quarks[2].seed;
      baryon_opB.smearing    = params.param.sink_quark_smearing;
      baryon_opB.orderings.resize(num_orderings);
      baryon_opB.perms.resize(num_orderings);

      push(xml_out, "OperatorB");

      // Sanity check
      if ( toBool(baryon_opB.seed_l == baryon_opB.seed_m) )
      {
	QDPIO::cerr << "baryon op seeds are the same" << endl;
	QDP_abort(1);
      }

      // Sanity check
      if ( toBool(baryon_opB.seed_l == baryon_opB.seed_r) )
      {
	QDPIO::cerr << "baryon op seeds are the same" << endl;
	QDP_abort(1);
      }

      // Sanity check
      if ( toBool(baryon_opB.seed_m == baryon_opB.seed_r) )
      {
	QDPIO::cerr << "baryon op seeds are the same" << endl;
	QDP_abort(1);
      }


      // Construct operator B
      QDPIO::cout << "Build annihilation operator" << endl;

      // Loop over all files of operators
      for(int f=0; f < coeffs.size(); ++f)
      {
	// Loop over each operator within a file
	for(int c=0; c < coeffs[f].ops.size(); ++c)
	{
	  QDPIO::cout << "Annihilation operator: f=" << f << "  op= " << c << endl;

	  // Loop over all orderings and build the operator
	  swiss.reset();
	  swiss.start();

	  for(int ord=0; ord < baryon_opB.orderings.size(); ++ord)
	  {
	    QDPIO::cout << "Annihilation operator: ordering = " << ord << endl;
	  
	    baryon_opB.perms[ord] = perms[ord];
     
	    const int n0 = perms[ord][0];
	    const int n1 = perms[ord][1];
	    const int n2 = perms[ord][2];

	    // The operator must hold all the dilutions
	    BaryonOperator_t::Orderings_t& ordering = baryon_opB.orderings[ord];
	    ordering.dilutions.resize(quarks[n0].dilutions.size(), 
				      quarks[n1].dilutions.size(), 
				      quarks[n2].dilutions.size());

	    for(int i=0; i < quarks[n0].dilutions.size(); ++i)
	    {
	      for(int j=0; j < quarks[n1].dilutions.size(); ++j)
	      {
		for(int k=0; k < quarks[n2].dilutions.size(); ++k)
		{
		  // The correlator
		  LatticeComplex bar = zero;

		  // Loop over the rows/terms within an operator
		  for(int l=0; l < coeffs[f].ops[c].op.size(); ++l)
		  {
		    const OperCoeffs_t::CoeffTerms_t::CoeffTerm_t& term_q = coeffs[f].ops[c].op[l];
	      
		    KeySmearedDispColorVector_t key0;
		    key0.displacement = term_q.quark[n0].displacement;
		    key0.spin         = term_q.quark[n0].spin;
	      
		    KeySmearedDispColorVector_t key1;
		    key1.displacement = term_q.quark[n1].displacement;
		    key1.spin         = term_q.quark[n1].spin;
	      
		    KeySmearedDispColorVector_t key2;
		    key2.displacement = term_q.quark[n2].displacement;
		    key2.spin         = term_q.quark[n2].spin;
	      
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
		      colorContract(disp_q0.quarks[n0].dilutions[i].source,
				    disp_q1.quarks[n1].dilutions[j].source,
				    disp_q2.quarks[n2].dilutions[k].source);
		    
		    bar += term_q.coeff * b_oper;
		  } // end for l

		  // Slow fourier-transform
		  multi2d<DComplex> hsum(phases.sft(bar));

		  // Unpack into separate momentum and correlator
		  ordering.dilutions(i,j,k).corrs.resize(phases.numMom());
		  for(int sink_mom_num=0; sink_mom_num < phases.numMom(); ++sink_mom_num) 
		  {
		    ordering.dilutions(i,j,k).corrs[sink_mom_num].mom  = phases.numToMom(sink_mom_num);
		    ordering.dilutions(i,j,k).corrs[sink_mom_num].corr = hsum[sink_mom_num];
		  }
		} // end for k
	      } // end for j
	    } // end for i
	  } // end for ord

	  swiss.stop();

	  QDPIO::cout << "Annihilation operator construction: file= " << f 
		      << "  operator= " << c 
		      << "  time= "
		      << swiss.getTimeInSeconds() 
		      << " secs" << endl;

	  // Write the meta-data and the binary for this operator
	  swiss.start();
	  XMLBufferWriter     record_xml;
	  BinaryBufferWriter  record_bin;

	  write(record_xml, "BaryonAnnihilationOperator", baryon_opB);
	  write(record_bin, baryon_opB);

	  write(qdp_file, record_xml, record_bin);
	  swiss.stop();

	  QDPIO::cout << "Annihilation operator writing: file= " << f 
		      << "  operator= " << c 
		      << "  time= "
		      << swiss.getTimeInSeconds() 
		      << " secs" << endl;
	} // end for c (operator within a coeff file)
      } // end for f (coeff file)

      pop(xml_out); // OperatorB

      swatch.stop();

      QDPIO::cout << "Annihilation operator computed: time= "
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;


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
