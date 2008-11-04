// $Id: inline_stoch_baryon_w.cc,v 3.14 2008-11-04 18:43:57 edwards Exp $
/*! \file
 * \brief Inline measurement of stochastic baryon operator
 *
 */

#include "handle.h"
#include "meas/inline/hadron/inline_stoch_baryon_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/sources/dilutezN_source_const.h"
#include "meas/sources/zN_src.h"
#include "meas/hadron/barspinmat_w.h"
#include "meas/hadron/baryon_operator_aggregate_w.h"
#include "meas/hadron/baryon_operator_factory_w.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlineStochBaryonEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineStochBaryon(InlineStochBaryonParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "STOCH_BARYON";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= BaryonOperatorEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
  }



  // Operator parameters
  void read(XMLReader& xml, const string& path, InlineStochBaryonParams::Prop_t::Operator_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "soln_files", input.soln_files);
  }


  // Operator parameters
  void write(XMLWriter& xml, const string& path, const InlineStochBaryonParams::Prop_t::Operator_t& input)
  {
    push(xml, path);
    write(xml, "soln_files", input.soln_files);
    pop(xml);
  }


  // Propagator parameters
  void read(XMLReader& xml, const string& path, InlineStochBaryonParams::Prop_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "operator_file", input.op_file);
    read(inputtop, "operator", input.op);
  }


  // Propagator parameters
  void write(XMLWriter& xml, const string& path, const InlineStochBaryonParams::Prop_t& input)
  {
    push(xml, path);

    write(xml, "operator_file", input.op_file);
    write(xml, "operator", input.op);

    pop(xml);
  }


  // Reader for input parameters
  void read(XMLReader& xml, const string& path, InlineStochBaryonParams::Param_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    switch (version) 
    {
    case 2:
      /**************************************************************************/
      break;

    default :
      /**************************************************************************/

      QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
      QDP_abort(1);
    }

    read(paramtop, "mom2_max", param.mom2_max);
    
    {
      XMLReader xml_tmp(paramtop, "BaryonOperator");
      std::ostringstream os;
      xml_tmp.print(os);
      read(xml_tmp, "BaryonOperatorType", param.baryon_operator_type);
      param.baryon_operator = os.str();
    }

  }


  // Reader for input parameters
  void write(XMLWriter& xml, const string& path, const InlineStochBaryonParams::Param_t& param)
  {
    push(xml, path);

    int version = 1;

    write(xml, "version", version);
    write(xml, "mom2_max", param.mom2_max);

    pop(xml);
  }


  //! Propagator parameters
  void read(XMLReader& xml, const string& path, InlineStochBaryonParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "Prop", input.prop);
  }

  //! Propagator parameters
  void write(XMLWriter& xml, const string& path, const InlineStochBaryonParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "Prop", input.prop);

    pop(xml);
  }


  // Param stuff
  InlineStochBaryonParams::InlineStochBaryonParams()
  { 
    frequency = 0; 
  }

  InlineStochBaryonParams::InlineStochBaryonParams(XMLReader& xml_in, const std::string& path) 
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
  InlineStochBaryonParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    // Parameters for source construction
    Chroma::write(xml_out, "Param", param);

    // Write out the output propagator/source configuration info
    Chroma::write(xml_out, "NamedObject", named_obj);

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

    int   j_decay;
    Seed  seed;
    multi1d<QuarkSolution_t>  dilutions;
  };


  //! Baryon operator
  struct BaryonOperator_t
  {

    //! Baryon operator
    struct BaryonOperatorInsertion_t
    {
      //! Possible operator index
      struct BaryonOperatorIndex_t
      {
	//! Baryon operator element
	struct BaryonOperatorElement_t
	{
	  multi2d<DComplex> elem;              /*!< time slice and momenta number */
	};

	multi1d<BaryonOperatorElement_t> ind;
      };
    
      multi3d<BaryonOperatorIndex_t> op;    /*!< hybrid list indices */

    };

    multi1d< multi1d<int> > perms;   /*!< Permutations of quark enumeration */

    Seed          seed_l;            /*!< Id of left quark */
    Seed          seed_m;            /*!< Id of middle quark */
    Seed          seed_r;            /*!< Id of right quark */

    int           mom2_max;
    int           j_decay;
    multi1d<BaryonOperatorInsertion_t> orderings;  /*!< Array is over quark orderings */
  };


  //! BaryonOperator binary writer
  void write(BinaryWriter& bin, 
	     const BaryonOperator_t::BaryonOperatorInsertion_t::BaryonOperatorIndex_t::BaryonOperatorElement_t& p){
    write(bin, p.elem);
  }

  //! BaryonOperator binary writer
  void write(BinaryWriter& bin, 
	     const BaryonOperator_t::BaryonOperatorInsertion_t::BaryonOperatorIndex_t& p){
    write(bin, p.ind);
  }

  //! BaryonOperator binary writer
  void write(BinaryWriter& bin, 
	     const BaryonOperator_t::BaryonOperatorInsertion_t& p){
    write(bin, p.op);
  }

  //! BaryonOperator binary writer
  void write(BinaryWriter& bin, const BaryonOperator_t& param){
    write(bin, param.seed_l);
    write(bin, param.seed_m);
    write(bin, param.seed_r);
    write(bin, param.mom2_max);
    write(bin, param.j_decay);
    write(bin, param.perms);
    write(bin, param.orderings);
  }
  
  //! BaryonOperator header writer
  void write(XMLWriter& xml, const string& path, const BaryonOperator_t& param)
  {
    if( path != "." )
      push(xml, path);

    int version = 1;
    write(xml, "version", version);
    write(xml, "mom2_max", param.mom2_max);
    write(xml, "j_decay", param.j_decay);
    write(xml, "seed_l", param.seed_l);
    write(xml, "seed_m", param.seed_m);
    write(xml, "seed_r", param.seed_r);
    write(xml, "perms", param.perms);

    if( path != "." )
      pop(xml);
  }



  //--------------------------------------------------------------
  // Function call
  void 
  InlineStochBaryon::operator()(unsigned long update_no,
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
  InlineStochBaryon::func(unsigned long update_no,
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
      QDPIO::cerr << InlineStochBaryonEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineStochBaryonEnv::name << ": map call failed: " << e 
		  << endl;
      QDP_abort(1);
    }
    const multi1d<LatticeColorMatrix>& u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

    push(xml_out, "stoch_baryon");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlineStochBaryonEnv::name << ": Stochastic Baryon Operator" << endl;

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.write(xml_out, "Input");

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
    StopWatch swatch;
    swatch.reset();
    swatch.start();

    multi1d<QuarkSourceSolutions_t>  quarks(params.named_obj.prop.op.size());
    QDPIO::cout << "num_quarks= " << params.named_obj.prop.op.size() << endl;

    try
    {
      QDPIO::cout << "quarks.size= " << quarks.size() << endl;
      for(int n=0; n < quarks.size(); ++n)
      {
	QDPIO::cout << "Attempt to read solutions for source number=" << n << endl;
	quarks[n].dilutions.resize(params.named_obj.prop.op[n].soln_files.size());

	QDPIO::cout << "dilutions.size= " << quarks[n].dilutions.size() << endl;
	for(int i=0; i < quarks[n].dilutions.size(); ++i)
	{
	  XMLReader file_xml, record_xml;

	  QDPIO::cout << "reading file= " << params.named_obj.prop.op[n].soln_files[i] << endl;
	  QDPFileReader from(file_xml, params.named_obj.prop.op[n].soln_files[i], QDPIO_SERIAL);
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

	  if (quarks[n].dilutions[i].source_header.source.id != DiluteZNQuarkSourceConstEnv::getName())
	  {
	    QDPIO::cerr << "Expected source_type = " << DiluteZNQuarkSourceConstEnv::getName() << endl;
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

	    quarks[0].j_decay = srcParams.j_decay;

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
    // Baryon operators
    //
    int j_decay = quarks[0].j_decay;

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
    // Create the baryon operator object
    //
    std::istringstream  xml_op(params.param.baryon_operator);
    XMLReader  optop(xml_op);
    const string operator_path = "/BaryonOperator";
	
    Handle< BaryonOperator<LatticeFermion> >
      baryonOperator(TheWilsonBaryonOperatorFactory::Instance().createObject(params.param.baryon_operator_type,
									     optop,
									     operator_path,
									     u));

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

    // Operator A
    swatch.start();
    BaryonOperator_t  baryon_opA;
    baryon_opA.mom2_max    = params.param.mom2_max;
    baryon_opA.j_decay     = j_decay;
    baryon_opA.seed_l      = quarks[0].seed;
    baryon_opA.seed_m      = quarks[1].seed;
    baryon_opA.seed_r      = quarks[2].seed;
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
    try
    {
      for(int ord=0; ord < baryon_opA.orderings.size(); ++ord)
      {
	QDPIO::cout << "Operator A: ordering = " << ord << endl;

	baryon_opA.perms[ord] = perms[ord];
      
	// Operator construction
	const QuarkSourceSolutions_t& q0 = quarks[perms[ord][0]];
	const QuarkSourceSolutions_t& q1 = quarks[perms[ord][1]];
	const QuarkSourceSolutions_t& q2 = quarks[perms[ord][2]];

	baryon_opA.orderings[ord].op.resize(q0.dilutions.size(),
					    q1.dilutions.size(),
					    q2.dilutions.size());

	
	for(int i=0; i < q0.dilutions.size(); ++i)
	{
	  for(int j=0; j < q1.dilutions.size(); ++j)
	  {
	    for(int k=0; k < q2.dilutions.size(); ++k)
	    {
	      multi1d<LatticeComplex> bar = (*baryonOperator)(q0.dilutions[i].source,
							      q1.dilutions[j].source,
							      q2.dilutions[k].source,
							      MINUS);

	      baryon_opA.orderings[ord].op(i,j,k).ind.resize(bar.size());
	      for(int l=0; l < bar.size(); ++l)
		baryon_opA.orderings[ord].op(i,j,k).ind[l].elem = phases.sft(bar[l]);
	      
	    } // end for k
	  } // end for j
	} // end for i
      } // end for ord
    } // end try
    catch(const std::string& e) 
    {
      QDPIO::cerr << ": Caught Exception creating source operator: " << e << endl;
      QDP_abort(1);
    }
    catch(...)
    {
      QDPIO::cerr << ": Caught generic exception creating source operator" << endl;
      QDP_abort(1);
    }

    pop(xml_out); // OperatorA

    swatch.stop();

    QDPIO::cout << "Operator A computed: time= "
		<< swatch.getTimeInSeconds() 
		<< " secs" << endl;


    // Operator B
    swatch.start();
    BaryonOperator_t  baryon_opB;
    baryon_opB.mom2_max    = params.param.mom2_max;
    baryon_opB.j_decay     = j_decay;
    baryon_opB.seed_l      = quarks[0].seed;
    baryon_opB.seed_m      = quarks[1].seed;
    baryon_opB.seed_r      = quarks[2].seed;
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
    try
    {
      for(int ord=0; ord < baryon_opB.orderings.size(); ++ord)
      {
	QDPIO::cout << "Operator B: ordering = " << ord << endl;

	baryon_opB.perms[ord] = perms[ord];
      
	// Operator construction
	const QuarkSourceSolutions_t& q0 = quarks[perms[ord][0]];
	const QuarkSourceSolutions_t& q1 = quarks[perms[ord][1]];
	const QuarkSourceSolutions_t& q2 = quarks[perms[ord][2]];

	baryon_opB.orderings[ord].op.resize(q0.dilutions.size(),
					    q1.dilutions.size(),
					    q2.dilutions.size());
	
	for(int i=0; i < q0.dilutions.size(); ++i)
	{
	  for(int j=0; j < q1.dilutions.size(); ++j)
	  {
	    for(int k=0; k < q2.dilutions.size(); ++k)
	    {
	      multi1d<LatticeComplex> bar = (*baryonOperator)(q0.dilutions[i].soln,
							      q1.dilutions[j].soln,
							      q2.dilutions[k].soln,
							      PLUS);

	      baryon_opB.orderings[ord].op(i,j,k).ind.resize(bar.size());
	      for(int l=0; l < bar.size(); ++l)
		baryon_opB.orderings[ord].op(i,j,k).ind[l].elem = phases.sft(bar[l]);

	    } // end for k
	  } // end for j
	} // end for i
      } // end for ord
    } // end try
    catch(const std::string& e) 
    {
      QDPIO::cerr << ": Caught Exception creating sink operator: " << e << endl;
      QDP_abort(1);
    }
    catch(...)
    {
      QDPIO::cerr << ": Caught generic exception creating sink operator" << endl;
      QDP_abort(1);
    }

    pop(xml_out); // OperatorB

    swatch.stop();

    QDPIO::cout << "Operator B computed: time= "
		<< swatch.getTimeInSeconds() 
		<< " secs" << endl;


    // Save the operators
    // ONLY SciDAC output format is supported!
    swatch.start();
    {
      XMLBufferWriter file_xml;
      push(file_xml, "baryon_operator");
      file_xml << params.param.baryon_operator;
      write(file_xml, "Config_info", gauge_xml);
      pop(file_xml);

      QDPFileWriter to(file_xml, params.named_obj.prop.op_file, 
		       QDPIO_SINGLEFILE, QDPIO_SERIAL, QDPIO_OPEN);

      // Write the scalar data
      {
	XMLBufferWriter record_xml;
	BinaryBufferWriter  record_bin;
	write(record_xml, "SourceBaryonOperator", baryon_opA);
	write(record_bin, baryon_opA);
	write(to, record_xml, record_bin);
      }

      // Write the scalar data
      {
	XMLBufferWriter record_xml;
	BinaryBufferWriter  record_bin;
	write(record_xml, "SinkBaryonOperator", baryon_opB);
	write(record_bin, baryon_opB);
	write(to, record_xml, record_bin);
      }

      close(to);
    }

    swatch.stop();

    QDPIO::cout << "Operators written: time= "
		<< swatch.getTimeInSeconds() 
		<< " secs" << endl;

    // Close the namelist output file XMLDAT
    pop(xml_out);     // StochBaryon

    snoop.stop();
    QDPIO::cout << InlineStochBaryonEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << InlineStochBaryonEnv::name << ": ran successfully" << endl;

    END_CODE();
  } 

}
