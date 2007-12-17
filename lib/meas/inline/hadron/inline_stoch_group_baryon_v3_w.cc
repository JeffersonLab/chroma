// $Id: inline_stoch_group_baryon_v3_w.cc,v 1.3 2007-12-17 22:44:03 edwards Exp $
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
#include "meas/hadron/dilution_operator_aggregate.h"
#include "meas/hadron/dilution_operator_factory.h"
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
    //! Number of quarks to be used in this construction
    const int N_quarks = 3;


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
      param.dilutions             = readXMLArrayGroup(paramtop, "Dilutions", "DilutionType");
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

      for(int i=0; i < param.correlators.size(); ++i)
	xml << input.dilutions[i].xml;

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


    //! Read named objects 
    void read(XMLReader& xml, const string& path, InlineStochGroupBaryonEnv::Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "operators_file", input.operators_file);
      read(inputtop, "Quark_ids", input.quark_ids);
    }

    //! Write named objects
    void write(XMLWriter& xml, const string& path, const InlineStochGroupBaryonEnv::Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "operators_file", input.operators_file);
      write(xml, "Quark_ids", input.quark_ids);

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
	success &= DilutionOperatorEnv::registerAll();
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



    //--------------------------------------------------------------
    //! 3-quark operator structure
    struct ThreeQuarkOps_t
    {
      struct ThreeQuarkOp_t
      {
	struct QuarkInfo_t
	{
	  int  displacement;    /*!< Orig plus/minus 1-based directional displacements */
	  int  spin;            /*!< 1-based spin index */
	};

	multi1d<QuarkInfo_t>  quark;    /*!< Displacement and spin for each quark */
	std::string name;               /*!< Name of the 3-quark operator */
      };
	
      multi1d<ThreeQuarkOp_t> ops; /*!< 3-quark ops within a file */
    };


    //----------------------------------------------------------------------------
    // Structure holding structures
    //! The key for smeared and displaced color vectors
    struct KeySmearedDispColorVector_t
    {
      int  t0;              /*!< Time of source */
      int  displacement;    /*!< Orig plus/minus 1-based directional displacements */
      int  spin;            /*!< 1-based spin index */
    };


    //! Support for the keys of smeared and displaced color vectors
    bool operator<(const KeySmearedDispColorVector_t& a, const KeySmearedDispColorVector_t& b)
    {
      multi1d<int> lga(3);
      lga[0] = a.displacement;
      lga[1] = a.spin;
      lga[2] = a.t0;

      multi1d<int> lgb(3);
      lgb[0] = b.displacement;
      lgb[1] = b.spin;
      lgb[2] = b.t0;

      return (lga < lgb);
    }


    //----------------------------------------------------------------------------
    //! The smeared and displaced objects
    class SmearedDispObjects
    {
    public:
      //! Constructor from a quark smearing
      SmearedDispObjects(int disp_length,
			 multi1d< Handle< DilutionOperator<LatticeFermion> > > dils,
			 Handle< QuarkSmearing<LatticeFermion> > src, 
			 Handle< QuarkSmearing<LatticeFermion> > snk) :
	displacement_length(disp_length), dilution_operators(dils), 
	sourceQuarkSmearing(src), sinkQuarkSmearing(snk) 
	{
	  srce_vecs.resize(N_quarks);
	  soln_vecs.resize(N_quarks);
	}

      //! Destructor
      ~SmearedDispObjects() {}

      //! Accessor
      virtual LatticeColorVector getDilutedSource(int quark_num, 
						  int t0,
						  DilutionOperator<LatticeFermion>::const_iterator dil, 
						  const KeySmearedDispColorVector_t& key);

      //! Accessor
      virtual LatticeColorVector getDilutedSolution(int quark_num, 
						    int t0,
						    DilutionOperator<LatticeFermion>::const_iterator dil, 
						    const KeySmearedDispColorVector_t& key);

    protected:
      //! Displace an object
      virtual const LatticeColorVector&
      displaceObject(map<KeySmearedDispColorVector_t, SmearedDispColorVector_t>& disp_quark,
		     const KeySmearedDispColorVector_t& key,
		     const multi1d<LatticeColorMatrix>& u_smr,
		     const LatticeFermion& ferm)
	
    private:
      //! Displacement length
      int displacement_length;

      //! Dilutions
      multi1d< Handle< DilutionOperator > > dilution_operators;

      //! Quark smearing
      Handle< QuarkSmearing<LatticeFermion> > sourceQuarkSmearing;

      //! Quark smearing
      Handle< QuarkSmearing<LatticeFermion> > sinkQuarkSmearing;

      //! The value of the map
      struct SmearedDispColorVector_t
      {
	LatticeColorVector  vec;
      }
      multi1d< map<KeySmearedDispColorVector_t, SmearedDispColorVector_t> > srce_vecs;
      multi1d< map<KeySmearedDispColorVector_t, SmearedDispColorVector_t> > soln_vecs;
    };


    //! Accessor
    LatticeColorVector 
    SmearedDispObjects::getDilutedSource(int quark_num, 
					 int t0,
					 DilutionOperator<LatticeFermion>::const_iterator dil, 
					 const KeySmearedDispColorVector_t& key,
					 const multi1d<LatticeColorMatrix>& u_smr)
    {
      if (! dilution_operators[n]->hasTimeSupport(dil,t0))
      {
	QDPIO::cerr << __func__ << ": no time support" << endl;
	QDP_abort(1);
      }
		
      // NOTE: The smearing should be in another map as well
      LatticeFermion src(dilution_operators[quark_num]->dilutedSource(dil_ptr,t0));
      (*sourceQuarkSmearing)(src, u_smr);

      //multiply by gamma_4 as well here
      SpinMatrix mat = rotate_mat * Gamma(8);
	    
      return displaceObject(srce_vecs[quark_num],key,u_smr,LatticeFermion(mat * src));
    }

    //! Accessor
    LatticeColorVector 
    SmearedDispObjects::getDilutedSolution(int quark_num, 
					   int t0,
					   DilutionOperator<LatticeFermion>::const_iterator dil, 
					   const KeySmearedDispColorVector_t& key,
					   const multi1d<LatticeColorMatrix>& u_smr)
    {
      // NOTE: The smearing should be in another map as well
      LatticeFermion soln(dilutio_operators[quark_num]->dilutedSolution(dil,t0));
      (*sinkQuarkSmearing)(soln, u_smr);

      return displaceObject(soln_vecs[quark_num],key,u_smr,LatticeFermion(rotate_mat * soln));
    }


    //! Accessor
    const LatticeColorVector&
    SmearedDispObjects::displaceObject()(map<KeySmearedDispColorVector_t, SmearedDispColorVector_t>& disp_quark,
					 const KeySmearedDispColorVector_t& key,
					 const multi1d<LatticeColorMatrix>& u_smr,
					 const LatticeFermion& ferm)
    {
      // If no entry, then create a displaced version of the quark
      if (disp_quark.find(key) == disp_quark.end())
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
	  disp_quark.insert(std::make_pair(key, disp_empty));

	  // Sanity check - the entry better be there
	  if (disp_quark.find(key) == disp_quark.end())
	  {
	    QDPIO::cerr << __func__ 
			<< ": internal error - could not insert empty key in map"
			<< endl;
	    QDP_abort(1);
	  }		      
	}

	// Modify the previous empty entry
	SmearedDispColorVector_t& disp_q = disp_quark.find(key)->second;
	      
	snoop.reset();
	snoop.start();
	
	// Pull out the appropriate spin component, then displace it
	disp_q.vec = peekSpin(ferm, key.spin-1);

	if (key.displacement > 0)
	{
	  int disp_dir = key.displacement - 1;
	  int disp_len = displacement_length;
	  displacement(u_smr, vec, disp_len, disp_dir);
	}
	else if (key.displacement < 0)
	{
	  int disp_dir = -key.displacement - 1;
	  int disp_len = -displacement_length;
	  displacement(u_smr, vec, disp_len, disp_dir);
	}
	      
	snoop.stop();

	QDPIO::cout << "Displaced Quarks: Spin = "<<key.spin<<" Disp = "
		    << key.displacement <<" Time = "<<snoop.getTimeInSeconds() <<" sec"<<endl;
	// Insert
	disp_quark.insert(std::make_pair(key, disp_q));
      } // if find in map

      // The key now must exist in the map, so return the vector
      SmearedDispColorVector_t& disp_q = disp_quark.find(key)->second;

      return disp_q.vec;
    }


    //----------------------------------------------------------------------------
    //! Baryon operator
    struct BaryonOperator_t
    {
      //! Baryon operator time slices corresponding to location of operator source
      struct TimeSlices_t
      {
	//! Quark orderings within a baryon operator
	struct Orderings_t
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
	  std::list<Dilutions_t> dilutions; /*!< Hybrid list indices */
	};
	  
	multi1d<int> perm;                  /*!< This particular permutation of quark orderings */
	multi1d<TimeSlices_t> time_sources; /*!< Time slices of the lattice that are used */
      };

      multi1d< multi1d<int> > perms;   /*!< Permutations of quark enumeration */

      GroupXML_t    smearing;          /*!< String holding quark smearing xml */

      Seed          seed_l;            /*!< Id of left quark */
      Seed          seed_m;            /*!< Id of middle quark */
      Seed          seed_r;            /*!< Id of right quark */

      std::string   id;                /*!< Tag/ID used in analysis codes */

      int           mom2_max;          /*!< |\vec{p}|^2 */
      int           decay_dir;         /*!< Direction of decay */
      multi1d<Orderings_t> orderings;  /*!< Array is over quark orderings */
    };


    //----------------------------------------------------------------------------
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
      write(bin, param.time_sources);
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
		 const std::string& ops_file)
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
	qqq.quark.resize(N_quarks);

	// Read 1-based spin
	multi1d<int> spin(N_quarks);
	reader >> spin[0] >> spin[1] >> spin[2];

	// Read 1-based displacement
	multi1d<int> displacement(N_quarks);
	reader >> displacement[0] >> displacement[1] >> displacement[2];

	// Insert for each quark
	for(int i=0; i < qqq.quark.size(); ++i)
	{
	  qqq.quark[i].spin = spin[i];
	  qqq.quark[i].displacement = displacement[i];
	}

      } //n
	
      reader.close();

      END_CODE();
    } //void readOps


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
      // Construct the dilution operator for each of the 3 quarks
      // 
      if (params.param.dilutions.size() != N_quarks)
      {
	QDPIO::cerr << name << ": expecting 3 quark dilutions" << endl;
	QDP_abort(1);
      }

      multi1d< Handle< DilutionOperator > > dilution_operators(N_quarks);  /*!< Here is the big (dil) pickle */

      try
      {
	// Loop over the 3 quark dilution operators
	for(int n=0; n < params.param.dilutions.size(); ++n)
	{
	  const GroupXML_t& dil_xml = params.param.dilutions[n];

	  std::istringstream  xml_d(dil_xml.xml);
	  XMLReader  diltop(xml_d);
	  QDPIO::cout << "Dilution type = " << dil_xml.id << endl;
	
	  dilution_operators[n] = TheFermDilutionOperatorFactory::Instance().createObject(
	    dil_xml.id,
	    diltop,
	    dil_xml.path);
	}
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception constructing dilution operator: " << e << endl;
	QDP_abort(1);
      }


      //
      // Initialize the slow Fourier transform phases
      //
      int decay_dir = dilution_operators[0]->getDecayDir();

      SftMom phases(params.param.mom2_max, false, decay_dir);
    
      // Sanity check - if this doesn't work we have serious problems
      if (phases.numSubsets() != QDP::Layout::lattSize()[decay_dir])
      {
	QDPIO::cerr << name << ": number of time slices not equal to that in the decay direction: " 
		    << QDP::Layout::lattSize()[decay_dir]
		    << endl;
	QDP_abort(1);
      }

		
      //
      // Another sanity check. The seeds of all the quarks must be different
      // THESE CAN NEVER EVER EVER BE THE SAME!! I MEAN IT (RGE)!!
      //
      for(int n=1; n < quarks.size(); ++n)
      {
	if ( toBool(dilution_operators[n]->getSeed() == dilution_operators[0]->getSeed()) ) 
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

      readOps(qqq_oplist, params.named_obj.operators_file.ops_file, params.param.displacement_length);


      //
      // Smear all the quark sources up front, overwrite the originals to save space,
      // and rotate them to the Dirac spin basis
      //
      Handle< QuarkSmearing<LatticeFermion> > sourceQuarkSmearing;

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
      Handle< QuarkSmearing<LatticeFermion> > sinkQuarkSmearing;
      try
      {
	QDPIO::cout << "Smear all the quark solutions up front" << endl;

	std::istringstream  xml_s(params.param.sink_quark_smearing.xml);
	XMLReader  smeartop(xml_s);
	
	sinkQuarkSmearing = 
	  TheFermSmearingFactory::Instance().createObject(params.param.sink_quark_smearing.id,
							  smeartop,
							  params.param.sink_quark_smearing.path);
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

      // The object holding the smeared and displaced spin components
      SmearedDispObjects disp_quarks(params.param.displacement_length,
				     dilution_operators,
				     

      //
      // Baryon operators
      //
      if (quarks.size() != N_quarks)
      {
	QDPIO::cerr << "expecting 3 quarks but have num quarks= " << quarks.size() << endl;
	QDP_abort(1);
      }

      //
      // Find the time slices that the source operator will be nonzero
      //
      std::list<int> participating_time_sources;
      for(int t0=0; t0 < QDP::Layout::lattSize()[decay_dir]; ++t0)
      {
	bool nonzero = true;

	for(int n=0; n < dilutions.size(); ++n)
	{
	  Handle< DilutionOperator<LatticeFermion> >& dilution = dilutions[n];

	  for(DilutionOperator<LatticeFermion>::const_iterator dil_ptr= dilution.begin(); 
	      dil_ptr != dilution.end(); 
	      ++dil_ptr)
	  {
	    if (dilution->hasTimeSupport(dil_ptr,t0))
	    {
	      nonzero = false;
	      break;
	    }
	  }
	}

	// Only if all dilutions have some support on this time slice is it enabled
	if (nonzero)
	  participating_time_sources.push_back(t0);
      }


      //
      // Permutations of quarks within an operator
      //
      const std::string& pstring = params.named_obj.quark_ids; 
      int num_orderings = 1; 
			
      //Should think of a cleverer algorithm for n quarks 
      if  (pstring.size() != N_quarks)
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
	multi1d<int> p(N_quarks);

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
      creat_oper.decay_dir   = decay_dir;
      creat_oper.seed_l      = dilutions[0]->getSeed();
      creat_oper.seed_m      = dilutions[1]->getSeed();
      creat_oper.seed_r      = dilutions[2]->getSeed();
      creat_oper.smearing    = params.param.source_quark_smearing;
      creat_oper.perms       = perms;
      creat_oper.orderings.resize(num_orderings);

      // Annihilation operator
      BaryonOperator_t  annih_oper;
      annih_oper.mom2_max    = params.param.mom2_max;
      annih_oper.decay_dir   = decay_dir;
      annih_oper.seed_l      = dilutions[0]->getSeed();
      annih_oper.seed_m      = dilutions[1]->getSeed();
      annih_oper.seed_r      = dilutions[2]->getSeed();
      annih_oper.smearing    = params.param.sink_quark_smearing;
      annih_oper.perms       = perms;
      annih_oper.orderings.resize(num_orderings);

      // Construct creation operator
      QDPIO::cout << "Build operators" << endl;
			
      // Loop over each operator 
      for(int l=0; l < qqq_oplist.ops.size(); ++l)
      {
	QDPIO::cout << "Elemental operator: op = " << l << endl;

	push(xml_out, "BaryonOperator");

	creat_oper.id = qqq_oplist.ops[l].name;
	annih_oper.id = qqq_oplist.ops[l].name;

	write(xml_out, "Name", creat_oper.id);

	// Loop over all orderings and build the operator
	swiss.reset();
	swiss.start();

	for(int ord = 0; ord < creat_oper.orderings.size(); ++ord)
	{
	  QDPIO::cout << "Ordering = " << ord << endl;

	  creat_oper.orderings[ord].perm = perms[ord];
	  creat_oper.orderings[ord].time_sources.resize(participating_time_sources.size());
		
	  annih_oper.orderings[ord].perm = perms[ord];
	  annih_oper.orderings[ord].time_sources.resize(participating_time_sources.size());
     
	  const int n0 = perms[ord][0];
	  const int n1 = perms[ord][1];
	  const int n2 = perms[ord][2];

	  // The operator must hold all the dilutions
	  // We know that all time slices match. However, not all time slices of the
	  // lattice maybe used
	  for(std::list<int>::const_iterator t_ptr= participating_time_sources.begin(); 
	      t_ptr != participating_time_sources.end(); 
	      ++t_ptr)
	  {
	    int t0 = *t_ptr;

	    // Creation operator
	    BaryonOperator_t::Orderings_t::TimeSources_t& cop = creat_oper.orderings[ord].time_sources[t];
	    
	    for(DilutionOperator<LatticeFermion>::const_iterator dil0= dilutions[n0].begin(); 
		dil0 != dilutions[n0].end(); 
		++dil0)
	    {
	      // Skip if a zero entry
	      if (dilutions[n0]->hasTimeSupport(dil0,t0))
		continue;

	      for(DilutionOperator<LatticeFermion>::const_iterator dil1= dilutions[n1].begin(); 
		  dil1 != dilutions[n1].end(); 
		  ++dil1)
	      {
		// Skip if a zero entry		// Skip if a zero entry
		if (dilutions[n1]->hasTimeSupport(dil1,t0))
		  continue;

		for(DilutionOperator<LatticeFermion>::const_iterator dil2= dilutions[n2].begin(); 
		    dil2 != dilutions[n2].end(); 
		    ++dil2)
		{
		  // Skip if a zero entry
		  if (dilutions[n2]->hasTimeSupport(dil2,t0))
		    continue;

		  // From now on we know the source operator is not zero
		  cop.t0 = t0;
		  
		  // The keys for the spin and displacements for this particular elemental operator
		  multi1d<KeySmearedDispColorVector_t> keySmearedDispColorVector(N_quark);
		  for(int n=0; n < N_quarks; ++n)
		  {
		    keySmearedDispColorVector[n].t0           = t0;
		    keySmearedDispColorVector[n].displacement = qqq_oplist.ops[l].quark[n].displacement;
		    keySmearedDispColorVector[n].spin         = qqq_oplist.ops[l].quark[n].spin;
		  }
	  
		  // Contract over color indices with antisym tensor.
		  // There is a potential optimization here - the colorcontract of
		  // the first two quarks could be pulled outside the innermost dilution
		  // loop.
		  // NOTE: the creation operator only lives on a few time slice, so restrict
		  // the operation to that time slice
				
		  LatticeComplex c_oper;
		  c_oper[phases.getSet()[t0]] =
		    colorContract(disp_quarks.dilutedSource(n0,dil0,keySmearedDispColorVector[0]),
				  disp_quarks.dilutedSource(n1,dil1,keySmearedDispColorVector[1]),
				  disp_quarks.dilutedSource(n2,dil2,keySmearedDispColorVector[2]));

		  // Slow fourier-transform
		  // We can restrict what the FT routine requires to a subset.
		  multi2d<DComplex> c_sum(phases.sft(c_oper, t0));

		  // Unpack into separate momentum and correlator
		  BaryonOperator_t::Orderings_t::TimeSlices_t::Dilutions_t cop_particular_dilution;
		  cop_particular_dilution.mom_projs.resize(phases.numMom());

		  for(int mom_num = 0 ; mom_num < phases.numMom() ; ++mom_num) 
		  {
		    cop_particular_dilution.mom_projs[mom_num].mom = phases.numToMom(mom_num);
		      
		    cop_particular_dilution.mom_projs[mom_num].op.resize(1);
		    cop_particular_dilution.mom_projs[mom_num].op[ 0 ] = c_sum[mom_num][ t0 ];
		  }

		  cop.dilutions.push_back(cop_particular_dilution);
		} // end for dil2
	      } // end for dil1
	    } // end for dil0
	    

	    // Annihilation operator
	    BaryonOperator_t::Orderings_t::TimeSources_t& aop = annih_oper.orderings[ord].time_sources[t];

	    for(DilutionOperator<LatticeFermion>::const_iterator dil0= dilutions[n0].begin(); 
		dil0 != dilutions[n0].end(); 
		++dil0)
	    {
	      // Skip if a zero entry
	      if (dilutions[n0]->hasTimeSupport(dil0,t0))
		continue;

	      for(DilutionOperator<LatticeFermion>::const_iterator dil1= dilutions[n1].begin(); 
		  dil1 != dilutions[n1].end(); 
		  ++dil1)
	      {
		// Skip if a zero entry		// Skip if a zero entry
		if (dilutions[n1]->hasTimeSupport(dil1,t0))
		  continue;

		for(DilutionOperator<LatticeFermion>::const_iterator dil2= dilutions[n2].begin(); 
		    dil2 != dilutions[n2].end(); 
		    ++dil2)
		{
		  // Skip if a zero entry
		  if (dilutions[n2]->hasTimeSupport(dil2,t0))
		    continue;

		  // From now on we know the source operator is not zero
		  aop.t0 = t0;

		  // The keys for the spin and displacements for this particular elemental operator
		  multi1d<KeySmearedDispColorVector_t> keySmearedDispColorVector(N_quark);
		  for(int n=0; n < N_quarks; ++n)
		  {
		    keySmearedDispColorVector[n].t0           = t0;
		    keySmearedDispColorVector[n].displacement = qqq_oplist.ops[l].quark[n].displacement;
		    keySmearedDispColorVector[n].spin         = qqq_oplist.ops[l].quark[n].spin;
		  }
	  
		  // Contract over color indices with antisym tensor.
		  // There is a potential optimization here - the colorcontract of
		  // the first two quarks could be pulled outside the innermost dilution
		  // loop.
		  LatticeComplex a_oper = 
		    colorContract(disp_quarks.dilutedSolution(n0,dil0,keySmearedDispColorVector[0]),
				  disp_quarks.dilutedSolution(n1,dil1,keySmearedDispColorVector[1]),
				  disp_quarks.dilutedSolution(n2,dil2,keySmearedDispColorVector[2]));
		  
		  // Slow fourier-transform
		  // We can restrict what the FT routine requires to a subset.
		  multi2d<DComplex> a_sum(phases.sft(a_oper));

		  // A new dilution holder for each iteration
		  BaryonOperator_t::Orderings_t::TimeSlices_t::Dilutions_t aop_particular_dilution;
		  BaryonOperator_t::Orderings_t::TimeSlices_t::Dilutions_t& aop_dilution = 
		    aop.dilutions.push_back(aop_particular_dilution);	

		  // Unpack into separate momentum and correlator
		  aop_dilution.mom_projs.resize(phases.numMom());

		  for(int mom_num = 0 ; mom_num < phases.numMom() ; ++mom_num) 
		  {
		    aop_dilution.mom_projs[mom_num].mom = phases.numToMom(mom_num);
		    aop_dilution.mom_projs[mom_num].op = a_sum[mom_num];
		  }

		} // end for dil2
	      } // end for dil1
	    } // end for dil0
	    
	  } // end for t
	} // end for ord

	swiss.stop();

		
	QDPIO::cout << "Operator construction: operator= " << l 
		    << "  time= "
		    << swiss.getTimeInSeconds() 
		    << " secs" << endl;

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

	QDPIO::cout << "Operator writing: operator = " << 
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
