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
#include "meas/hadron/dilution_scheme_aggregate.h"
#include "meas/hadron/dilution_scheme_factory.h"
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
    void read(XMLReader& xml, const std::string& path, InlineStochGroupBaryonEnv::Params::Param_t& param)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      multi1d< multi1d<int> > temp;
      switch (version) 
      {
      case 1:
					
	QDPIO::cerr << "version 1 no longer supported. " << std::endl;
	exit(0);
	/*
	  read(paramtop, "mom2_max", param.mom2_max);


				
	  param.moms.resize(1,3);
				
	  param.moms[0][0] = 0;
	  param.moms[0][1] = 0;
	  param.moms[0][2] = 0;
					
	*/

	break;

      case 2:

	read(paramtop, "moms" , temp);

	param.mom2_max = 0;
	param.moms.resize(temp.size(), temp[0].size() );

	for (int i = 0 ; i < temp.size() ; ++i) 
	  param.moms[i] = temp[i];


	break;
				
      default :

	QDPIO::cerr << "Input parameter version " << version << " unsupported." << std::endl;
	QDP_abort(1);
      }

      read(paramtop, "displacement_length", param.displacement_length);

      param.quark_smearing = readXMLGroup(paramtop, "QuarkSmearing", "wvf_kind");
      param.link_smearing         = readXMLGroup(paramtop, "LinkSmearing", "LinkSmearingType");
      param.quark_dils             = readXMLArrayGroup(paramtop, "QuarkDilutions", "DilutionType");


    }


    // Writer for input parameters
    void write(XMLWriter& xml, const std::string& path, const InlineStochGroupBaryonEnv::Params::Param_t& param)
    {
      push(xml, path);

      int version = 1;

      write(xml, "version", version);
      write(xml, "mom2_max", param.mom2_max);
      write(xml, "displacement_length", param.displacement_length);
      xml << param.quark_smearing.xml;
      xml << param.link_smearing.xml;

      push(xml, "QuarkDilutions");

      for(int i=0; i < param.quark_dils.size(); ++i)
	xml << param.quark_dils[i].xml;

      pop(xml);

      pop(xml);
    }


    //Reader for 3-quark operator file
    void read(XMLReader& xml, const std::string& path, InlineStochGroupBaryonEnv::Params::NamedObject_t::ThreeQuarkOpsFile_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "ops_file", input.ops_file);
      read(inputtop, "id", input.id);
    }


    // Writer for 3-quark operator file
    void write(XMLWriter& xml, const std::string& path, const InlineStochGroupBaryonEnv::Params::NamedObject_t::ThreeQuarkOpsFile_t& input)
    {
      push(xml, path);
      write(xml, "ops_file", input.ops_file);
      write(xml, "id", input.id);
      pop(xml);
    }


    //! Read named objects 
    void read(XMLReader& xml, const std::string& path, InlineStochGroupBaryonEnv::Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "operators_file", input.operators_file);
      read(inputtop, "Quark_ids", input.quark_ids);
    }

    //! Write named objects
    void write(XMLWriter& xml, const std::string& path, const InlineStochGroupBaryonEnv::Params::NamedObject_t& input)
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
	success &= DilutionSchemeEnv::registerAll();
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
	QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
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

    //! Write quark 
    void write(XMLWriter& xml, const std::string& path, 
	       const ThreeQuarkOps_t::ThreeQuarkOp_t::QuarkInfo_t& input)
    {
      push(xml, path);

      write(xml, "Displacement", input.displacement);
      write(xml, "Spin", input.spin);

      pop(xml);
    }

    //! Write three quark op 
    void write(XMLWriter& xml, const std::string& path, 
	       const ThreeQuarkOps_t::ThreeQuarkOp_t& input)
    {
      push(xml, path);

      write(xml, "Quarks", input.quark);
      write(xml, "Name", input.name);

      pop(xml);
    }
	


    //---------------------------------------------------------------------------


    //! The key for smeared quarks 
    struct KeySmearedQuark_t
    {
      int  t0;              /*!< Time of source */
      int dil;              /*!< dilution component per timeslice */

    };


    //! Support for the keys of smeared quarks 
    bool operator<(const KeySmearedQuark_t& a, const KeySmearedQuark_t& b)
    {
      multi1d<int> lga(2);
      lga[0] = a.t0;
      lga[1] = a.dil;

      multi1d<int> lgb(2);
      lgb[0] = b.t0;
      lgb[1] = b.dil;

      return (lga < lgb);
    }


    struct SmearedQuark_t
    {
      LatticeFermion quark;
    };

    //----------------------------------------------------------------------------
    //! The key for smeared and displaced color vectors
    struct KeySmearedDispColorVector_t
    {
      int  t0;              /*!< Time of source */
      int dil;              /*!< dilution component per timeslice */

      int  displacement;    /*!< Orig plus/minus 1-based directional displacements */
      int  spin;            /*!< 1-based spin index */
    };


    //! Support for the keys of smeared and displaced color vectors
    bool operator<(const KeySmearedDispColorVector_t& a, const KeySmearedDispColorVector_t& b)
    {
      multi1d<int> lga(4);
      lga[0] = a.displacement;
      lga[1] = a.spin;
      lga[2] = a.t0;
      lga[3] = a.dil;

      multi1d<int> lgb(4);
      lgb[0] = b.displacement;
      lgb[1] = b.spin;
      lgb[2] = b.t0;
      lgb[3] = b.dil;

      return (lga < lgb);
    }


    //! The value of the std::map
    struct SmearedDispColorVector_t
    {
      multi1d<LatticeComplex> vec;
    };


    //----------------------------------------------------------------------------
    //! The smeared and displaced objects
    class SmearedDispObjects
    {
    public:
      //! Constructor from smeared std::map 
      SmearedDispObjects(int disp_length,
			 multi1d< Handle< DilutionScheme<LatticeFermion> > > dil_quarks,	
			 Handle< QuarkSmearing<LatticeFermion> > qsmr,
			 const multi1d<LatticeColorMatrix> & u_smr) :
	displacement_length(disp_length),diluted_quarks(dil_quarks),
	quarkSmearing(qsmr), u(u_smr)
	{
	  smeared_src_maps.resize(N_quarks);
	  smeared_soln_maps.resize(N_quarks);
			
	  disp_src_maps.resize(N_quarks);
	  disp_soln_maps.resize(N_quarks);
				
	}

      //! Destructor
      ~SmearedDispObjects() {}

      //! Accessor
      virtual multi1d<LatticeComplex> getDispSource(int quark_num, 						   
						    const KeySmearedDispColorVector_t& key);

      //! Accessor
      virtual multi1d<LatticeComplex> getDispSolution(int quark_num,  
						      const KeySmearedDispColorVector_t& key);

    protected:
      
      //! Displace an object
      virtual const multi1d<LatticeComplex>&	displaceObject(
	std::map<KeySmearedDispColorVector_t, SmearedDispColorVector_t>& disp_quark_map,
	const KeySmearedDispColorVector_t& key,
	const LatticeFermion& smrd_q);
			
      //! Smear sources and solutions
      virtual const LatticeFermion&
      smearSource(int qnum , const KeySmearedQuark_t & key);
					 
      virtual const LatticeFermion&
      smearSolution(int qnum , const KeySmearedQuark_t & key);


    private:
     
      multi1d< Handle< DilutionScheme<LatticeFermion> > > diluted_quarks;	
			
      Handle < QuarkSmearing<LatticeFermion> > quarkSmearing;

		
      //!Gauge field 
      const multi1d<LatticeColorMatrix> & u;
			
      //! Displacement length
      int displacement_length;

			
      //! Maps of smeared color vectors 
      multi1d< std::map<KeySmearedQuark_t, SmearedQuark_t> > smeared_src_maps;
      multi1d< std::map<KeySmearedQuark_t, SmearedQuark_t> > smeared_soln_maps;

      
      //!Maps of smeared displaced color vectors 
      multi1d< std::map<KeySmearedDispColorVector_t, SmearedDispColorVector_t> > disp_src_maps;
      multi1d< std::map<KeySmearedDispColorVector_t, SmearedDispColorVector_t> > disp_soln_maps;
    };

	
    const LatticeFermion&
    SmearedDispObjects::smearSource(int qnum , 
				    const KeySmearedQuark_t & key)
    {
      //
      // The spin basis matrix to goto Dirac
      //
      SpinMatrix rotate_mat(adj(DiracToDRMat()));

      std::map<KeySmearedQuark_t, SmearedQuark_t> & qmap = smeared_src_maps[qnum];

      //If entry is not in std::map create it
      if ( qmap.find(key) == qmap.end() )
      {

	// Insert an empty entry and then modify it. This saves on
	// copying the data around
	{
	  SmearedQuark_t smrd_empty;
	  qmap.insert(std::make_pair(key, smrd_empty));

	  // Sanity check - the entry better be there
	  if ( qmap.find(key) == qmap.end() )
	  {
	    QDPIO::cerr << __func__ 
			<< ": internal error - could not insert empty key in std::map"
			<< std::endl;
	    QDP_abort(1);
	  }		      
	}

	// Modify the previous empty entry
	SmearedQuark_t& smrd_q = qmap.find(key)->second;
	     
	StopWatch snoop;

	snoop.reset();
	snoop.start();
	
	smrd_q.quark = diluted_quarks[qnum]->dilutedSource(key.t0, key.dil);

	(*quarkSmearing)(smrd_q.quark, u);

	SpinMatrix mat = rotate_mat * Gamma(8);
				
	smrd_q.quark = mat * smrd_q.quark;
	
	      
	snoop.stop();

	QDPIO::cout << " Smeared Sources: Quark = "<< qnum <<" t0 = "
		    << key.t0 <<" dil = "<< key.dil << " Time = "<< snoop.getTimeInSeconds() <<" sec"<<std::endl;

	// Insert
	qmap.insert(std::make_pair(key, smrd_q));
      
      } // if find in std::map

      // The key now must exist in the std::map, so return the smeared quark field

      return (qmap.find(key)->second).quark ;
    }

    	
		
    const LatticeFermion&
    SmearedDispObjects::smearSolution(int qnum , 
				      const KeySmearedQuark_t & key)
    {
      //
      // The spin basis matrix to goto Dirac
      //
      SpinMatrix rotate_mat(adj(DiracToDRMat()));

      std::map<KeySmearedQuark_t, SmearedQuark_t> & qmap = smeared_soln_maps[qnum];


      //If entry is not in std::map create it
      if ( qmap.find(key) == qmap.end() )
      {

	// Insert an empty entry and then modify it. This saves on
	// copying the data around
	{
	  SmearedQuark_t smrd_empty;
	  qmap.insert(std::make_pair(key, smrd_empty));

	  // Sanity check - the entry better be there
	  if ( qmap.find(key) == qmap.end() )
	  {
	    QDPIO::cerr << __func__ 
			<< ": internal error - could not insert empty key in std::map"
			<< std::endl;
	    QDP_abort(1);
	  }		      
	}

	// Modify the previous empty entry
	SmearedQuark_t& smrd_q = qmap.find(key)->second;
	     
	StopWatch snoop;
	snoop.reset();
	snoop.start();
	
	smrd_q.quark = diluted_quarks[qnum]->dilutedSolution(key.t0, key.dil);

	(*quarkSmearing)(smrd_q.quark, u);

	smrd_q.quark = rotate_mat * smrd_q.quark;
	
	snoop.stop();

	QDPIO::cout << " Smeared Sinks: Quark = "<< qnum <<" t0 = "
		    << key.t0 <<" dil = "<< key.dil << " Time = "<< snoop.getTimeInSeconds() <<" sec"<<std::endl;
	
	// Insert
	qmap.insert(std::make_pair(key, smrd_q));
      
      } // if find in std::map

      // The key now must exist in the std::map, so return the smeared quark field

      return (qmap.find(key)->second).quark ;
    }

    //! Accessor
    multi1d<LatticeComplex>
    SmearedDispObjects::getDispSource(int quark_num, 
				      const KeySmearedDispColorVector_t& key)
    {
		
      //Get Smeared quark 
      KeySmearedQuark_t smr_key;
      smr_key.t0 = key.t0;
      smr_key.dil = key.dil;

      const LatticeFermion& smrd_q = smearSource(quark_num, smr_key);
					
      multi1d<LatticeComplex> vec;

      //Check if any displacement is needed
      if (displacement_length == 0) 
      {
	vec.resize(Nc);
	LatticeColorVector colvec = peekSpin( smrd_q, key.spin - 1);

	for (int c = 0 ; c < Nc ; ++c)
	  vec[c] = peekColor( colvec, c);

      }
      else
      {
	vec = displaceObject( disp_src_maps[quark_num] , key , smrd_q);
      }

      return vec;
		
    }

    //! Accessor
    multi1d<LatticeComplex>
    SmearedDispObjects::getDispSolution(int quark_num, 
					const KeySmearedDispColorVector_t& key)
    {
	
      //Get Smeared quark 
      KeySmearedQuark_t smr_key;
      smr_key.t0 = key.t0;
      smr_key.dil = key.dil;

      const LatticeFermion& smrd_q = smearSolution(quark_num, smr_key);
				
      multi1d<LatticeComplex> vec;

      //Check if any displacement is needed
      if (displacement_length == 0) 
      {
	vec.resize(Nc);
	LatticeColorVector colvec = peekSpin( smrd_q, key.spin - 1);

	for (int c = 0 ; c < Nc ; ++c)
	  vec[c] = peekColor( colvec, c);

      }
      else
      {
	vec = displaceObject( disp_soln_maps[quark_num] , key , smrd_q);
      }

      return vec;
    }


    //! Accessor
    const multi1d<LatticeComplex> &
    SmearedDispObjects::displaceObject(
      std::map<KeySmearedDispColorVector_t, SmearedDispColorVector_t>& disp_quark_map,
      const KeySmearedDispColorVector_t& key,
      const LatticeFermion& smrd_q)
    {

      StopWatch snoop;

      // If no entry, then create a displaced version of the quark
      if (disp_quark_map.find(key) == disp_quark_map.end())
      {
	//	      std::cout << __func__ 
	//		   << ": n=" << n
	//		   << " l=" << l
	//		   << " i=" << i 
	//		   << " disp=" << term.quark[i].displacement
	//		   << " len=" << term.quark[i].disp_len
	//		   << " dir=" << term.quark[i].disp_dir
	//		   << std::endl;



	// Insert an empty entry and then modify it. This saves on
	// copying the data around
	{
	  SmearedDispColorVector_t disp_empty;

	  snoop.reset();
	  snoop.start();

	  disp_quark_map.insert(std::make_pair(key, disp_empty));

	  snoop.stop();

	  QDPIO::cout<<"Inserted key in std::map: time = "<< snoop.getTimeInSeconds() << "secs"<<std::endl;

	  // Sanity check - the entry better be there
	  if (disp_quark_map.find(key) == disp_quark_map.end())
	  {
	    QDPIO::cerr << __func__ 
			<< ": internal error - could not insert empty key in std::map"
			<< std::endl;
	    QDP_abort(1);
	  }		      
	}

	// Modify the previous empty entry
	SmearedDispColorVector_t& disp_q = disp_quark_map.find(key)->second;

	snoop.reset();
	snoop.start();

	//Chroma uses a zero-based spin convention
	LatticeColorVector vec = peekSpin(smrd_q,  key.spin - 1);

	if (key.displacement > 0)
	{
	  int disp_dir = key.displacement - 1;
	  int disp_len = displacement_length;
	  displacement(u, vec, disp_len, disp_dir);
	}
	else if (key.displacement < 0)
	{
	  int disp_dir = -key.displacement - 1;
	  int disp_len = -displacement_length;
	  displacement(u, vec, disp_len, disp_dir);
	}

	snoop.stop();

	QDPIO::cout << "Displaced Quarks: Spin = "<<key.spin<<" Disp = "
		    << key.displacement <<" Time = "<<snoop.getTimeInSeconds() <<" sec"<<std::endl;

	disp_q.vec.resize(Nc);

	for(int i = 0 ; i < Nc ; ++i ) 
	{
	  disp_q.vec[i] = peekColor(vec, i);
	}

      } // if find in std::map

      snoop.reset();
      snoop.start();

      // The key now must exist in the std::map, so return the std::vector
      SmearedDispColorVector_t& disp_q = disp_quark_map.find(key)->second;

      snoop.stop(); 

      //QDPIO::cout << "Retrieved entry from std::map : time = "<< snoop.getTimeInSeconds() << "secs "<<std::endl;

      return disp_q.vec;
    }


    //----------------------------------------------------------------------------
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
	      multi1d<int>       mom;         /*!< D-1 momentum of this operator */
	      multi1d<DComplex>  op;          /*!< Momentum projected operator */
	    };

	    multi1d<Mom_t> mom_projs;         /*!< Holds momentum projections of the operator */

	  };

	  multi1d<int> perm;                  /*!< This particular permutation of quark orderings */

	  multi3d<Dilutions_t> dilutions;     /*!< Hybrid list indices */
	};

	multi1d<Orderings_t> orderings;  			/*!< Array is over quark orderings */

	int t0; 														/*!< Actual time corresponding to this timeslice */
      };

      multi1d< multi1d<int> > perms;   /*!< Permutations of quark enumeration */

      GroupXML_t    smearing;          /*!< String holding quark smearing xml */

      Seed          seed_l;            /*!< Id of left quark */
      Seed          seed_m;            /*!< Id of middle quark */
      Seed          seed_r;            /*!< Id of right quark */


      GroupXML_t    dilution_l;        /*!< Dilution scheme of left quark */ 
      GroupXML_t    dilution_m;        /*!< Dilution scheme of middle quark */ 
      GroupXML_t    dilution_r;        /*!< Dilution scheme of right quark */ 

      std::string   id;                /*!< Tag/ID used in analysis codes */

      int           mom2_max;          /*!< |\vec{p}|^2 */
      int           decay_dir;         /*!< Direction of decay */

      multi1d<TimeSlices_t> time_slices; /*!< Time slices of the lattice that are used */
    };


    //----------------------------------------------------------------------------
    //! BaryonOperator header writer
    void write(XMLWriter& xml, const std::string& path, const BaryonOperator_t& param)
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);
      write(xml, "id", param.id);
      write(xml, "mom2_max", param.mom2_max);
      write(xml, "decay_dir", param.decay_dir);
      write(xml, "seed_l", param.seed_l);
      write(xml, "seed_m", param.seed_m);
      write(xml, "seed_r", param.seed_r);
      write(xml, "perms", param.perms);

      push(xml, "dilution_l");
      xml << param.dilution_l.xml; 
      pop(xml);

      push(xml, "dilution_m");
      xml << param.dilution_m.xml; 
      pop(xml);
      
      push(xml, "dilution_r");
      xml << param.dilution_r.xml; 
      pop(xml);
			
      xml <<  param.smearing.xml;
    	
      pop(xml);
    }



    //! BaryonOperator binary writer
    void write(BinaryWriter& bin, const BaryonOperator_t::TimeSlices_t::Orderings_t::Dilutions_t::Mom_t& param)
    {
      write(bin, param.mom);
      write(bin, param.op);
    }

    //! BaryonOperator binary writer
    void write(BinaryWriter& bin, const BaryonOperator_t::TimeSlices_t::Orderings_t::Dilutions_t& param)
    {
      write(bin, param.mom_projs);
    }

    //! BaryonOperator binary writer
    void write(BinaryWriter& bin, const BaryonOperator_t::TimeSlices_t::Orderings_t& param)
    {
      write(bin, param.dilutions);
      write(bin, param.perm);
    }

    //! BaryonOperator binary writer
    void write(BinaryWriter& bin, const BaryonOperator_t::TimeSlices_t& param)
    {
      write(bin, param.orderings);
      write(bin, param.t0);
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
      write(bin, param.time_slices);
    }


    //! Read 3-quark operators file, assign correct displacement length
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
	std::string xml_file = makeXMLFileName(params.xml_file, update_no);

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
		    << std::endl;
	QDP_abort(1);
      }
      catch (const std::string& e) 
      {
	QDPIO::cerr << InlineStochGroupBaryonEnv::name << ": std::map call failed: " << e 
		    << std::endl;
	QDP_abort(1);
      }
      const multi1d<LatticeColorMatrix>& u = 
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

      push(xml_out, "StochGroupBaryon");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << InlineStochGroupBaryonEnv::name << ": Stochastic Baryon Operator" << std::endl;

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
      // Construct the dilution scheme for each of the 3 quarks
      // 
      if (params.param.quark_dils.size() != N_quarks)
      {
	QDPIO::cerr << name << ": expecting 3 quark dilutions" << std::endl;
	QDP_abort(1);
      }

      multi1d< Handle< DilutionScheme<LatticeFermion> > > diluted_quarks(N_quarks);  /*!< Here is the big (dil) pickle */

      try
      {
	// Loop over the 3 quark dilutions
	for(int n = 0; n < params.param.quark_dils.size(); ++n)
	{
	  const GroupXML_t& dil_xml = params.param.quark_dils[n];

	  std::istringstream  xml_d(dil_xml.xml);
	  XMLReader  diltop(xml_d);
	  QDPIO::cout << "Dilution type = " << dil_xml.id << std::endl;
	
	  diluted_quarks[n] = TheFermDilutionSchemeFactory::Instance().createObject(
	    dil_xml.id, diltop, dil_xml.path);
	}
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception constructing dilution scheme: " << e << std::endl;
	QDP_abort(1);
      }
//------------------------------------------------------------------------
//Sanity checks	

      //The participating timeslices must match for each quark
      //grab info from first quark to prime the work

      multi1d<int> participating_timeslices( diluted_quarks[0]->getNumTimeSlices() );

      for (int t0 = 0 ; t0 < participating_timeslices.size() ; ++t0)
      {
	participating_timeslices[t0] = diluted_quarks[0]->getT0(t0);
      }

      for (int n = 1 ; n < N_quarks ; ++n)
      {
	if ( diluted_quarks[n]->getNumTimeSlices() != participating_timeslices.size() )
	{
	  QDPIO::cerr << name << " : Quarks do not contain the same number of dilution timeslices: Quark " 
		      << n << std::endl; 

	  QDP_abort(1);
	}

	for (int t0 = 0 ; t0 < participating_timeslices.size() ; ++t0)
	{
	  if  ( diluted_quarks[n]->getT0(t0) != participating_timeslices[t0] )
	  {
	    QDPIO::cerr << name << " : Quarks do not contain the same participating timeslices: Quark "<<
	      n << " timeslice "<< t0 << std::endl;

	    QDP_abort(1);
	  }
	}
      }
		
      //Another Sanity check, the three quarks must all be 
      //inverted on the same cfg
      for (int n = 1 ; n < N_quarks ; ++n)
      {
	if (diluted_quarks[0]->getCfgInfo() != diluted_quarks[n]->getCfgInfo())
	{
	  QDPIO::cerr << name 
		      << " : Quarks do not contain the same cfg info, quark "<< n << std::endl;
	
	  QDP_abort(1);
	}
			
      }

      //Also ensure that the cfg on which the inversions were performed 
      //is the same as the cfg that we are using
      {	
	std::string cfgInfo; 

	//Really ugly way of doing this(Robert Heeeelp!!)
	XMLBufferWriter top;
	write(top, "Config_info", gauge_xml);
	XMLReader from(top);
	XMLReader from2(from, "/Config_info");
	std::ostringstream os;
	from2.print(os);

	cfgInfo = os.str();

	if (cfgInfo != diluted_quarks[0]->getCfgInfo())
	{
	  QDPIO::cerr << name 
		      << " : Quarks do not contain the same cfg info as the gauge field."
		      << "gauge: XX"<<cfgInfo<<"XX quarks: XX"<<diluted_quarks[0]->getCfgInfo()<<"XX"<<  std::endl;


	  QDP_abort(1);
	}
      }

      //
      // Initialize the slow Fourier transform phases
      //
      int decay_dir = diluted_quarks[0]->getDecayDir();

      //SftMom phases(params.param.mom2_max, false, decay_dir);
      //Changed this to cut down on the size of the files created 
      //Stupid way to do this.....
      //phases1 is the instance to be used if version one has been selected
//likewise for phases2
      //Idea: make a constructor for SftMom that takes a bit of 
      //xml and handles the problem. 

			
      SftMom phases(params.param.moms, decay_dir);

      // Sanity check - if this doesn't work we have serious problems
      if (phases.numSubsets() != QDP::Layout::lattSize()[decay_dir])
      {
	QDPIO::cerr << name << ": number of time slices not equal to that in the decay direction: " 
		    << QDP::Layout::lattSize()[decay_dir]
		    << std::endl;
	QDP_abort(1);
      }

		
      // Another sanity check. The seeds of all the quarks must be different
      // and thier decay directions must be the same 
      for(int n = 1 ; n < diluted_quarks.size(); ++n)
      {
	if ( toBool( diluted_quarks[n]->getSeed() == diluted_quarks[0]->getSeed() ) ) 
	{
	  QDPIO::cerr << name << ": error, quark seeds are the same" << std::endl;
	  QDP_abort(1);
	}

	if ( toBool( diluted_quarks[n]->getDecayDir() != diluted_quarks[0]->getDecayDir() ) )
	{
	  QDPIO::cerr << name << ": error, quark decay dirs do not match" <<std::endl;
	  QDP_abort(1);
	}

      }

      //
      // Smear the gauge field if needed
      //
      multi1d<LatticeColorMatrix> u_smr = u;

      try
      {
	std::istringstream  xml_l(params.param.link_smearing.xml);
	XMLReader  linktop(xml_l);
	QDPIO::cout << "Link smearing type = " << params.param.link_smearing.id << std::endl;
	
	
	Handle< LinkSmearing >
	  linkSmearing(TheLinkSmearingFactory::Instance().createObject(params.param.link_smearing.id,
								       linktop, params.param.link_smearing.path));

	(*linkSmearing)(u_smr);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception link smearing: " << e << std::endl;
	QDP_abort(1);
      }

      MesPlq(xml_out, "Smeared_Observables", u_smr);

      //Used for testing purposes 	
      multi1d<int> orig(4);
      for (int ind = 0 ; ind < 4 ; ++ind)
      {
	orig[ind] = 0;
      }

      //
      // Read operator coefficients
      //
      QDPIO::cout << "Reading 3-quark operators" << std::endl;
      ThreeQuarkOps_t qqq_oplist; 

      readOps(qqq_oplist, params.named_obj.operators_file.ops_file);

      //
      // Create the quark smearing factory 
      //
      Handle< QuarkSmearing<LatticeFermion> > quarkSmearing;

      try
      {
	QDPIO::cout << "Create quark smearing object" << std::endl;

	// Create the quark smearing object
	std::istringstream  xml_s(params.param.quark_smearing.xml);
	XMLReader  smeartop(xml_s);
	
	quarkSmearing =
	  TheFermSmearingFactory::Instance().createObject(params.param.quark_smearing.id,
							  smeartop, params.param.quark_smearing.path);
			
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << ": Caught Exception creating quark smearing object: " << e << std::endl;
	QDP_abort(1);
      }
      catch(...)
      {
	QDPIO::cerr << ": Caught generic exception creating smearing object" << std::endl;
	QDP_abort(1);
      }

      //
      // Baryon operators
      //

      //
      // Permutations of quarks within an operator
      //
      const std::string& pstr = params.named_obj.quark_ids; 
      int num_orderings = 1; 
			
      //Should think of a cleverer algorithm for n quarks 
      if  (pstr.size() != N_quarks)
      {
	QDPIO::cerr << "Invalid size for 'quark_ids'. Must be 3 but is " << pstr.size() << std::endl;
	QDP_abort(1);
      }

      //If 2 identical quarks, different one must be in the fisrt position
      if (  ( (pstr[0] == pstr[2]) && (pstr[1] != pstr[2]) ) ||
	    ( (pstr[0] == pstr[1]) && (pstr[1] != pstr[2]) ) )
      {
	QDPIO::cerr << "Invalid format for 'quark_ids'. Identical q's must be last 2 entries.: "
		    << pstr << std::endl;
	QDP_abort(1);
      }
			
      //Check that the kappas of the supposed identical quarks are the same.
      if ( (pstr[0] != pstr[1]) && (pstr[1] == pstr[2]) )
      {
	num_orderings = 2;
      
	if ( toBool(diluted_quarks[0]->getKappa() == diluted_quarks[1]->getKappa()) ||
	     toBool(diluted_quarks[1]->getKappa() != diluted_quarks[2]->getKappa()) )
	{
	  QDPIO::cerr << "quark_id's do not correspond to the correct identical quarks"
		      << std::endl;
	  QDP_abort(1);
	}
      }
			
      if  (pstr[0] == pstr[2]) 
      {
	num_orderings = 6;
      
	if ( toBool(diluted_quarks[0]->getKappa() != diluted_quarks[1]->getKappa()) ||
	     toBool(diluted_quarks[0]->getKappa() != diluted_quarks[2]->getKappa()) )
	{
					
	  QDPIO::cerr << "quark_id's do not correspond to the correct identical quarks"
		      << std::endl;
	  QDP_abort(1);
	}
      }
				
      QDPIO::cout << "Num Ordering = " << num_orderings << std::endl;

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

      //We make all source operators before we make all sink operators to 
      //save on memory. 

      for(int t0 = 0; t0 < participating_timeslices.size() ; ++t0)
      {
	StopWatch watch;
	//Make the source operators 
	{

	  // The object holding the smeared and displaced color std::vector std::maps  
	  SmearedDispObjects smrd_disp_srcs(params.param.displacement_length,
					    diluted_quarks, quarkSmearing, u_smr );

	  // Creation operator
	  BaryonOperator_t  creat_oper;
	  creat_oper.mom2_max    = 0;
	  creat_oper.decay_dir   = decay_dir;
	  creat_oper.seed_l      = diluted_quarks[0]->getSeed();
	  creat_oper.seed_m      = diluted_quarks[1]->getSeed();
	  creat_oper.seed_r      = diluted_quarks[2]->getSeed();
	  creat_oper.dilution_l  = params.param.quark_dils[0];
	  creat_oper.dilution_m  = params.param.quark_dils[1];
	  creat_oper.dilution_r  = params.param.quark_dils[2];
	  creat_oper.smearing    = params.param.quark_smearing;
	  creat_oper.perms       = perms;
	  creat_oper.time_slices.resize( 1 ); //Only a single t0 per file 

	  // Construct creation operator
	  // Loop over each operator 
	  for(int l=0; l < qqq_oplist.ops.size(); ++l)
	  {
	    QDPIO::cout << "Elemental operator: op = " << l << std::endl;

	    push(xml_out, "BaryonOperator");

	    creat_oper.id = qqq_oplist.ops[l].name;

	    write(xml_out, "Name", creat_oper.id);

	    // Loop over all orderings and build the operator
	    swiss.reset();
	    swiss.start();

	    // The keys for the spin and displacements for this particular elemental operator
	    multi1d<KeySmearedDispColorVector_t> keySmearedDispColorVector(N_quarks);

	    for(int n = 0 ; n < N_quarks ; ++n)
	    {
	      keySmearedDispColorVector[n].displacement = qqq_oplist.ops[l].quark[n].displacement;
	      keySmearedDispColorVector[n].spin         = qqq_oplist.ops[l].quark[n].spin;
	    }


	    creat_oper.time_slices[0].t0 = participating_timeslices[t0];
	    creat_oper.time_slices[0].orderings.resize(num_orderings);

	    for(int ord = 0; ord < num_orderings ; ++ord)
	    {
	      QDPIO::cout << "Ordering = " << ord << std::endl;

	      creat_oper.time_slices[0].orderings[ord].perm = perms[ord];

	      const int n0 = perms[ord][0];
	      const int n1 = perms[ord][1];
	      const int n2 = perms[ord][2];

	      // The operator must hold all the dilutions

	      // Creation operator
	      BaryonOperator_t::TimeSlices_t::Orderings_t& cop = creat_oper.time_slices[0].orderings[ord];

	      cop.dilutions.resize(diluted_quarks[n0]->getDilSize(t0), diluted_quarks[n1]->getDilSize(t0),
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

		  const multi1d<LatticeComplex> &q0 = smrd_disp_srcs.getDispSource(n0, 
										   keySmearedDispColorVector[0]); 

		  const multi1d<LatticeComplex> &q1 = smrd_disp_srcs.getDispSource(n1, 
										   keySmearedDispColorVector[1]);

		  /*QDPIO::cout<<"q0[0] testval= "<< peekSite(q0[0], orig)
		    << std::endl; 

		    QDPIO::cout<<"q1[0] testval= "<< peekSite(q1[0], orig)
		    << std::endl; 

		  */

		  watch.reset();
		  watch.start();
		  //For the source, restrict this operation to a subset
		  makeDiquark( diquark, q0 , q1, phases.getSet()[ participating_timeslices[t0] ] ); 
		  watch.stop();

		  /*QDPIO::cout<< " Made diquark : time = " << 
		    watch.getTimeInSeconds() << "secs" << std::endl;
		  */		

		  for(int k = 0 ; k < diluted_quarks[n2]->getDilSize(t0) ; ++k)	
		  {

		    keySmearedDispColorVector[2].dil = k;

		    // Contract over color indices with antisym tensor.
		    // NOTE: the creation operator only lives on a time slice, so restrict
		    // the operation to that time slice

		    LatticeComplex c_oper;

		    const multi1d<LatticeComplex> &q2 = smrd_disp_srcs.getDispSource(n2, 
										     keySmearedDispColorVector[2]);

		    /*QDPIO::cout<<"q2[0] testval= "<< peekSite(q2[0], orig)
		      << std::endl; 
		    */
		    watch.reset();
		    watch.start();

		    makeColorSinglet( c_oper, diquark, q2, phases.getSet()[ 
					participating_timeslices[t0] ] );

		    watch.stop();

		    /*QDPIO::cout<< "Made Color singlet : time =  " <<  
		      watch.getTimeInSeconds() << "secs" << std::endl;
		    */	
		    /*QDPIO::cout << "testval = " << peekSite(c_oper, orig) 
		      << std::endl;
		    */	
		    // Slow fourier-transform
		    // We can restrict what the FT routine requires to a subset.
		    watch.reset();
		    watch.start();

		    /*
		      multi2d<DComplex> c_sum(phases.sft(c_oper, 
		      participating_timeslices[t0] ));
		    */

		    multi2d<DComplex> c_sum;
		    int num_mom;

		    c_sum = phases.sft(c_oper, participating_timeslices[t0]);
		    num_mom = phases.numMom();
			
		    watch.stop();

		    /*QDPIO::cout << " Spatial sums completed: time = " << 
		      watch.getTimeInSeconds() << "secs " << std::endl;
		    */
		    // Unpack into separate momentum and correlator

		    cop.dilutions(i,j,k).mom_projs.resize(num_mom);

		    for(int mom_num = 0 ; mom_num < num_mom ; ++mom_num) 
		    {
		      cop.dilutions(i,j,k).mom_projs[mom_num].mom = 
			params.param.moms[mom_num];

		      cop.dilutions(i,j,k).mom_projs[mom_num].op.resize(1);

		      cop.dilutions(i,j,k).mom_projs[mom_num].op[ 0 ] = c_sum[mom_num][ 
			participating_timeslices[t0] ];
		    }

		  } // end for k
		} // end for j
	      } // end for i
	    }//end ord 

	    swiss.stop();


	    QDPIO::cout << "Source operator construction: operator= " << l 
			<< "  time= "
			<< swiss.getTimeInSeconds() 
			<< " secs" << std::endl;

	    QDPIO::cout << "Source operator testval(t0 = " << 
	      participating_timeslices[t0] << ") = " << 
	      creat_oper.time_slices[0].orderings[0].dilutions(0,0,0).mom_projs[0].op[0];



	    //Hard code the elemental op name for now 
	    std::stringstream cnvrt;
	    cnvrt <<  creat_oper.id  << "_t" << participating_timeslices[t0] << "_src.lime";

	    std::string filename;

	    filename = cnvrt.str(); 

	    // Write the meta-data and the binary for this operator
	    swiss.reset();
	    swiss.start();
	    {
	      XMLBufferWriter     src_record_xml, file_xml;
	      BinaryBufferWriter  src_record_bin;

	      push(file_xml, "SourceBaryonOperator");
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

		//	QDPIO::cout<< "t0 = " << t0 << " dil = "<< dil <<
		//	" srdhdr = XX"<<diluted_quarks[0]->getSourceHeader(t0,dil) << std::endl;
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

	      pop(file_xml); //Quark Sinks  
	      pop(file_xml);//SourceBaryonOp

	      QDPFileWriter qdp_file(file_xml, filename,     // are there one or two files???
				     QDPIO_SINGLEFILE, QDPIO_SERIAL, QDPIO_OPEN);


	      write(src_record_xml, "BaryonCreationOperator", creat_oper);
	      write(src_record_bin, creat_oper);

	      write(qdp_file, src_record_xml, src_record_bin);

	    }
	    swiss.stop();

	    QDPIO::cout << "Source Operator writing: operator = " << 
	      l	<< "  time= "
			<< swiss.getTimeInSeconds() 
			<< " secs" << std::endl;

	    pop(xml_out); // BaryonOperator 

	  } // end for l (operator )

	} //End Make creation operator



	//Make Annilation Operator
	{

	  // The object holding the smeared and displaced color std::vector std::maps  
	  SmearedDispObjects smrd_disp_snks(params.param.displacement_length,
					    diluted_quarks, quarkSmearing, u_smr );


	  // Annihilation operator
	  BaryonOperator_t  annih_oper;
	  annih_oper.mom2_max    = 0;
	  annih_oper.decay_dir   = decay_dir;
	  annih_oper.seed_l      = diluted_quarks[0]->getSeed();
	  annih_oper.seed_m      = diluted_quarks[1]->getSeed();
	  annih_oper.seed_r      = diluted_quarks[2]->getSeed();
	  annih_oper.dilution_l  = params.param.quark_dils[0];
	  annih_oper.dilution_m  = params.param.quark_dils[1];
	  annih_oper.dilution_r  = params.param.quark_dils[2];
	  annih_oper.smearing    = params.param.quark_smearing;
	  annih_oper.perms       = perms;
	  annih_oper.time_slices.resize( 1 );

	  // Construct annihilation operator
	  QDPIO::cout << "Building Sink operators" << std::endl;

	  // Loop over each operator 
	  for(int l=0; l < qqq_oplist.ops.size(); ++l)
	  {
	    QDPIO::cout << "Elemental operator: op = " << l << std::endl;

	    annih_oper.id = qqq_oplist.ops[l].name;

	    // Loop over all orderings and build the operator
	    swiss.reset();
	    swiss.start();

	    // The keys for the spin and displacements for this particular elemental operator
	    multi1d<KeySmearedDispColorVector_t> keySmearedDispColorVector(N_quarks);

	    for(int n = 0 ; n < N_quarks ; ++n)
	    {
	      keySmearedDispColorVector[n].displacement = qqq_oplist.ops[l].quark[n].displacement;
	      keySmearedDispColorVector[n].spin         = qqq_oplist.ops[l].quark[n].spin;
	    }

	    annih_oper.time_slices[0].t0 = participating_timeslices[t0];
	    //annih_oper.time_slices[0].orderings.resize(num_orderings);
	    annih_oper.time_slices[0].orderings.resize(1);


	    int ord = 0;
	    //for(int ord = 0 ; ord < num_orderings ; ++ord)
	    {
	      QDPIO::cout << "Ordering = " << ord << std::endl;

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
		  //	<< std::endl; 

		  //QDPIO::cout<<"q1[0] testval= "<< peekSite(q1[0], orig)
		  //	<< std::endl; 


		  watch.reset();
		  watch.start();

		  makeDiquark( diquark, q0 , q1, all ); 

		  watch.stop();
		  /*QDPIO::cout << "Made diquark: time = " << 
		    watch.getTimeInSeconds() << "secs " << std::endl;
		  */

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
		    //<< std::endl;

		    watch.reset();
		    watch.start();

		    makeColorSinglet( a_oper, diquark, q2, all);

		    watch.stop();

		    /*
		      QDPIO::cout <<	"Made Color Singlet: time = " <<
		      watch.getTimeInSeconds() << "secs" << std::endl;
		    */
		    /*QDPIO::cout << "testval = " << peekSite(a_oper, orig) 
		      << std::endl;
		    */

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
		      watch.getTimeInSeconds() << "secs" << std::endl;
		    */		
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
			<< " secs" << std::endl;

	    QDPIO::cout << "Sink op testval( t0 = " << 
	      participating_timeslices[t0] << ") = " << 
	      annih_oper.time_slices[0].orderings[0].dilutions(0,0,0).mom_projs[0].op[0] 
			<< std::endl;

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
			<< "  time= " << swiss.getTimeInSeconds() << " secs" << std::endl;

	  } // end for l (operator )

	} //End Make annihilation operator

      } //end t0 
      // Close the namelist output file XMLDAT
      pop(xml_out);     // StochBaryon

      snoop.stop();
      QDPIO::cout << InlineStochGroupBaryonEnv::name << ": total time = " 
		  << snoop.getTimeInSeconds() 
		  << " secs" << std::endl;

      QDPIO::cout << InlineStochGroupBaryonEnv::name << ": ran successfully" << std::endl;

      END_CODE();
    } // func

  } // namespace InlineStochGroupBaryonEnv

  /*! @} */  // end of group hadron

} // namespace Chroma
