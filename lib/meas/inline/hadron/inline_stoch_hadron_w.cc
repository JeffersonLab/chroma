// $Id: inline_stoch_hadron_w.cc,v 1.8 2008-04-24 05:41:51 kostas Exp $
/*! \file
 * \brief Inline measurement of stochastic hadron operator (mesons and baryons).
 *
 */

#include "handle.h"
#include "meas/inline/hadron/inline_stoch_hadron_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/smear/quark_smearing_factory.h"
#include "meas/smear/quark_smearing_aggregate.h"
#include "meas/smear/displacement.h"
#include "meas/sources/source_smearing_aggregate.h"
#include "meas/sources/source_smearing_factory.h"
#include "meas/sources/dilutezN_source_const.h"
#include "meas/sources/zN_src.h"
#include "meas/sinks/sink_smearing_aggregate.h"
#include "meas/sinks/sink_smearing_factory.h"
#include "meas/hadron/barspinmat_w.h"
#include "meas/hadron/baryon_operator_aggregate_w.h"
#include "meas/hadron/baryon_operator_factory_w.h"
#include "meas/hadron/dilution_scheme_aggregate.h"
#include "meas/hadron/dilution_scheme_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma{ 
  namespace InlineStochHadronEnv{ 
    namespace{
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMeas(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "STOCH_HADRON";

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
  
  /** OLD STUFF
  // Read solution files
  void read(XMLReader& xml, const string& path, InlineStochHadronParams::Flavor_t::TimeSlice_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "t", input.t);
    read(inputtop, "dilution_files", input.dilution_files);
  }

  // Write solution files
  void write(XMLWriter& xml, const string& path, const InlineStochHadronParams::Flavor_t::TimeSlice_t& input)
  {
    push(xml, path);
    write(xml, "t", input.t);
    write(xml, "dilution_files", input.dilution_files);
    pop(xml);s
  }

 // Read solution files
  void read(XMLReader& xml, const string& path, InlineStochHadronParams::Flavor_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "time_slices", input.time_slices);
  }

  // Write solution files
  void write(XMLWriter& xml, const string& path, const InlineStochHadronParams::Flavor_t& input)
  {
    push(xml, path);
    write(xml, "time_slices", input.time_slices);
    pop(xml);s
  }

  **/

  // Reader for input parameters
    void read(XMLReader& xml, const string& path, Params::Param_t& param){
      XMLReader paramtop(xml, path);
      
      int version;
      read(paramtop, "version", version);
      
      switch (version) 
	{
	case 1:
	  /************************************************************/
	  read(paramtop, "mom2_max", param.mom2_max);
	  
	  param.smearing = readXMLArrayGroup(paramtop, "Smearing", "wvf_kind");
	  param.displace = readXMLArrayGroup(paramtop, "Displacement","DisplacementType");
	  param.link_smear = readXMLGroup(paramtop, "LinkSmearing", "LinkSmearingType");

	  param.ops = readXMLArrayGroup(paramtop, "HadronOperators", "Type");
	  param.quarks = readXMLArrayGroup(paramtop, "Quarks", "DilutionType");
	  
	  break;
	  
	default :
	  /**************************************************************/

	  QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
	  QDP_abort(1);
	}
      
      
    }


  // Writter for input parameters
    void write(XMLWriter& xml, const string& path, const Params::Param_t& param){
      push(xml, path);
      
      int version = 1;
      
      write(xml, "version", version);
      write(xml, "mom2_max", param.mom2_max);
      xml << param.link_smear.xml;
      
      push(xml,"Smearing");
      for( int t(0);t<param.smearing.size();t++){
	push(xml,"elem");
	xml<<param.smearing[t].xml;
	pop(xml);
      }
      pop(xml);
      
      push(xml,"Displacement");
      for( int t(0);t<param.displace.size();t++){
	push(xml,"elem");
	xml<<param.displace[t].xml;
	pop(xml);
      }
      pop(xml);
      
      push(xml,"Quarks");
      for( int t(0);t<param.quarks.size();t++){
	push(xml,"elem");
	xml<<param.quarks[t].xml;
	pop(xml);
      }
      pop(xml);
      
      push(xml,"HadronOperators");
      for( int t(0);t<param.ops.size();t++){
	push(xml,"elem");
	xml<<param.ops[t].xml;
	pop(xml);
      }
      pop(xml);
      
      
      pop(xml);
    }


    //! Gauge field parameters
    void read(XMLReader& xml, const string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);
      
      read(inputtop, "gauge_id", input.gauge_id);
    }
    
    //! Gauge field parameters
    void write(XMLWriter& xml, const string& path, const Params::NamedObject_t& input){
      push(xml, path);
      
      write(xml, "gauge_id", input.gauge_id);
      pop(xml);
    }
    
    
    // Param stuff
    Params::Params(){ 
      frequency = 0;
      param.mom2_max = 0 ;
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


    void Params::write(XMLWriter& xml_out, const std::string& path) 
    {
      push(xml_out, path);
      
      // Parameters for source construction
      InlineStochHadronEnv::write(xml_out, "Param", param);
      
      // Write out the output propagator/source configuration info
      InlineStochHadronEnv::write(xml_out, "NamedObject", named_obj);

      pop(xml_out);
    }



    class Key{
    public:
      multi1d<int> k;
      Key(){
	k.resize(3);
	k[0]=k[1]=k[2]=0;
      }
      Key(int i,int j, int l){
	k.resize(3);
	k[0]=i ;  k[1]=j ; k[2]=l ;
      }
      Key& set(int i,int j, int l){
	k[0]=i ; k[1]=j ; k[2]=l ;
	return *this ;
      }
      Key(const Key& klidi){
	k.resize(3) ;
	k[0] = klidi.k[0] ;
	k[1] = klidi.k[1] ;
	k[2] = klidi.k[2] ;
      }
      //     int operator[](const int i){return k[i];} const 
      
      ~Key(){} ;
    };
    
    bool operator<(const Key& a, const Key& b){
      return (a.k<b.k) ;
    }
    
    //! Baryon operator
    struct HadronOperator_t{
      //! Baryon operator time slices
      struct TimeSlice_t{
	struct Mom_t{
	  struct Permut_t{
	    //! Baryon operator dilutions
	    struct Dilutions_t{
	      multi3d<DComplex> d;
	      //vector<vector<vector<complex<double> > > > d ;
	    } ;
	    multi1d<Dilutions_t> s;
	    Permut_t(){s.resize(4);}
	  } ;
	  map<Key,Permut_t> p ;
	} ;
	map<Key,Mom_t> m;
      } ;
      map<int, TimeSlice_t> t;
      
      GroupXML_t    smearing;       /*!< String holding quark smearing xml */
      
    multi1d<Seed> seed  ;          /*!< Id of quarks */
      
      std::string   id;                /*!< Tag/ID used in analysis codes */
      
      int           mom2_max;          /*!< |\vec{p}|^2 */
      int           time_dir;         /*!< Direction of decay */
    };

    //! BaryonOperator header writer
    void write(XMLWriter& xml, const string& path, const HadronOperator_t& param)
    {
      push(xml, path);
      
      int version = 1;
      write(xml, "version", version);
      write(xml, "mom2_max", param.mom2_max);
      write(xml, "time_dir", param.time_dir);
      write(xml, "seed", param.seed);
      xml <<  param.smearing.xml;
      
      pop(xml);
    }
    
    //! Key binary writer
    void write(BinaryWriter& bin, const Key& klidi){
      write(bin, klidi.k);
    }

  /***** *****************
  //! BaryonOperator binary writer
  void write(BinaryWriter& bin, const BaryonOperator_t::TimeSlice_t::Mom_t::Permut_t::Dilutions_t& dil){
    write(bin, dil.d);
  }
  
  //! BaryonOperator binary writer
  void write(BinaryWriter& bin, const BaryonOperator_t::TimeSlice_t::Mom_t::Permut_t& p){
    write(bin, p.s);
  }
  
  //! BaryonOperator binary writer
  void write(BinaryWriter& bin, const BaryonOperator_t::TimeSlice_t::Mom_t& mm){
    map<Key,BaryonOperator_t::TimeSlice_t::Mom_t::Permut_t>::iterator it;
    for(it = mm.p.begin();it != mm.p.end();it++){
      write(bin, it->first);
      write(bin, it->second);
    }
  }
  
  //! BaryonOperator binary writer
  void write(BinaryWriter& bin, const BaryonOperator_t::TimeSlice_t& tt){
    map<Key,BaryonOperator_t::TimeSlice_t::Mom_t>::iterator it;
    for(it = tt.m.begin();it != tt.m.end();it++){
      write(bin, it->first);
      write(bin, it->second);
    }
  }

  //! BaryonOperator binary writer
  void write(BinaryWriter& bin, const BaryonOperator_t& bo){
    int version = 1;
    write(bin, param.mom2_max);
    write(bin, param.time_dir);
    write(bin, param.seed);

    map<int,BaryonOperator_t::TimeSlice_t>::iterator it;
    for(bo.begin();it != bo.end();it++){
      write(bin, it->first);
      write(bin, it->second);
    }
  }
  
  ****************/

  //--------------------------------------------------------------
  // Function call
  //  void 
  //InlineStochHadron::operator()(unsigned long update_no,
  //				XMLWriter& xml_out) 
    void InlineMeas::operator()(unsigned long update_no,
				XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != ""){
	string xml_file = makeXMLFileName(params.xml_file, update_no);
	
	push(xml_out, "stoch_hadron");
	write(xml_out, "update_no", update_no);
	write(xml_out, "xml_file", xml_file);
	pop(xml_out);
	
	XMLFileWriter xml(xml_file);
	func(update_no, xml);
      }
      else{
	func(update_no, xml_out);
      }
    }
    
    
    // Function call
    //void 
    //InlineStochHadron::func(unsigned long update_no,
    //			  XMLWriter& xml_out) 
    void InlineMeas::func(unsigned long update_no,
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
	  QDPIO::cerr << InlineStochHadronEnv::name << ": caught dynamic cast error" 
		      << endl;
	  QDP_abort(1);
	}
      catch (const string& e) 
	{
	  QDPIO::cerr << name << ": map call failed: " << e 
		      << endl;
	  QDP_abort(1);
	}
      const multi1d<LatticeColorMatrix>& u = 
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
      
      push(xml_out, "stoch_hadron");
      write(xml_out, "update_no", update_no);
      
      QDPIO::cout << name << ": Stochastic Hadron Operator" << endl;
      
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
      
      int N_quarks = params.param.quarks.size() ;
      if(N_quarks<2){
	QDPIO::cout << name << ": Need at least 2" << endl;
	QDP_abort(1);
      }
      bool BuildMesons(false);
      if(N_quarks>1){
	QDPIO::cout << name << ": Can build mesons" << endl;
	BuildMesons = true ;
      }
      bool BuildBaryons(false);
      if(N_quarks>2){
	QDPIO::cout << name << ": Can build baryons" << endl;
	BuildBaryons = true ;
      }
      
      multi1d< Handle< DilutionScheme<LatticeFermion> > > quarks(N_quarks); 
      
      try{
	// Loop over quark dilutions
	for(int n(0); n < params.param.quarks.size(); ++n){
	  const GroupXML_t& dil_xml = params.param.quarks[n];
	  std::istringstream  xml_d(dil_xml.xml);
	  XMLReader  diltop(xml_d);
	  QDPIO::cout << "Dilution type = " << dil_xml.id << endl;
	  quarks[n] = 
	    TheFermDilutionSchemeFactory::Instance().createObject(dil_xml.id, 
								  diltop, 
								  dil_xml.path);
	}
      }
      catch(const std::string& e){
	QDPIO::cerr << name << ": Caught Exception constructing dilution scheme: " << e << endl;
	QDP_abort(1);
      }
      
      
      //-------------------------------------------------------------------
      //Sanity checks	

      //The participating timeslices must match for each quark
      //grab info from first quark to prime the work
      
      multi1d<int> participating_timeslices(quarks[0]->getNumTimeSlices());
      
      for (int t0 = 0 ; t0 < participating_timeslices.size() ; ++t0){
	participating_timeslices[t0] = quarks[0]->getT0(t0);
      }
      
      for (int n = 1 ; n < N_quarks ; ++n){
	if(quarks[n]->getNumTimeSlices() != participating_timeslices.size()){
	  QDPIO::cerr << name ;
	  QDPIO::cerr << " : Quarks do not contain the same number";
	  QDPIO::cerr << "of dilution timeslices: Quark " << n << endl; 
	  QDP_abort(1);
	}
	
	for (int t0 = 0 ; t0 < participating_timeslices.size() ; ++t0){
	  if(quarks[n]->getT0(t0) != participating_timeslices[t0]){
	    QDPIO::cerr << name << " : Quarks do not contain the same";
	    QDPIO::cerr << "participating timeslices: Quark ";
	    QDPIO::cerr << n << " timeslice "<< t0 << endl;
	    QDP_abort(1);
	  }
	}
      }
		
      //Another Sanity check, the three quarks must all be 
      //inverted on the same cfg
      for (int n = 1 ; n < N_quarks ; ++n){
	if (quarks[0]->getCfgInfo() != quarks[n]->getCfgInfo()){
	  QDPIO::cerr << name << " : Quarks do not contain the same cfg info";
	  QDPIO::cerr << ", quark "<< n << endl;
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
	
	if (cfgInfo != quarks[0]->getCfgInfo()){
	  QDPIO::cerr << name << " : Quarks do not contain the same";
	  QDPIO::cerr << " cfg info as the gauge field." ;
	  QDPIO::cerr << "gauge: XX"<<cfgInfo<<"XX quarks: XX" ;
	  QDPIO::cerr << quarks[0]->getCfgInfo()<<"XX"<<  endl;
	  QDP_abort(1);
	}
      }

      //
      // Initialize the slow Fourier transform phases
      //
      int decay_dir = quarks[0]->getDecayDir();
      
      SftMom phases(params.param.mom2_max, false, decay_dir);
      
      // Sanity check - if this doesn't work we have serious problems
      if (phases.numSubsets() != QDP::Layout::lattSize()[decay_dir]){
	QDPIO::cerr << name << ": number of time slices not equal to that";
	QDPIO::cerr << " in the decay direction: " 
		    << QDP::Layout::lattSize()[decay_dir]
		    << endl;
	QDP_abort(1);
      }

		
      // Another sanity check. The seeds of all the quarks must be different
      // and thier decay directions must be the same 
      for(int n = 1 ; n < quarks.size(); ++n){
	if(toBool(quarks[n]->getSeed()==quarks[0]->getSeed())){
	  QDPIO::cerr << name << ": error, quark seeds are the same" << endl;
	  QDP_abort(1);
	}

	if(toBool(quarks[n]->getDecayDir()!=quarks[0]->getDecayDir())){
	  QDPIO::cerr<<name<< ": error, quark decay dirs do not match" <<endl;
	  QDP_abort(1);
	}
      }
      
      // Smear the gauge field if needed
      //
      multi1d<LatticeColorMatrix> u_smr = u;

      try{
	std::istringstream  xml_l(params.param.link_smear.xml);
	XMLReader  linktop(xml_l);
	QDPIO::cout << "Link smearing type = " << params.param.link_smear.id ; 
	QDPIO::cout << endl;
	
	
	Handle<LinkSmearing> linkSmearing(TheLinkSmearingFactory::Instance().createObject(params.param.link_smear.id, linktop, params.param.link_smear.path));

	(*linkSmearing)(u_smr);
      }
      catch(const std::string& e){
	QDPIO::cerr<<name<<": Caught Exception link smearing: " << e << endl;
	QDP_abort(1);
      }

      MesPlq(xml_out, "Smeared_Observables", u_smr);

    /***** COMMENT OUT  *** IN HOPE IT WILL COMPILE

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

      QDPFileWriter to(file_xml, params.named_obj.prop.op_file,     // are there one or two files???
		       QDPIO_SINGLEFILE, QDPIO_SERIAL, QDPIO_OPEN);

      // Write the scalar data
      {
	XMLBufferWriter record_xml;
	write(record_xml, "SourceBaryonOperator", baryon_opA);
	write(to, record_xml, baryon_opA.serialize());
      }

      // Write the scalar data
      {
	XMLBufferWriter record_xml;
	write(record_xml, "SinkBaryonOperator", baryon_opB);
	write(to, record_xml, baryon_opB.serialize());
      }

      close(to);
    }
    ****************************************/
      swatch.stop();
      
      QDPIO::cout << "Operators written: time= "
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;
      
      // Close the namelist output file XMLDAT
      pop(xml_out);     // StochHadron
      
      snoop.stop();
      QDPIO::cout << name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;
      
      QDPIO::cout << name << ": ran successfully" << endl;
      
      END_CODE();
    } 
  }  // namespace InlineHadronEnv
}// namespace chroma
