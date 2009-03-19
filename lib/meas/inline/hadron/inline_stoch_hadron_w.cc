// $Id: inline_stoch_hadron_w.cc,v 1.24 2009-03-19 17:12:20 mcneile Exp $
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

#include "util/ferm/key_val_db.h"

namespace Chroma{ 
  namespace InlineStochHadronEnv{ 
    enum HadronType {MESON_SRC_SRC, MESON_SOL_SOL,MESON_SRC_SOL,BARYON_SRC, BARYON_SOL} ;

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
	  
	  param.smearing = readXMLGroup(paramtop, "Smearing", "wvf_kind");
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
      xml<<param.smearing.xml;
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
      //read(inputtop, "meson_DB_file", input.meson_db);
      //read(inputtop, "baryon_DB_file", input.baryon_db);
    }
    
    //! Gauge field parameters
    void write(XMLWriter& xml, const string& path, const Params::NamedObject_t& input){
      push(xml, path);
      write(xml, "gauge_id", input.gauge_id);
      //write(xml, "meson_DB_file", input.meson_db);
      //write(xml, "baryon_DB_file", input.baryon_db);
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


    void ParseMeson(MesonOp& m, const GroupXML_t& grpXML){
      QDPIO::cout<<"I am a meson!\n" ;
      std::istringstream  xml_l(grpXML.xml);
      XMLReader xmltop(xml_l);
      XMLReader xml(xmltop,"/elem");
      QDPIO::cout << "Meson state is = " <<grpXML.id ; 
      QDPIO::cout << endl;
      QDPIO::cout << " XML =" <<grpXML.xml ; 
      QDPIO::cout << endl;
      read(xml,"Gamma",m.g);
      read(xml,"File",m.file);
    }

    void ParseBaryon(BaryonOp& m, const GroupXML_t& grpXML){
      QDPIO::cout<<"I am a Baryon!\n" ;
      std::istringstream  xml_l(grpXML.xml);
      XMLReader  xmltop(xml_l);
      XMLReader xml(xmltop,"/elem");
      QDPIO::cout << "Baryon state is = " <<grpXML.id ; 
      QDPIO::cout << endl;
      QDPIO::cout << " XML =" <<grpXML.xml ; 
      QDPIO::cout << endl;
      read(xml,"Gamma",m.g);
      read(xml,"File",m.file);
    }

    void meson(DComplex& corr,
	       const int& g,
	       const LatticeComplex& phase,
	       const LatticeFermion& eta,
	       const LatticeFermion& chi,
	       const Subset& s){
      corr = sum(localInnerProduct(eta,Gamma(g)*chi)*phase,s) ;
    }

    
    void baryon(multi1d<DComplex>& d,
		const int& g,
		const LatticeComplex& phase,
		const LatticeFermion& eta1,
		const LatticeFermion& eta2,
		const LatticeFermion& eta3,
		const Subset& s){
      //QDPIO::cout<<"I am a baryon!\n" ;
      if ( Nc != 3 ){    /* Code is specific to Ns=4 and Nc=3. */
	QDPIO::cerr<<"baryon code only works for Nc=3 and Ns=4\n";
	QDP_abort(111) ;
      }
#if QDP_NC == 3

      START_CODE();

      d.resize(Ns) ;

      d = zero;

      // C gamma_5 = Gamma(5)
      SpinMatrix g_one = 1.0 ;
      SpinMatrix Cg5 = Gamma(g)*g_one ; //BaryonSpinMats::Cg5();

      for(int k=0; k < Ns; ++k)
	{
	  LatticeSpinMatrix di_quark = zero;

	  for(int j=0; j < Ns; ++j)
	    {
	      for(int i=0; i < Ns; ++i)
		{
		  // Contract over color indices with antisym tensors
		  LatticeComplex b_oper = colorContract(peekSpin(eta1, i),
							peekSpin(eta2, j),
							peekSpin(eta3, k));

		  pokeSpin(di_quark, b_oper, j, i);
		}
	    }

	  d[k] += sum(traceSpin(Cg5 * di_quark)*phase,s);
	}

      END_CODE();
#endif
    }

    class Key{
    public:
      multi1d<int> k;

      Key& set(int i,int j, int l,int s){
	k.resize(4);
	k[0]=i ; k[1]=j ; k[2]=l ; k[3]=s;
	return *this ;
      }

      Key& set(int i,int j, int l){
	k.resize(3);
	k[0]=i ; k[1]=j ; k[2]=l ;
	return *this ;
      }

      Key& set(int i,int j){
	k.resize(2);
	k[0]=i ; k[1]=j ; 
	return *this ;
      }

      Key(){
	k.resize(1);
	k[0]=0;
      }
      Key(int i,int j, int l){
	set(i,j,l);
      }
      Key(int i,int j, int l,int s){
	set(i,j,l,s);
      }
      Key(int i,int j){
	set(i,j);
      }
      
      Key(const Key& klidi){
	k.resize(klidi.k.size()) ;
	for(int i(0);i<k.size();i++) k[i] = klidi.k[i] ;
      }
      //     int operator[](const int i){return k[i];} const 
      
      ~Key(){} ;
    };

    
    bool operator<(const Key& a, const Key& b){
      return (a.k<b.k) ;
    }
    //! Key binary writer
    void write(BinaryWriter& bin, const Key& klidi){
      write(bin, klidi.k);
    }

    //! Key binary reader
    void read(BinaryReader& bin, Key& klidi){
      read(bin, klidi.k);
    }


    struct HadronKey{
      // for the moment ignore displacements
      int type ; // creation: 0 or anihilation: 1
      int t    ; // time slice 
      int t0   ; // source time slice
      //if size of qn is 3 then it's a baryon if it is 2 then it's a meson
      multi1d<int> qn ; // the quark noise id 
      multi1d<int> p ;
      int gamma ; // the meson gamma matrix or the diquark baryon gamma matrix
    } ;

    //! HadronKey binary writer
    void write(BinaryWriter& bin, const HadronKey& h){
      write(bin, h.type);
      write(bin, h.t);
      write(bin, h.t0);
      write(bin, h.qn);
      write(bin, h.p);
      write(bin, h.gamma);
    }

    //! HadronKey binary reader
    void read(BinaryReader& bin, HadronKey& h){
      read(bin, h.type);
      read(bin, h.t);
      read(bin, h.t0);
      read(bin, h.qn);
      read(bin, h.p);
      read(bin, h.gamma);
    }

    /*
      intented use: 
      HadronOperator foo ;
      foo.data(Key(1,2)) ; // should return the 1,2 dilution of a meson
      foo.data(Key(1,2,3,s); for baryons where s is the spin index... 
     */
    struct HadronOperator{
      map<Key,DComplex> data ; // the Key has size 4 for a baryon of 2 for a hadron
    } ;
    
     //! HadronKey binary writer
    void write(BinaryWriter& bin, HadronOperator& h){
      map<Key,DComplex>::iterator it;
      int count(0);
      for(it=h.data.begin();it!=h.data.end();it++){
	count++;
      }
      write(bin,count);
      for(it=h.data.begin();it!=h.data.end();it++){
	write(bin, it->first);
	write(bin, it->second);
      }
    }

    //! HadronKey binary reader
    void read(BinaryReader& bin, HadronOperator& h){
      int count ;
      read(bin,count);
      Key k ;
      for(int i(0);i<count;i++){
	read(bin, k);
	read(bin, h.data[k]);
      }
    }

    class MesonOpData{
    public:
      multi2d<DComplex> data ;
      MesonOpData(){}
      MesonOpData(int n ){
	data.resize(n,n) ;
      }
      void resize(int n){
	data.resize(n,n);
      }
    } ;
    class BaryonOpData{
    public:
      multi3d< multi1d<DComplex> > data ;
      BaryonOpData(){}
      BaryonOpData(int n ){
        data.resize(n,n,n) ;
      }
      void resize(int n){
	data.resize(n,n,n);
      } 
    } ;


    //! MesonOp reader
    void read(BinaryReader& bin, MesonOpData& param)
    {
      int n;
      read(bin, n);  
      param.data.resize(n,n);
  
      for(int i=0; i < n; ++i)
      {
	for(int j=0; j < n; ++j)
	{
	  read(bin, param.data(i,j));
	}
      }
    }

    //! MesonOp write
    void write(BinaryWriter& bin, const MesonOpData& param)
    {
      int n = param.data.size1();  // all sizes the same
      write(bin, n);
      for(int i=0; i < n; ++i)
      {
	for(int j=0; j < n; ++j)
	{
	  write(bin, param.data(i,j));
	}
      }
    }



    //! BaryonOp reader
    void read(BinaryReader& bin, BaryonOpData& param)
    {
      int n;
      read(bin, n);  
      param.data.resize(n,n,n);
  
      for(int i=0; i < n; ++i)
      {
	for(int j=0; j < n; ++j)
	{
	  for(int k=0; k < n; ++k)
	    {
	      read(bin, param.data(i,j,k));
	    }
	}
      }
    }

    //! BaryonOp write
    void write(BinaryWriter& bin, const BaryonOpData& param)
    {
      int n = param.data.size1();  // all sizes the same
      write(bin, n);
      for(int i=0; i < n; ++i)
      {
	for(int j=0; j < n; ++j)
	{
	  for(int k=0; k < n; ++k)
            {
	      write(bin, param.data(i,j,k));
	    }
	}
      }
    }


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
      
      XMLBufferWriter UserData_xml;

      push(UserData_xml, "StochHadron");
      proginfo(UserData_xml);

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
      params.write(UserData_xml, "Input");

      // Write out the config info
      write(xml_out, "Config_info", gauge_xml);
      write(UserData_xml,"Config_info",gauge_xml);
      
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
      int NumBarOrd(N_quarks*(N_quarks-1)*(N_quarks-2)) ;
      int NumMesOrd(N_quarks*N_quarks) ; // to allow for flavor singlets
      if(N_quarks<2){
	QDPIO::cout << name << ": Need at least 2" << endl;
	QDP_abort(1);
      }
      bool BuildMesons(NumMesOrd>0);
      if(BuildMesons){
	QDPIO::cout << name << ": Can build mesons" << endl;
	QDPIO::cout<<name<< ": Number of meson orderings: "<<NumMesOrd<<endl;
      }
      else
	NumMesOrd=0;

      bool BuildBaryons(NumBarOrd>0);
      if(BuildBaryons){
	QDPIO::cout << name << ": Can build baryons" << endl;
	QDPIO::cout<<name<< ": Number of baryon orderings: "<<NumBarOrd<<endl;
      }
      else
	NumBarOrd=0;
      
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

      write(UserData_xml,"decay_dir",decay_dir);
      push(UserData_xml,"Quarks");
      for(int k(0);k<quarks.size();k++){
	push(UserData_xml,"elem");
	//write(UserData_xml,"Seed",quarks[0]->getSeed());
	UserData_xml<<quarks[k]->getSeed();
	pop(UserData_xml);
      }
      pop(UserData_xml);
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

      // Parse the Hadron operators
      map<string, MesonOp> LocalMesonOps ;
      map<string, BaryonOp> LocalBaryonOps ;
      for(int o(0);o<params.param.ops.size();o++){
	QDPIO::cout<<"Found Hadron: "<<params.param.ops[o].id<<endl ;
	map<string, void (*)(MesonOp&,const GroupXML_t& )>::iterator 
	  iter=mesons.find(params.param.ops[o].id);
	if(iter != mesons.end()){
	  iter->second(LocalMesonOps[params.param.ops[o].id],
		       params.param.ops[o]);
	}
	else{// Maybe it's a baryon
	  map<string, void (*)(BaryonOp&, const GroupXML_t&)>::iterator 
	    it=baryons.find(params.param.ops[o].id);
	  if(it != baryons.end())
	    it->second(LocalBaryonOps[params.param.ops[o].id],
			  params.param.ops[o]);
	  
	  else{
	    QDPIO::cout<<" Operator: "<<params.param.ops[o].id ;
	    QDPIO::cout<<" is unkown " <<endl ;
	  }
	}
      }

      //We only do diagonal one smearing smearing
      // I could read off the smearing from the source header so that
      // it is guaranteed to be the same as the one in the source...
      // set up the smearing and then do it
      Handle< QuarkSmearing<LatticeFermion> >  Smearing ;
      try{
	std::istringstream  xml_l(params.param.smearing.xml);
	XMLReader  smrtop(xml_l);
	QDPIO::cout << "Quark smearing type = " <<params.param.smearing.id ; 
          QDPIO::cout << endl;
          
          Smearing = TheFermSmearingFactory::Instance().createObject(params.param.smearing.id, smrtop, params.param.smearing.path);
      }
      catch(const std::string& e){
	QDPIO::cerr <<name<< ": Caught Exception creating quark smearing object: " << e << endl;
	QDP_abort(1);
      }
      catch(...){
	QDPIO::cerr <<name<< ": Caught generic exception creating smearing object" << endl;
	QDP_abort(1);
      }

      multi1d< multi1d< multi1d<LatticeFermion> > > smearedSol(quarks.size());
      multi1d< multi1d< multi1d<LatticeFermion> > > src(quarks.size());
      for(int q(0);q< quarks.size() ;q++){
	smearedSol[q].resize(participating_timeslices.size());
	src[q].resize(participating_timeslices.size());
	for (int t0 = 0 ; t0 < participating_timeslices.size() ; ++t0){
	  smearedSol[q][t0].resize(quarks[q]->getDilSize(t0)) ;
	  src[q][t0].resize(quarks[q]->getDilSize(t0)) ;
	  for(int i = 0 ; i <  quarks[q]->getDilSize(t0) ; ++i){
	    smearedSol[q][t0][i] = quarks[q]->dilutedSolution(t0,i) ;
	    src[q][t0][i] = quarks[q]->dilutedSource(t0, i) ;
	    (*Smearing)(smearedSol[q][t0][i], u_smr);  
	  }
	}
      }
      // Solution vectors are now smeared
      		  
      QDPIO::cout<<"   Doing " << phases.numMom()<< "  momenta "<<endl ;
      pop(UserData_xml);//done with UserData_xml 
      //First do all the mesons
      //Make a loop over meson operators
      //source source creation
      {
	map<string, MesonOp>::iterator it ;
	for(it=LocalMesonOps.begin();it!=LocalMesonOps.end();it++){
	  MesonOp op = it->second ;
	  // DB storage
	  //BinaryFxStoreDB< SerialDBKey<HadronKey>, SerialDBData<HadronOperator > > 
	  BinaryStoreDB< SerialDBKey<HadronKey>, SerialDBData<MesonOpData > > qdp_db;
	  qdp_db.setMaxUserInfoLen(UserData_xml.str().size());
	  qdp_db.open(op.file, O_RDWR | O_CREAT, 0664);
	  qdp_db.insertUserdata(UserData_xml.str());

	  SerialDBKey<HadronKey> key ;
	  //HadronKey key ;
	  key.key().gamma = op.g ;
	  // loop over the momentum projection
	  for(int mom_num = 0 ; mom_num < phases.numMom() ; ++mom_num){
	    key.key().p = phases.numToMom(mom_num);
	    for (int t0 = 0 ; t0 < participating_timeslices.size() ; ++t0){
	      key.key().t0 = participating_timeslices[t0] ;
	      key.key().t = key.key().t0 ; // creation ops leave on one time slice only
	      key.key().qn.resize(2);
	      //first do the sources
	      key.key().type = MESON_SRC_SRC ;
	      QDPIO::cout<<"   Doing MESON_SRC_SRC"<<endl ;
	      for(int q(0);q< quarks.size() ;q++){
		key.key().qn[0]=q;
		for(int q1(0);q1< quarks.size() ;q1++)
		  if(q!=q1){
		    key.key().qn[1] = q1 ;
		    SerialDBData<MesonOpData > val ;
		    val.data().resize(quarks[q]->getDilSize(t0));
		    //SerialDBData< HadronOperator > val ;
		    QDPIO::cout<<"   quarks: "<<q<<" "<<q1 <<endl ;
		    for ( int d(0) ; d < quarks[q]->getDilSize(t0); d++){
		      //LatticeFermion quark_bar = quarks[q]->dilutedSource(t0,d);
		      LatticeFermion quark_bar = src[q][t0][d];
		      quark_bar = Gamma(Ns-1)*quark_bar ;
		      for ( int d1(0) ; d1 < quarks[q1]->getDilSize(t0); d1++){
			//QDPIO::cout<<" quark: "<<q<<" "<<q1 ;
			//QDPIO::cout<<" dilution: "<<d<<" "<<d1<<endl ;
			//LatticeFermion quark = quarks[q1]->dilutedSource(t0,d1);
			LatticeFermion quark = src[q1][t0][d1];
			DComplex cc ;
			//meson(cc,op.g,phases[mom_num],quark_bar,quark,   
			//      phases.getSet()[key.key().t]);
			//if(toBool(real(cc)!=0.0) && toBool(imag(cc)!=0.0))
			//  val.data().data[Key(d,d1)] = cc ;
			meson(val.data().data(d,d1),op.g,phases[mom_num],
			      quark_bar,quark, phases.getSet()[key.key().t]);
		      }// dilutions d
		    }// dilutions d1
		    qdp_db.insert(key, val);
		  } // quark 2
	      } // quark 1 
	      
	      key.key().type = MESON_SOL_SOL ;
	      QDPIO::cout<<"   Doing MESON_SOL_SOL"<<endl ;
	      for(int q(0);q< quarks.size() ;q++){
		key.key().qn[0]=q;
		for(int q1(q+1);q1< quarks.size() ;q1++){
		  key.key().qn[1] = q1 ;
		  QDPIO::cout<<"   quarks: "<<q<<" "<<q1<<endl ;
		  SerialDBData<MesonOpData > val ;
		  val.data().resize(quarks[q]->getDilSize(t0));
		  //SerialDBData< HadronOperator > val ;
		  for(int t(0);t<phases.numSubsets();t++){
		    key.key().t = t ;
		    for ( int d(0) ; d < quarks[q]->getDilSize(t0); d++){
		      LatticeFermion quark_bar = smearedSol[q][t0][d] ;
		      quark_bar = Gamma(Ns-1)*quark_bar ;
		      for ( int d1(0) ; d1 < quarks[q1]->getDilSize(t0); d1++){
			//QDPIO::cout<<" quark: "<<q<<" "<<q1 ;
			//QDPIO::cout<<" dilution: "<<d<<" "<<d1<<endl ;
			LatticeFermion quark =  smearedSol[q1][t0][d1] ;
			meson(val.data().data(d,d1),op.g,phases[mom_num],
			      quark_bar,quark, phases.getSet()[key.key().t]);
			//DComplex cc ;
			//meson(cc,op.g,phases[mom_num],quark_bar,quark,
			//      phases.getSet()[key.key().t]);
			//if(toBool(real(cc)!=0.0)&&toBool(imag(cc)!=0.0))
			//  val.data().data[Key(d,d1)] = cc ;
		      }// dilutions d
		    }// dilutions d1
		    qdp_db.insert(key, val);
		  } // loop over time
		} // quark 2
	      } // quark 1 
	      

	      key.key().type = MESON_SRC_SOL ;
	      QDPIO::cout<<"   Doing MESON_SRC_SOL"<<endl  ;
              for(int q(0);q< quarks.size() ;q++){
                key.key().qn[0]=q;
                for(int q1(0);q1< quarks.size() ;q1++){
                  key.key().qn[1] = q1 ;
		  QDPIO::cout<<"   quarks: "<<q<<" "<<q1<<endl ;
		  SerialDBData<MesonOpData > val ;
                  val.data().resize(quarks[q]->getDilSize(t0));
		  //SerialDBData< HadronOperator > val ;
		  for (int tt = 0 ; tt < participating_timeslices.size() ; ++tt){
		    int t = participating_timeslices[tt] ;
		    key.key().t = t ;
		    for ( int d(0) ; d < quarks[q]->getDilSize(tt); d++){
		      //LatticeFermion quark_bar = quarks[q]->dilutedSource(tt,d);
		      LatticeFermion quark_bar = src[q][tt][d] ;
		      for ( int d1(0) ; d1 < quarks[q1]->getDilSize(t0); d1++){
			//QDPIO::cout<<" quark: "<<q<<" "<<q1 ;
			//QDPIO::cout<<" dilution: "<<d<<" "<<d1<<endl ;
			LatticeFermion quark =  smearedSol[q1][t0][d1] ;
			//DComplex cc ;
			meson(val.data().data(d,d1),op.g,phases[mom_num],quark_bar,quark,
			      phases.getSet()[key.key().t]);
			//if(toBool(real(cc)!=0.0)&&toBool(imag(cc)!=0.0))
			//val.data().data[Key(d,d1)] = cc ;
		      }// dilutions d1                                   
		    }// dilutions d                                        
		    qdp_db.insert(key, val);
		  } // loop over time                                              
		} // quark 2                                                             
	      } // quark 1          
	      
	    }// t0
	  }//mom
	}// ops
      }// Done with Mesons

      //Now the Baryons
      {
	map<string, BaryonOp>::iterator it ;
	for(it=LocalBaryonOps.begin();it!=LocalBaryonOps.end();it++){
	  BaryonOp op = it->second ;
	  BinaryStoreDB< SerialDBKey<HadronKey>, SerialDBData<BaryonOpData > > qdp_db;
	  qdp_db.setMaxUserInfoLen(UserData_xml.str().size());
	  qdp_db.open(op.file, O_RDWR | O_CREAT, 0664);
	  qdp_db.insertUserdata(UserData_xml.str());

	  SerialDBKey<HadronKey> key ;
	  //HadronKey key ;
	  key.key().gamma = op.g ;
	  // loop over the momentum projection
	  for(int mom_num = 0 ; mom_num < phases.numMom() ; ++mom_num){
	    key.key().p = phases.numToMom(mom_num);
	    for (int t0 = 0 ; t0 < participating_timeslices.size() ; ++t0){
	      key.key().t0 = participating_timeslices[t0] ;
	      key.key().t = key.key().t0 ; // creation ops leave on one time slice only
	      key.key().qn.resize(3);
	      //first do the sources
	      key.key().type = BARYON_SRC ;
	      QDPIO::cout<<"   Doing BARYON_SRC "<<endl ;
	      for(int q0(0);q0< quarks.size() ;q0++){
		key.key().qn[0]=q0;
		for(int q1(0);q1< quarks.size() ;q1++)
		  if(q0!=q1){
		    key.key().qn[1]=q1;
		    for(int q2(0);q2< quarks.size() ;q2++)
		      if((q1!=q2)&&(q0!=q2)){		    
			key.key().qn[2] = q2 ;
			QDPIO::cout<<"   quarks: "<<q0<<" "<<q1<<" "<<q2<<endl ;
			//SerialDBData< HadronOperator > val ;
			SerialDBData< BaryonOpData > val ;
			val.data().resize(quarks[q0]->getDilSize(t0));
			for ( int d0(0) ; d0 < quarks[q0]->getDilSize(t0); d0++){
			  //LatticeFermion quark0 = quarks[q0]->dilutedSource(t0,d0);
			  LatticeFermion quark0 = src[q0][t0][d0];
			  for ( int d1(0) ; d1 < quarks[q1]->getDilSize(t0); d1++){
			    //LatticeFermion quark1 = quarks[q1]->dilutedSource(t0,d1);
			    LatticeFermion quark1 = src[q1][t0][d1];
			    for ( int d2(0) ; d2 < quarks[q2]->getDilSize(t0); d2++){
			      //QDPIO::cout<<" quarks: "<<q0<<" "<<q1<<" "<<q2 ;
			      //QDPIO::cout<<" dilution: "<<d0<<" "<<d1<<" "<<d2<<endl ;
			      //LatticeFermion quark2 = quarks[q2]->dilutedSource(t0,d2);
			      LatticeFermion quark2 = src[q2][t0][d2];
			      multi1d<DComplex> cc ;
			      baryon(cc,op.g,phases[mom_num],quark0,quark1,quark2,   
				     phases.getSet()[key.key().t]);
			      val.data().data(d0,d1,d2) = cc ;
			      //for(int s(0);s<Ns;s++)
			      //if(toBool(real(cc[s])!=0.0) && toBool(imag(cc[s])!=0.0))
			      //  val.data().data[Key(d0,d1,d2,s)] = cc[s] ;
			    }// dilutions d2
			  }//dilutions d1
			}// dilutions d0
			qdp_db.insert(key, val);
		      } // quark 2
		  } // quark 1 
	      } // quark 0 

	      
	      // Then next block needs work
	      key.key().type = BARYON_SOL ;
	      QDPIO::cout<<"   Doing BARYON_SOL "<<endl ;
	      for(int q0(0);q0< quarks.size() ;q0++){
		key.key().qn[0]=q0;
		for(int q1(q0+1);q1< quarks.size() ;q1++){
		  key.key().qn[1] = q1 ;
		  for(int q2(q1+1);q2< quarks.size() ;q2++){
		    key.key().qn[2] = q2 ;
		    QDPIO::cout<<"   quarks: "<<q0<<" "<<q1<<" "<<q2 <<endl ;
		    //SerialDBData< HadronOperator > val ;
		    SerialDBData< BaryonOpData > val ;
		    val.data().resize(quarks[q0]->getDilSize(t0));
		    for(int t(0);t<phases.numSubsets();t++){
		      key.key().t = t ;
		      for ( int d0(0) ; d0 < quarks[q0]->getDilSize(t0); d0++){
			LatticeFermion quark0 = smearedSol[q0][t0][d0] ;
			for ( int d1(0) ; d1 < quarks[q1]->getDilSize(t0); d1++){
			  LatticeFermion quark1 = smearedSol[q1][t0][d1] ;
			  for ( int d2(0) ; d2 < quarks[q2]->getDilSize(t0); d2++){
			    //QDPIO::cout<<" quarks: "<<q0<<" "<<q1<<" "<<q2 ;
			    //QDPIO::cout<<" dilution: "<<d0<<" "<<d1<<" "<<d2<<endl ;
			    LatticeFermion quark2 =  smearedSol[q2][t0][d2] ;
			    
			    multi1d<DComplex> cc ;
			    baryon(cc,op.g,phases[mom_num],quark0,quark1,quark2,
				   phases.getSet()[key.key().t]);
			    val.data().data(d0,d1,d2) = cc ;
			    //for(int s(0);s<Ns;s++)
			    //  if(toBool(real(cc[s])!=0.0) && toBool(imag(cc[s])!=0.0))
			    //  val.data().data[Key(d0,d1,d2,s)] = cc[s] ;
			  }// dilutions d0
			}// dilutions d1
		      }//dilutions d2
		      qdp_db.insert(key, val);
		    } // loop over time
		  } // quark 2
		} // quark 1 
	      }// quark 0
	      
	    }// t0
	  }//mom
	}// ops

      }// done with baryons
      
      
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

