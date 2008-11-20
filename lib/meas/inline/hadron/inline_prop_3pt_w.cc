// $Id: inline_prop_3pt_w.cc,v 1.10 2008-11-20 03:53:44 kostas Exp $
/*! \file
 * \brief Inline measurement 3pt_prop
 *
 */

#include "handle.h"
#include "meas/inline/hadron/inline_prop_3pt_w.h"
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
#include "util/info/unique_id.h"
#include "util/ferm/transf.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma{ 
  namespace InlineProp3ptEnv{ 
    namespace{
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMeas(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "PROPAGATOR_3PT";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	//success &= BaryonOperatorEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
  
    void read(XMLReader& xml, const string& path, Params::Operator_t& op){
      XMLReader paramtop(xml, path);
      
      read(paramtop, "gamma", op.gamma);
      read(paramtop, "p", op.p);
      read(paramtop, "t", op.t);
      if(paramtop.count("factor")>0)
	read(paramtop, "factor", op.f);
      else
	op.f = 1.0 ;
    }
    
    // Writter for input parameters                                           
    void write(XMLWriter& xml, const string& path, const Params::Operator_t& op){
    push(xml, path);

    write(xml, "p", op.p);
    write(xml, "t", op.t);
    write(xml, "gamma", op.gamma);
    write(xml, "factor", op.f);

    pop(xml);
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
	  read(paramtop,"op",param.op);
	  param.chi = readXMLArrayGroup(paramtop, "Quarks", "DilutionType");
	  
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
      write(xml, "op", param.op);
      
      push(xml,"Quarks");
      for( int t(0);t<param.chi.size();t++){
	push(xml,"elem");
	xml<<param.chi[t].xml;
	pop(xml);
      }
      pop(xml);
            
      pop(xml); // final pop
    }


    //! Gauge field parameters
  void read(XMLReader& xml, const string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);
      
      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "prop_id", input.prop_id);
      read(inputtop, "prop3pt_id", input.prop3pt_id);
    }
    
    //! Gauge field parameters
  void write(XMLWriter& xml, const string& path, const Params::NamedObject_t& input){
      push(xml, path);
      
      write(xml, "gauge_id", input.gauge_id);
      write(xml, "prop_id", input.prop_id);
      write(xml, "prop3pt_id", input.prop3pt_id);
      pop(xml);
    }
    
    
    // Param stuff
    Params::Params(){ 
      frequency = 0;
      param.op.gamma = 0 ;
      param.op.f = 1.0 ;
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
      InlineProp3ptEnv::write(xml_out, "Param", param);
      
      // Write out the output propagator/source configuration info
      InlineProp3ptEnv::write(xml_out, "NamedObject", named_obj);

      pop(xml_out);
    }

    

  //--------------------------------------------------------------
  // Function call
  //  void 
  //InlineProp3pt::operator()(unsigned long update_no,
  //				XMLWriter& xml_out) 
    void InlineMeas::operator()(unsigned long update_no,
				XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != ""){
	string xml_file = makeXMLFileName(params.xml_file, update_no);
	
	push(xml_out, "propagator_3pt");
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
    //InlineProp3pt::func(unsigned long update_no,
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
	  QDPIO::cerr << InlineProp3ptEnv::name << ": caught dynamic cast error" 
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
      
      push(xml_out, "prop_3pt");
      write(xml_out, "update_no", update_no);
      
      QDPIO::cout << name << ": Propagator for 3pt functions" << endl;
      
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
      
      int N_quarks = params.param.chi.size() ;

      multi1d< Handle< DilutionScheme<LatticeFermion> > > quarks(N_quarks);

      try{
	// Loop over quark dilutions
	for(int n(0); n < params.param.chi.size(); ++n){
	  const GroupXML_t& dil_xml = params.param.chi[n];
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
      

      int decay_dir = quarks[0]->getDecayDir();
      //
      // Initialize the slow Fourier transform phases
      //
      int mom2(0);
      for(int i(0);i<Nd-1;i++)
	mom2 += params.param.op.p[i]*params.param.op.p[i] ;
      
      //SftMom phases(params.param.mom2_max, false, decay_dir);
      SftMom phases(mom2, false, decay_dir);
      
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
 
      // Read the quark propagator and extract headers
      ChromaProp_t prop_header;
      PropSourceConst_t source_header;
      QDPIO::cout << "Attempt to read propagator info" << endl;
      try
	{
	  // Try the cast to see if this is a valid source
	  LatticePropagator& tt_prop =
	    TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop_id);
	  
	  //Snarf the source info. 
	  //This is will throw if the source_id is not there
	  XMLReader prop_file_xml, prop_record_xml;
	  TheNamedObjMap::Instance().get(params.named_obj.prop_id).getFileXML(prop_file_xml);
	  TheNamedObjMap::Instance().get(params.named_obj.prop_id).getRecordXML(prop_record_xml);
	  // Try to invert this record XML into a ChromaProp struct
	  // Also pull out the id of this source
	  {
	    read(prop_record_xml, "/Propagator/ForwardProp", prop_header);
	    read(prop_record_xml, "/Propagator/PropSource", source_header);
	  }
	}    
      catch (std::bad_cast){
	QDPIO::cerr << name << ": caught dynamic cast error" ;
	QDPIO::cerr << endl;
	QDP_abort(1);
      }
      catch (const string& e){
	QDPIO::cerr << name << ": error extracting prop_header: " << e << endl;
	QDP_abort(1);
      }
      const LatticePropagator& prop = 
	TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop_id);
 
      QDPIO::cout << "Propagator successfully read and parsed" << endl;

      //create the 3pt prop
      try{
	TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.prop3pt_id);
      }
      catch (std::bad_cast){
	QDPIO::cerr << name << ": caught dynamic cast error" << endl;
	QDP_abort(1);
      }
      catch (const string& e){
	QDPIO::cerr << name << ": error creating prop: " << e << endl;
	QDP_abort(1);
      }

      // Cast should be valid now
      LatticePropagator& prop_3pt = 
	TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop3pt_id);

      // Need to check if the propagator gauge field is the same as the rest...
      // NEED TO IMPLEMENT THIS

      //Print out here the operator details
      write(xml_out, "op", params.param.op);
      //write(xml_out, "t0", t0);
      {
	XMLBufferWriter top;
        push(top, "tt");
	write(top, "Operator", params.param.op);
        //write(top, "t0", t0);
	pop(top);
	XMLReader from(top);
	XMLReader from2(from, "/tt/Operator");
	std::ostringstream os;
	QDPIO::cout<<name<<" Operator: "<<endl ;
	from2.print(os);
	QDPIO::cout<<os.str()<<endl ;
	//QDPIO::cout<< " Time slice: "<<t0<<endl ;
	QDPIO::cout<<name<<" End Operator "<<endl ;
      }
 
      LatticePropagator tmp_prop = Gamma(params.param.op.gamma)*prop ;
      tmp_prop = params.param.op.f*tmp_prop;
          
      LatticeFermion ferm_3pt = zero ;
      LatticeFermion ferm ;
      LatticeComplex phase = phases[phases.momToNum(params.param.op.p)];
      int t0 = params.param.op.t;
      for(int s(0);s<Ns;s++)
	for(int c(0);c<Nc;c++){
	  QDPIO::cout<<" Doing quark color and spin: "<<c<<" "<<s <<endl ;
	  PropToFerm(tmp_prop,ferm,c,s) ;
	  ferm_3pt = zero;
	  int count(0);
	  for(int n(0);n<quarks.size();n++){
	    int i_t0(-100) ;
	    for (int tt(0) ; tt < quarks[n]->getNumTimeSlices() ; ++tt)
	      if(quarks[n]->getT0(tt) == t0) 
		i_t0 = tt;
	    if(i_t0!=-100){ //this quark contains the appropriate t
	      count++;

	      QDPIO::cout<<" Doing quark: "<<n <<endl ;
	      QDPIO::cout<<"   quark: "<<n <<" has "<<quarks[n]->getDilSize(i_t0);
	      QDPIO::cout<<" dilutions on time slice "<<t0<<endl ;
	      for(int i = 0 ; i <  quarks[n]->getDilSize(i_t0) ; ++i){
		QDPIO::cout<<"   Doing dilution : "<<i<<endl ;
	        LatticeComplex cc = 
		  phase*localInnerProduct(quarks[n]->dilutedSource(i_t0,i),ferm);
		ferm_3pt += sum(cc,phases.getSet()[t0])*quarks[n]->dilutedSolution(i_t0,i);
	      }
	      QDPIO::cout<<" Done with dilutions for quark: "<<n <<endl ;
	    }
	  }
	  if(count==0){
	    QDPIO::cerr<<name<< ": error, no appropriate time slice found" <<endl;
	    QDP_abort(1);
	    
	  }
	  ferm_3pt = ferm_3pt/Double(count) ; // compute the mean over noises
	  FermToProp(ferm_3pt,prop_3pt,c,s);
	  QDPIO::cout<<" Done with quark color and spin: "<<c<<" "<<s <<endl ;
	}
      
      // Sanity check - 
      // write out the propagator (pion) correlator in the Nd-1 direction
      {
	// Initialize the slow Fourier transform phases
	SftMom  phases(0, true, Nd-1);

	multi1d<Double> corr = sumMulti(localNorm2(prop_3pt),phases.getSet());

	push(xml_out, "Prop_correlator");
	write(xml_out, "prop_corr", corr);
	pop(xml_out);
      }


      // Save the propagator info
      try
	{
	  QDPIO::cout << "Start writing propagator info" << endl;

	  XMLBufferWriter file_xml;
	  push(file_xml, "propagator");
	  write(file_xml, "id", uniqueId());  // NOTE: new ID form
	  pop(file_xml);

	  XMLBufferWriter record_xml;
	  push(record_xml , "Propagator");
	  write(record_xml, "ForwardProp", prop_header);
	  write(record_xml, "PropSource", source_header);
	  write(record_xml, "Config_info", gauge_xml);
	  //write(record_xml, "t0",t0);
	  write(record_xml, "Operator",params.param.op);
	  pop(record_xml);
	  // Write the propagator xml info
	  TheNamedObjMap::Instance().get(params.named_obj.prop3pt_id).setFileXML(file_xml);
	  TheNamedObjMap::Instance().get(params.named_obj.prop3pt_id).setRecordXML(record_xml);

	  QDPIO::cout << "Propagator successfully updated" << endl;
	}
      catch (std::bad_cast){
	QDPIO::cerr << name << ": caught dynamic cast error" << endl;
	QDP_abort(1);
      }
      catch (const string& e){
	QDPIO::cerr << name << ": error extracting prop_header: " << e << endl;
	QDP_abort(1);
      }


      // Close the namelist output file XMLDAT
      pop(xml_out);     // Prop3pt
      
      snoop.stop();
      QDPIO::cout << name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;
      
      QDPIO::cout << name << ": ran successfully" << endl;
      
      END_CODE();
    } 
  }  // namespace InlineProp3ptEnv
}// namespace chroma
