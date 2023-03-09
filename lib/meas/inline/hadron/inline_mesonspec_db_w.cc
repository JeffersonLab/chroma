/*! \file
 * \brief Inline construction of hadron spectrum
 *
 * Spectrum calculations
 */


#include <sstream>

#include "util/ferm/key_val_db.h"
#include "util/ferm/key_hadron_2pt_corr.h"

#include "inline_mesonspec_db_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "io/param_io.h"
#include "io/qprop_io.h"
#include "meas/hadron/curcor2_w.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/smear/no_quark_displacement.h"
#include "util/ferm/meson_ops.h"

namespace Chroma 
{ 
  namespace InlineMesSpecEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMesSpec(InlineMesSpecParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "MESON_DB_SPECTRUM";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
  }

    
  //! Reader for parameters
  void read(XMLReader& xml, const std::string& path, InlineMesSpecParams::Param_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    switch (version) 
    {
    case 1:
      break;

    default:
      QDPIO::cerr << "Input parameter version " << version << " unsupported." << std::endl;
      QDP_abort(1);
    }

    //read(paramtop, "MesonP", param.MesonP);
    //read(paramtop, "CurrentP", param.CurrentP);
    //read(paramtop, "BaryonP", param.BaryonP);

    read(paramtop, "time_rev", param.time_rev);
    read(paramtop, "mom2_max", param.mom2_max);
    read(paramtop, "avg_equiv_mom", param.avg_equiv_mom);

    param.pz_only = false ;
    param.pz_max = -1;
    if( paramtop.count("pz_only")!=0 ){
      read(paramtop, "pz_only", param.pz_only);
      read(paramtop, "pz_max", param.pz_max);
    }

    read(paramtop, "ensemble", param.ensemble);
    read(paramtop, "States", param.states);
  }


  //! Writer for parameters
  void write(XMLWriter& xml, const std::string& path, const InlineMesSpecParams::Param_t& param)
  {
    push(xml, path);

    int version = 1;
    write(xml, "version", version);

    //write(xml, "MesonP", param.MesonP);
    //write(xml, "CurrentP", param.CurrentP);
    //write(xml, "BaryonP", param.BaryonP);

    write(xml, "time_rev", param.time_rev);

    write(xml, "mom2_max", param.mom2_max);
    write(xml, "avg_equiv_mom", param.avg_equiv_mom);
    write(xml, "pz_only", param.pz_only);
    write(xml, "pz_max", param.pz_max);
    write(xml, "ensemble", param.ensemble);
    write(xml, "States", param.states);

    pop(xml);
  }

  

  //! Propagator input
  void read(XMLReader& xml, const std::string& path, InlineMesSpecParams::NamedObject_t::Props_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "up_id", input.up_id);
    if(inputtop.count("down_id")!=0)
      read(inputtop, "down_id", input.down_id);
    else
      input.down_id = "NULL";
    if(inputtop.count("strange_id")!=0)
      read(inputtop, "strange_id", input.strange_id);
    else
      input.strange_id = "NULL";
    if(inputtop.count("charm_id")!=0)
      read(inputtop, "charm_id", input.charm_id);
    else
      input.charm_id = "NULL";
    if(inputtop.count("bottom_id")!=0)
      read(inputtop, "bottom_id", input.bottom_id);
    else
      input.bottom_id = "NULL";
    

    read(inputtop, "src", input.src);
    read(inputtop, "snk", input.snk);
  }

  //! Propagator output
  void write(XMLWriter& xml, const std::string& path, const InlineMesSpecParams::NamedObject_t::Props_t& input)
  {
    push(xml, path);

    write(xml, "up_id", input.up_id);
    write(xml, "down_id", input.down_id);
    write(xml, "strange_id", input.strange_id);
    write(xml, "charm_id", input.charm_id);
    write(xml, "bottom_id", input.bottom_id);
    
    write(xml, "src", input.src);
    write(xml, "snk", input.snk);

    pop(xml);
  }


  //! Propagator input
  void read(XMLReader& xml, const std::string& path, InlineMesSpecParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id" , input.gauge_id);
    read(inputtop, "props"    , input.props);
  }

  //! Propagator output
  void write(XMLWriter& xml, const std::string& path, const InlineMesSpecParams::NamedObject_t& input)
  {
    push(xml, path);
    write(xml, "gauge_id" , input.gauge_id);
    write(xml, "props"    , input.props);
    pop(xml);
  }


  // Param stuff
  InlineMesSpecParams::InlineMesSpecParams()
  { 
    frequency = 0; 
  }

  InlineMesSpecParams::InlineMesSpecParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;

      // Parameters for source construction
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
  InlineMesSpecParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    Chroma::write(xml_out, "Param", param);
    Chroma::write(xml_out, "NamedObject", named_obj);
    QDP::write(xml_out, "xml_file", xml_file);

    pop(xml_out);
  }


  // Anonymous namespace
  namespace 
  {
    //! Useful structure holding sink props
    class SinkPropContainer_t{
    public:
      ForwardProp_t prop_header;
      std::string quark_propagator_id;
      Real Mass;

      bool exists ;

      multi1d<int> bc; 
    
      std::string source_type;
      std::string source_disp_type;
      std::string sink_type;
      std::string sink_disp_type;

      SinkPropContainer_t(){
	exists=false ;
      }
              //! Read a sink prop
      void readSinkProp(const std::string& id)
      {
	if(id=="NULL"){
	  return ;
	}
	try
	  {
	    // Try a cast to see if it succeeds
	    const LatticePropagator& foo = 
	      TheNamedObjMap::Instance().getData<LatticePropagator>(id);
	    
	    // Snarf the data into a copy
	    quark_propagator_id = id;
	    
	    exists = true ;
	    // Snarf the prop info. This is will throw if the prop_id is not there
	    XMLReader prop_file_xml, prop_record_xml;
	    TheNamedObjMap::Instance().get(id).getFileXML(prop_file_xml);
	    TheNamedObjMap::Instance().get(id).getRecordXML(prop_record_xml);
	    
	    // Try to invert this record XML into a ChromaProp struct
	    // Also pull out the id of this source
	    {
	      std::string xpath;
	      read(prop_record_xml, "/SinkSmear", prop_header);
	      
	      read(prop_record_xml, "/SinkSmear/PropSource/Source/SourceType", source_type);
	      xpath = "/SinkSmear/PropSource/Source/Displacement/DisplacementType";
	      if (prop_record_xml.count(xpath) != 0)
		read(prop_record_xml, xpath, source_disp_type);
	      else
		source_disp_type = NoQuarkDisplacementEnv::getName();
	      
	      read(prop_record_xml, "/SinkSmear/PropSink/Sink/SinkType", sink_type);
	      xpath = "/SinkSmear/PropSink/Sink/Displacement/DisplacementType";
	      if (prop_record_xml.count(xpath) != 0)
		read(prop_record_xml, xpath, sink_disp_type);
	      else
		sink_disp_type = NoQuarkDisplacementEnv::getName();
	    }
	  }
	catch( std::bad_cast ) 
	  {
	    QDPIO::cerr << InlineMesSpecEnv::name << ": caught dynamic cast error" 
			<< std::endl;
	    QDP_abort(1);
	  }
	catch (const std::string& e) 
	  {
	    QDPIO::cerr << InlineMesSpecEnv::name << ": error message: " << e 
			<< std::endl;
	    QDP_abort(1);
	  }
	
	
	// Derived from input prop
	// Hunt around to find the mass
	// NOTE: this may be problematic in the future if actions are used with no
	// clear def. of a Mass
	QDPIO::cout << "Try action and mass" << std::endl;
	Mass = getMass(prop_header.prop_header.fermact);
	
	// Only baryons care about boundaries
	// Try to find them. If not present, assume dirichlet.
	// This turns off any attempt to time reverse which is the
	// only thing that the BC are affecting.
	bc.resize(Nd);
	bc = 0;
	
	try
	  {
	    bc = getFermActBoundary(prop_header.prop_header.fermact);
	  }
	catch (const std::string& e) 
	  {
	    QDPIO::cerr << InlineMesSpecEnv::name 
			<< ": caught exception. No BC found in these headers. Will assume dirichlet: " << e 
			<< std::endl;
	  }
	
	QDPIO::cout << "FermAct = " << prop_header.prop_header.fermact.id << std::endl;
	QDPIO::cout << "Mass = " << Mass << std::endl;
      }
      
    };





    void sanity_check_props(const SinkPropContainer_t& p1, 
			    const SinkPropContainer_t& p2){
      
      int j_decay = p1.prop_header.source_header.j_decay;
      if (p2.prop_header.source_header.j_decay != 
	  p1.prop_header.source_header.j_decay){
	QDPIO::cout<< "Error!! j_decay must be the same for all propagators"
		   << std::endl
		   << "     p2.j_decay= "<<p2.prop_header.source_header.j_decay
		   <<std::endl
		   << "     p1.j_decay= "<<p1.prop_header.source_header.j_decay
		   <<std::endl
		   << " Continuing hoping you know what you are doing!"
		   <<std::endl ;
	
	QDP_abort(1);
      }
      if (p2.prop_header.source_header.t_source != 
	  p1.prop_header.source_header.t_source)
	{
	  QDPIO::cerr<<"Error!! t_source must be the same for all propagators " << std::endl;
	  //	  QDPIO::cout<<"     p2.t_source= "<<p2.prop_header.source_header.t_source
	  //	     <<std::endl
	  //	     << "    p1.t_source= "<<p1.prop_header.source_header.t_source
	  //	     <<std::endl ;
	  QDP_abort(1);
	}
      if (p1.source_type != p2.source_type)
	{
	  QDPIO::cerr << "Warning: source_types not the same in a pair " << std::endl;
	  QDPIO::cerr << " Hope you know what your are doing... " << std::endl;
	  //	  QDP_abort(1);
	}
      if (p1.sink_type != p2.sink_type)
	{
	  QDPIO::cerr << "Error!! sink_type must be the same in a pair " << std::endl;
	  QDP_abort(1);
	}
      
      if (p2.bc[j_decay] != p1.bc[j_decay])
	{
	  QDPIO::cerr << "Error!! bc must be the same for all propagators " << std::endl;
	  QDP_abort(1);
	}
    }




    //! Useful structure holding sink props
    struct AllSinkProps_t
    {
      int j_decay;
      int t0 ; 
      multi1d<int> t_srce ;
      int bc_spec ;
      std::map<std::string,SinkPropContainer_t>  prop;
      std::map<std::string,LatticePropagator> rprop;

      //! Read all sinks
      AllSinkProps_t(const InlineMesSpecParams::NamedObject_t::Props_t& p){

	QDPIO::cout<<"Attempt to parse forward propagator= "<<p.up_id<<std::endl;
	prop["up"].readSinkProp(p.up_id);
	rprop[p.up_id]=TheNamedObjMap::Instance().getData<LatticePropagator>(p.up_id);
	QDPIO::cout<<"up quark  propagator successfully parsed" << std::endl;
	j_decay = prop["up"].prop_header.source_header.j_decay;
	t0      = prop["up"].prop_header.source_header.t_source;
	//t_srce  = prop["up"].prop_header.source_header.getTSrce();
	bc_spec = prop["up"].bc[j_decay];
	QDPIO::cout<<"Attempt to parse forward propagator= " <<p.down_id<<std::endl;
	//Always need a down quark 
	prop["down"].readSinkProp(p.down_id);
	QDPIO::cout << "down quark propagator successfully parsed" << std::endl;
	// get the source from the down because I may use the up for the
	// seqop
	t_srce  = prop["down"].prop_header.source_header.getTSrce();
	QDPIO::cout<<"   t_srce = ";
	for(int i(0);i<t_srce.size();i++)
	  QDPIO::cout<<" "<<t_srce[i] ;
	QDPIO::cout<<std::endl ;

	if(rprop.find(p.down_id) == rprop.end()){
	  QDPIO::cout<<__func__<<": Need to convert prop id: "
		     <<p.down_id<<std::endl;
	  rprop[p.down_id]=(TheNamedObjMap::Instance().getData<LatticePropagator>(p.down_id));
	}
	QDPIO::cout<<"Attempt to parse forward propagator= "<<p.strange_id<<std::endl;
	prop["strange"].readSinkProp(p.strange_id);
	if(p.strange_id != "NULL"){
	  QDPIO::cout <<"strange quark propagator successfully parsed" << std::endl;
	  if(rprop.find(p.strange_id) == rprop.end()){
	     QDPIO::cout<<__func__<<": Need to convert prop id: "
			<<p.strange_id<<std::endl;
	    rprop[p.strange_id]=(TheNamedObjMap::Instance().getData<LatticePropagator>(p.strange_id));
	  }
	}

	QDPIO::cout<<"Attempt to parse forward propagator= "<<p.charm_id<<std::endl;
	prop["charm"].readSinkProp(p.charm_id);
	if(p.charm_id != "NULL"){
	  QDPIO::cout << "charm quark propagator successfully parsed" << std::endl;
	  if(rprop.find(p.charm_id) == rprop.end()){
	    QDPIO::cout<<__func__<<": Need to convert prop id: "
		       <<p.charm_id<<std::endl;
	    rprop[p.charm_id]=(TheNamedObjMap::Instance().getData<LatticePropagator>(p.charm_id));
	  }
	}

	QDPIO::cout<<"Attempt to parse forward propagator= "<<p.bottom_id<<std::endl;
	prop["bottom"].readSinkProp(p.bottom_id);
	if(p.bottom_id != "NULL"){
	  QDPIO::cout << "bottom quark propagator successfully parsed" << std::endl;
	  if(rprop.find(p.bottom_id) == rprop.end()){
	    QDPIO::cout<<__func__<<": Need to convert prop id: "
		       <<p.bottom_id<<std::endl;
	    rprop[p.bottom_id]=(TheNamedObjMap::Instance().getData<LatticePropagator>(p.bottom_id));
	  }
	}
	
	 
	sanity_check_props(prop["up"], prop["down"]) ;
	if(prop["strange"].exists)
	  sanity_check_props(prop["up"], prop["strange"]) ;
	if(prop["charm"].exists)
	  sanity_check_props(prop["up"], prop["charm"]) ;
	if(prop["bottom"].exists)
	  sanity_check_props(prop["up"], prop["bottom"]) ;

      }


      std::string sink(const std::string& flavor){
	return prop[flavor].sink_type ;
      }

      std::string source(const std::string& flavor){
	return prop[flavor].source_type ;
      }

      const LatticePropagator& prop_ref(const std::string& flavor){

	return  rprop[prop[flavor].quark_propagator_id] ;
      
      }


    };


    

    void xml_print_prop_info(XMLWriter& xml_out, 
			     const SinkPropContainer_t& p,
			     const Set& s){
      
      // Sanity check - write out the norm2 of the forward prop 
      // Use this for any possible verification
      if(p.exists){
	push(xml_out, p.quark_propagator_id);
	{
	  const LatticePropagator& prop = 
	    TheNamedObjMap::Instance().getData<LatticePropagator>(p.quark_propagator_id);
	  
	  QDPIO::cout << "propagator_id = " << p.quark_propagator_id << std::endl;
	  QDPIO::cout << "Source_type   = " << p.source_type << std::endl;
	  QDPIO::cout << "Sink_type     = " << p.sink_type << std::endl;
	  
	  write(xml_out, "correlator",sumMulti(localNorm2(prop),s));
	  write(xml_out, "source_type", p.source_type);
	  write(xml_out, "source_disp_type", p.source_disp_type);
	  write(xml_out, "sink_type", p.sink_type);
	  write(xml_out, "sink_disp_type", p.sink_disp_type);
	  write(xml_out, "Mass", p.Mass);

	}
	pop(xml_out);
      }
    }



  } // namespace anonymous


  // Function call
  void 
  InlineMesSpec::operator()(unsigned long update_no,
			    XMLWriter& xml_out) 
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      std::string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "messpec");
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


  // Real work done here
  void 
  InlineMesSpec::func(unsigned long update_no,
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
      QDPIO::cerr << InlineMesSpecEnv::name << ": caught dynamic cast error" 
		  << std::endl;
      QDP_abort(1);
    }
    catch (const std::string& e) 
    {
      QDPIO::cerr << InlineMesSpecEnv::name << ": std::map call failed: " << e 
		  << std::endl;
      QDP_abort(1);
    }
    const multi1d<LatticeColorMatrix>& u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

    push(xml_out, "barspec");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << " MESSPEC: Spectroscopy for Wilson-like fermions" << std::endl;
    QDPIO::cout << std::endl << "     Gauge group: SU(" << Nc << ")" << std::endl;
    QDPIO::cout << "     volume: " << Layout::lattSize()[0];
    for (int i=1; i<Nd; ++i) {
      QDPIO::cout << " x " << Layout::lattSize()[i];
    }
    QDPIO::cout << std::endl;

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.write(xml_out, "Input");

    // Write out the config info
    write(xml_out, "Config_info", gauge_xml);

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 1);
    pop(xml_out);


    // First calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);

    //open the database
    //
    BinaryStoreDB< SerialDBKey<KeyHadron2PtCorr_t>, SerialDBData<multi1d<ComplexD> > > qdp_db;


    // Now loop over the various fermion pairs
    for(int p=0; p < params.named_obj.props.size(); ++p)
    {
      const InlineMesSpecParams::NamedObject_t::Props_t named_obj = params.named_obj.props[p];

      push(xml_out, "elem");

      AllSinkProps_t all_sinks(params.named_obj.props[p]);

      // Derived from input prop
      const multi1d<int>& t_srce = all_sinks.t_srce ;
      const int& j_decay         = all_sinks.j_decay ;
      const int& t0              = all_sinks.t0 ;

      const int& bc_spec         = all_sinks.bc_spec ;
     

      // Initialize the slow Fourier transform phases
      SftMom phases(params.param.mom2_max, t_srce, params.param.avg_equiv_mom,
                    j_decay);

      // Keep a copy of the phases with NO momenta
      SftMom phases_nomom(0, true, j_decay);

      write(xml_out, "t0", t0);

      // Save prop input
      push(xml_out, "Forward_prop_headers");
      write(xml_out, "up_prop", all_sinks.prop["up"].prop_header);
      write(xml_out, "down_prop", all_sinks.prop["down"].prop_header);
      if(all_sinks.prop["strange"].exists)
	write(xml_out, "strange_prop", all_sinks.prop["strange"].prop_header);
      if(all_sinks.prop["charm"].exists)
	write(xml_out, "charm_prop", all_sinks.prop["charm"].prop_header);
      if(all_sinks.prop["bottom"].exists)
	write(xml_out, "bottom_prop", all_sinks.prop["bottom"].prop_header);

      pop(xml_out);

      push(xml_out, "Forward_Propagator_Properties"                     );
      std::map<std::string,SinkPropContainer_t>::iterator it;
      for(it=all_sinks.prop.begin();it!=all_sinks.prop.end();it++)
	xml_print_prop_info(xml_out, it->second      , phases.getSet()  );
      
      pop(xml_out);

      int Nt = Layout::lattSize()[j_decay];

      LatticeComplex latC ;
      StopWatch tictoc;
      tictoc.reset();
      tictoc.start();
      for(int s(0);s<params.param.states.size();s++){
	QDPIO::cout<<"Doing state "<<params.param.states[s].name<<std::endl;
	//Key& qnums = params.param.states[s].q ;
	MesSpec::Key qnums   ;
	XMLReader wf_xml(params.param.states[s].wavefunc_file) ;
	multi1d<std::string> flavor ;
	multi1d<MesonOps::Operators_t> ops ;
	//	multi1d<SpinWF> quarkWF ;
	read(wf_xml,"/State/key",qnums);
	read(wf_xml,"/State/flavor",flavor);
	read(wf_xml,"/State/Operators",ops);
	


	QDPIO::cout<<"   Flavor structure: " ;
	for(int k(0);k<flavor.size();k++)
	  QDPIO::cout<<"  "<<flavor[k] ;
	QDPIO::cout<<std::endl;
	multi1d<std::string> prop_id(flavor.size()) ;
	for(int k(0);k<flavor.size();k++){
	  prop_id[k] = flavor[k] ;
	}
	// Open the file, and write the meta-data and the binary for this state
	if (! qdp_db.fileExists(params.param.states[s].db)){
	  XMLBufferWriter file_xml;
	  
	  push(file_xml, "DBMetaData");
	  write(file_xml, "id", std::string("hadron2Pt"));
	  write(file_xml, "lattSize", QDP::Layout::lattSize());
	  write(file_xml, "decay_dir", j_decay);
	  write(file_xml, "State", params.param.states[s].name);
	  proginfo(file_xml);    // Print out basic program info
	  write(file_xml, "Params", params.param);
	  write(file_xml, "Config_info", gauge_xml);
	  pop(file_xml);
	  
	  std::string file_str(file_xml.str());
	  qdp_db.setMaxUserInfoLen(file_str.size());	  
	  qdp_db.open(params.param.states[s].db, O_RDWR | O_CREAT, 0664);
	  qdp_db.insertUserdata(file_str);
	}
	else
	  {
	    qdp_db.open(params.param.states[s].db, O_RDWR, 0664);
	  }


	// References for use later 
	const LatticePropagator& q1 = all_sinks.prop_ref(prop_id[0]) ;
	const LatticePropagator& q2 = all_sinks.prop_ref(prop_id[1]) ;

	// the anti-quarks is the first
	LatticePropagator aq = Gamma(Ns*Ns-1)*adj(q1)*Gamma(Ns*Ns-1);


	KeyHadron2PtCorr_t key  ;

	std::ostringstream os ;
	os.precision(4);
	os<<std::fixed ;
	os<<"m"<<all_sinks.prop[prop_id[0]].Mass<<"_";
	os<<"m"<<all_sinks.prop[prop_id[1]].Mass;
	key.mass        = os.str();
	key.ensemble    = params.param.ensemble;
	key.num_vecs = Nc;  /*!< Number of color std::vectors */

	key.src_smear    =  params.named_obj.props[p].src ;
	key.src_spin  =  0 ; // not used params.param.states[s].spin ;
	key.src_lorentz = qnums.serialize() ;
	//key.src_lorentz =  key.src_spin  ;

	key.snk_smear    =  params.named_obj.props[p].snk ;
	key.snk_spin  =  0 ; // not used anymore params.param.states[s].spin ;
	key.snk_lorentz.resize(0); // this is also not needed
	//= key.src_lorentz ; // soure and sink have the same qnums
	//key.snk_lorentz =  key.snk_spin  ;

	//loop over momenta goes here
	for(int oi(0);oi<ops.size();oi++){ //sink
	  int snk(ops[oi].gamma);
	  for(int oj(0);oj<ops.size();oj++){//source

	    int src(ops[oj].gamma);

	    key.src_name    =  params.param.states[s].name + "-" + 
	                       ops[oj].name;
	    key.snk_name    =  params.param.states[s].name + "-" + 
	                       ops[oi].name;
	    
	    latC=trace(Gamma(snk)*q2*Gamma(src)*aq);

	    multi2d<ComplexD> hsum;
	  
	    // I  hard coded z to be the direction 2 ( 0 1 2 3 : x y z t)
	    int z_dir = 2 ; 
	    int Nmom = 2*params.param.pz_max + 1 ;
	    if( params.param.pz_max > Layout::lattSize()[z_dir]){
	      QDPIO::cout<<"   OOPS! pz_max > Lz "<<std::endl  ;
	      QDP_abort(121);
	    }
	    if(params.param.pz_only){
	      QDPIO::cout<<"   Doing momentum only in z"<<std::endl ;
	      hsum.resize( Nmom, phases.getSet().numSubsets() ) ;
	      for (int k(0) ; k<Nmom; k++){
		int mom = k - params.param.pz_max ;
		LatticeReal pz_z = (Layout::latticeCoordinate(z_dir) - 
				    t_srce[z_dir])
		  * mom * twopi / Real(Layout::lattSize()[z_dir]);
		LatticeComplex pp = cmplx(cos(pz_z),sin(pz_z));
		hsum[k] =  sumMulti(pp*latC, phases.getSet()) ;
	      }
	    }	      
	    else
	      hsum = phases.sft(latC) ;

	    
	    for(int mom(0);mom<hsum.size2();mom++){
	      if(params.param.pz_only){
		key.mom.resize(3) ;
		key.mom[0]=key.mom[1]=0;
		key.mom[2]=mom-params.param.pz_max ;
	      }
	      else
		key.mom = phases.numToMom(mom);    /*<! Momentum  */

	      SerialDBKey<KeyHadron2PtCorr_t> K;
	      K.key() = key ;
	      SerialDBData<multi1d<ComplexD> >  V ;
	      V.data().resize(Nt);
	      for(int t(0);t<Nt;t++){
		int t_eff = (t - t0 + Nt) % Nt;
		// Mesons do not need to flip the sign
		V.data()[t_eff] =  hsum[mom][t];
	      }//loop over time
	      qdp_db.insert(K,V);
	    }// loop over momenta

	  }// loop  over source ops
	}// loop over sink ops
	qdp_db.close() ;
      }// loop over states
      tictoc.stop();
      QDPIO::cout << InlineMesSpecEnv::name << ": contraction time = "
		  << tictoc.getTimeInSeconds() 
		  << " secs" << std::endl;

      pop(xml_out); //elem
    }// loop over propagator groups
    pop(xml_out);  // barspec

    snoop.stop();
    QDPIO::cout << InlineMesSpecEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << std::endl;

    QDPIO::cout << InlineMesSpecEnv::name << ": ran successfully" << std::endl;

    END_CODE();
  } 

};
