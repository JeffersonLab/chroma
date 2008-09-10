// $Id: inline_qqqNucNuc_w.cc,v 3.17 2008-09-10 21:27:50 kostas Exp $
/*! \file
 * \brief The QQQ and QQBAR object calculation
 *
 */

#include "inline_qqqNucNuc_w.h"
#include "meas/hadron/qqq_w.h"
#include "meas/hadron/qqbar_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
//#include "io/param_io.h"
#include "io/qprop_io.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlineQQQNucNucEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineQQQNucNuc(InlineQQQNucNucParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "QQQ_NUCNUC";

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
  void read(XMLReader& xml, const string& path, InlineQQQNucNucParams::Param_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    switch (version) 
    {
    case 2:
      break;
    default:
      QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
      QDP_abort(1);
    }
    read(paramtop, "max_p2", param.max_p2);
    param.doVectorMesons = false ;
    param.doDecupletBar  = false ;
    if(paramtop.count("doVectorMesons") != 0 )
      read(paramtop, "doVectorMesons", param.doVectorMesons);
    if(paramtop.count("doDecupletBar") != 0 )
      read(paramtop, "doDecupletBar", param.doDecupletBar);

  }


  //! Writer for parameters
  void write(XMLWriter& xml, const string& path, const InlineQQQNucNucParams::Param_t& param)
  {
    push(xml, path);

    write(xml, "max_p2", param.max_p2);
    write(xml, "doVectorMesons", param.doVectorMesons);
    write(xml, "doDecupletBar", param.doDecupletBar);

    pop(xml);
  }


  //! Propagator input
  void read(XMLReader& xml, const string& path, InlineQQQNucNucParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "prop_ids", input.prop_ids);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlineQQQNucNucParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "prop_ids", input.prop_ids);

    pop(xml);
  }


  // Param stuff
  InlineQQQNucNucParams::InlineQQQNucNucParams() { frequency = 0; }

  InlineQQQNucNucParams::InlineQQQNucNucParams(XMLReader& xml_in, const std::string& path) 
  {
    try{
      XMLReader paramtop(xml_in, path);
      
      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;
      
      // Parameters for source construction
      read(paramtop, "Param", param);
      
      // Read in the output propagator/source configuration info
      read(paramtop, "NamedObject", named_obj);
      
      // Possible alternate qqq output file
      if (paramtop.count("qqq_file") != 0) 
	read(paramtop, "qqq_file", qqq_file);
      else // default qqq_file
	qqq_file = "DONTDO_qqq" ;

      // Possible alternate qqq output file
      if (paramtop.count("qqbar_file") != 0) 
	read(paramtop, "qqbar_file", qqbar_file);
      else // default qqq_file
	qqbar_file = "DONTDO_qqbar" ;
      
      // Possible alternate XML file pattern
      if (paramtop.count("xml_file") != 0) 
	read(paramtop, "xml_file", xml_file);
    }
    catch(const std::string& e){
      QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  void
  InlineQQQNucNucParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    Chroma::write(xml_out, "Param", param);
    Chroma::write(xml_out, "NamedObject", named_obj);
    QDP::write(xml_out, "qqq_file", qqq_file);
    QDP::write(xml_out, "xml_file", xml_file);

    pop(xml_out);
  }


  // Function call
  void 
  InlineQQQNucNuc::operator()(unsigned long update_no,
			      XMLWriter& xml_out) 
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != ""){
      string xml_file = makeXMLFileName(params.xml_file, update_no);
      
      push(xml_out, "qqqNucNuc_w");
      write(xml_out, "update_no", update_no);
      write(xml_out, "xml_file", xml_file);
      pop(xml_out);

      XMLFileWriter xml(xml_file);
      func(update_no, xml);
    }
    else
      func(update_no, xml_out);
  }


  // Real work done here
  void InlineQQQNucNuc::func(unsigned long update_no,
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
      QDPIO::cerr << InlineQQQNucNucEnv::name << ": caught dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineQQQNucNucEnv::name << ": map call failed: " << e 
		  << endl;
      QDP_abort(1);
    }
    const multi1d<LatticeColorMatrix>& u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

    push(xml_out, "qqqNucNuc_w");
    write(xml_out, "update_no", update_no);
    QDPIO::cout << " QQQNucNuc: Spectroscopy for Wilson fermions" << endl;
    
    // Type of sink smearing


    /*
     * Sanity checks
     */
    switch( params.named_obj.prop_ids.size()){
    case 1: //only up down
      break ;
    case 2: //only up down (0) and strange (1)
      break ;
    case 3:
      // 0 ---> up down only: Proton// neutron (up-down degenerate)
      // 1 ---> up down strange: Lambda 
      // 3 ---> up down strange charm:  charmed baryons and mesons
      break ;
    default:
      QDPIO::cerr << "OOOPS!! Don't know what to do with all theses propagators.... " << endl;
      QDP_abort(1);
    }

    QDPIO::cout << endl << "     Gauge group: SU(" << Nc << ")" << endl;


    QDPIO::cout << "     volume: " << Layout::lattSize()[0];
    for (int i=1; i<Nd; ++i) {
      QDPIO::cout << " x " << Layout::lattSize()[i];
    }
    QDPIO::cout << endl;
    

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.write(xml_out, "Input");
    

    // Write out the config info
    write(xml_out, "Config_info", gauge_xml);

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 1);
    pop(xml_out);

    /* I AM HERE */

    // First calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);
    
    
    multi1d<ForwardProp_t> quark_header(params.named_obj.prop_ids.size());
    multi1d<LatticePropagator> qprop(params.named_obj.prop_ids.size());
    multi1d<Real> Mass(params.named_obj.prop_ids.size());
    multi1d<string> sink_types(params.named_obj.prop_ids.size());
    multi2d<int> bc(params.named_obj.prop_ids.size(), 4); 

    // Now read the propagators we need
    for (int loop(0); loop < params.named_obj.prop_ids.size(); ++loop)
    {
      //	readQprop(prop_file_xml, prop_xml, qprop[loop],
      //	  params.prop.prop_files[loop], QDPIO_SERIAL);
	
      // Snarf the data into a copy
      qprop[loop] =
	TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop_ids[loop]);
	
      // Snarf the source info. 
      // This is will throw if the source_id is not there
      XMLReader prop_file_xml,prop_xml ;
      TheNamedObjMap::Instance().get(params.named_obj.prop_ids[loop]).getFileXML(prop_file_xml);
      TheNamedObjMap::Instance().get(params.named_obj.prop_ids[loop]).getRecordXML(prop_xml);
      // Try to invert this record XML into a ChromaProp struct
      // Also pull out the id of this source
      try
      {
	read(prop_xml, "/SinkSmear", quark_header[loop]);
	read(prop_xml, "/SinkSmear/PropSink/Sink/SinkType", sink_types[loop]);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << "Error extracting forward_prop header: " << e << endl;
	QDP_abort(1);
      }
	
      // Hunt around to find the mass and the boundary conditions
      // NOTE: this may be problematic in the future if actions are used with no
      // clear def. of a Mass
      QDPIO::cout << "Try action and mass" << endl;
      Mass[loop] = getMass(quark_header[loop].prop_header.fermact);
      bc[loop]   = getFermActBoundary(quark_header[loop].prop_header.fermact);
	
      QDPIO::cout << "FermAct = " << quark_header[loop].prop_header.fermact.path << endl;
      QDPIO::cout << "Mass = " << Mass[loop] << endl;
      QDPIO::cout << "boundary = "
		  << bc[loop][0]<<" "
		  << bc[loop][1]<<" "
		  << bc[loop][2]<<" "
		  << bc[loop][3]<< endl;
    }

    // Derived from input prop
    int j_decay = quark_header[0].source_header.j_decay;
    multi1d<int> t_srce = quark_header[0].source_header.getTSrce();
    //int t0     = t_source[j_decay];
    int t_source = quark_header[0].source_header.t_source;
    int t0 = t_source ;
    int bc_spec = bc[0][j_decay] ;
    for (int loop(0); loop < params.named_obj.prop_ids.size(); ++loop)
    {
      if(sink_types[loop]!=sink_types[0]){
	QDPIO::cerr << "Error!! sink types must be same for all propagators " << endl;
	QDP_abort(1);
      }
      if(quark_header[loop].source_header.j_decay!=j_decay){
	QDPIO::cerr << "Error!! j_decay must be the same for all propagators " << endl;
	QDP_abort(1);
      }
      if(bc[loop][j_decay]!=bc_spec){
	QDPIO::cerr << "Error!! bc must be the same for all propagators " << endl;
	QDP_abort(1);
      }
      for(int d(0);d<Nd;d++)
	if(quark_header[loop].source_header.t_source!=t_source){
	  QDPIO::cerr << "Error!! t_source must be the same for all propagators " << endl;
	  QDP_abort(1);
	}
    }
  
    // Establish the sink smearing type
    string sink_type;
    if (sink_types[0] == "POINT_SINK")
      sink_type = "POINT";
    else if (sink_types[0] == "SHELL_SINK")
      sink_type = "SHELL";
    else
    {
      QDPIO::cerr << InlineQQQNucNucEnv::name << ": unknown sink type = " << sink_types[0] << endl;
      QDP_abort(1);
    }

    // phases with  momenta
    SftMom phases(params.param.max_p2, t_srce, false, j_decay);
    
    
    push(xml_out,"Propagator_info") ;
    write(xml_out, "Masses", Mass);
    write(xml_out, "t_source", t_source);
    push(xml_out, "Propagator");
    for (int loop=0; loop < params.named_obj.prop_ids.size(); ++loop)
    {
      push(xml_out, "elem");
      write(xml_out, "ForwardProp",  quark_header[loop]);
      multi1d<Double> qp_corr= sumMulti(localNorm2(qprop[loop]),phases.getSet());
      push(xml_out, "Qprop_correlator");
      write(xml_out, "qp_corr", qp_corr);
      pop(xml_out); //Qprop_correlator
      pop(xml_out); //elem
    }
    pop(xml_out); //Propagator
    pop(xml_out); //Propagator_info
    
    XMLBufferWriter file_xml;
    push(file_xml, "qqqNucNuc_w");
    push(file_xml, "Output_version");
    write(file_xml, "out_version", 1);
    pop(file_xml);
    proginfo(file_xml);    // Print out basic program info
    params.write(file_xml, "Input");
    // Write out the config info
    write(xml_out, "Config_info", gauge_xml);
    push(file_xml,"Propagator_info") ;
    write(file_xml, "Masses", Mass);
    write(file_xml, "t_source", t_source);
    push(file_xml, "Propagator");
    for (int loop=0; loop < params.named_obj.prop_ids.size(); ++loop)
    {
      push(file_xml, "elem");
      write(file_xml, "ForwardProp",  quark_header[loop]);
      multi1d<Double> qp_corr= sumMulti(localNorm2(qprop[loop]),phases.getSet());
      push(file_xml, "Qprop_correlator");
      write(file_xml, "qp_corr", qp_corr);
      pop(file_xml); //Qprop_correlator
      pop(file_xml); //elem
    }
    pop(file_xml); //Propagator
    pop(file_xml); //Propagator_info
    write(file_xml, "MomNum", phases.numMom());
    pop(file_xml); //qqqNucNuc_w

    // Write the scalar data
    QDPFileWriter qqqto;
    QDPFileWriter qqbarto;
    if(params.qqq_file != "DONTDO_qqq")
      qqqto.open(file_xml,params.qqq_file, QDPIO_SINGLEFILE, QDPIO_SERIAL);
    //QDPIO_OPEN);

    if(params.qqbar_file != "DONTDO_qqbar")
      qqbarto.open(file_xml,params.qqbar_file, QDPIO_SINGLEFILE, QDPIO_SERIAL);
    // QDPIO_OPEN);
  

    multi2d<ThreeQuarks> qqq(phases.numMom(),phases.numSubsets()); 
    multi2d<DPropagator> qqbar(phases.numMom(),phases.numSubsets()); 
  

    {
      if(params.qqq_file != "DONTDO_qqq"){
	compute_qqq(qqq, qprop[0],qprop[0],qprop[0],phases,t0, bc_spec);
	write_qqq(qqqto, qqq, phases, "nucleon",sink_type);
	
	if(params.param.doDecupletBar)
	  for(int k(1);k<Nd;k++){
	    ostringstream tag ;
	    tag<<"Delta_"<<k;
	    compute_qqq(qqq, k,qprop[0],qprop[0],qprop[0],phases,t0, bc_spec);
	    write_qqq(qqqto, qqq, phases, tag.str(), sink_type);
	  }

      }

      if(params.qqbar_file != "DONTDO_qqbar"){
	compute_qqbar(qqbar, qprop[0],qprop[0],phases,t0 );
	write_qqbar(qqbarto, qqbar, phases, "pion",sink_type);
	if(params.param.doVectorMesons)
	  for(int k(0);k<Nd-1;k++){
	    ostringstream tag ;
	    tag<<"rho_"<<k;
	    compute_qqbar(qqbar, (1<<k),qprop[0],qprop[0],phases,t0 );
	    write_qqbar(qqbarto, qqbar, phases, tag.str(),sink_type);
	  }
      }

      if(params.named_obj.prop_ids.size()>1){
	if(params.qqq_file != "DONTDO_qqq"){
	  compute_qqq(qqq, qprop[0],qprop[0],qprop[1],phases,t0, bc_spec);
	  write_qqq(qqqto, qqq, phases, "lambda",sink_type);
	  compute_qqq(qqq, qprop[1],qprop[0],qprop[0],phases,t0, bc_spec);
	  write_qqq(qqqto, qqq, phases, "sigma",sink_type);
	  compute_qqq(qqq, qprop[0],qprop[1],qprop[1],phases,t0, bc_spec);
	  write_qqq(qqqto, qqq, phases, "xi",sink_type);

	  if(params.param.doDecupletBar)
	    for(int k(1);k<Nd;k++){
	      ostringstream tag ;
	      tag<<"Omega_"<<k;
	      compute_qqq(qqq,k,qprop[1],qprop[1],qprop[1],phases,t0, bc_spec);
	      write_qqq(qqqto, qqq, phases, tag.str(), sink_type);
	    }

	}
	
	if(params.qqbar_file != "DONTDO_qqbar"){
	  compute_qqbar(qqbar, qprop[0],qprop[1],phases,t0 );
	  write_qqbar(qqbarto, qqbar, phases, "kaon",sink_type);
	  compute_qqbar(qqbar, qprop[1],qprop[0],phases,t0 );
	  write_qqbar(qqbarto, qqbar, phases, "kaonbar",sink_type);
	  
	  if(params.param.doVectorMesons)
	    for(int k(0);k<Nd-1;k++){
	      ostringstream tag ;
	      tag<<"kaonst_"<<k;
	      compute_qqbar(qqbar, (1<<k), qprop[0],qprop[1],phases,t0 );
	      write_qqbar(qqbarto, qqbar, phases, tag.str(),sink_type);
	      compute_qqbar(qqbar, (1<<k), qprop[1],qprop[0],phases,t0 );
	      write_qqbar(qqbarto, qqbar, phases, "bar"+tag.str(),sink_type);
	    }
	}
      } // number of props > 1
      // do strange and charmed states 20' states
      if(params.named_obj.prop_ids.size()>2){
	if(params.qqq_file != "DONTDO_qqq"){
	  compute_qqq(qqq, qprop[0],qprop[0],qprop[2],phases,t0, bc_spec);
	  write_qqq(qqqto, qqq, phases, "lambda_c",sink_type);
	  compute_qqq(qqq, qprop[2],qprop[0],qprop[0],phases,t0, bc_spec);
	  write_qqq(qqqto, qqq, phases, "sigma_c",sink_type);
	  compute_qqq(qqq, qprop[0],qprop[1],qprop[2],phases,t0, bc_spec);
	  write_qqq(qqqto, qqq, phases, "xi_c",sink_type);
	  compute_qqq(qqq, qprop[2],qprop[1],qprop[0],phases,t0, bc_spec);
	  write_qqq(qqqto, qqq, phases, "xip_c",sink_type);
	  compute_qqq(qqq, qprop[2],qprop[1],qprop[1],phases,t0, bc_spec);
	  write_qqq(qqqto, qqq, phases, "omega_c",sink_type);
	  compute_qqq(qqq, qprop[0],qprop[2],qprop[2],phases,t0, bc_spec);
	  write_qqq(qqqto, qqq, phases, "xi_cc",sink_type);
	  compute_qqq(qqq, qprop[1],qprop[2],qprop[2],phases,t0, bc_spec);
	  write_qqq(qqqto, qqq, phases, "omega_cc",sink_type);
	}

	if(params.qqbar_file != "DONTDO_qqbar"){
	  compute_qqbar(qqbar, qprop[0],qprop[2],phases,t0 );
	  write_qqbar(qqbarto, qqbar, phases, "D",sink_type);
	  compute_qqbar(qqbar, qprop[2],qprop[0],phases,t0 );
	  write_qqbar(qqbarto, qqbar, phases, "Dbar",sink_type);
	  compute_qqbar(qqbar, qprop[1],qprop[2],phases,t0 );
	  write_qqbar(qqbarto, qqbar, phases, "Ds",sink_type);
	  compute_qqbar(qqbar, qprop[2],qprop[1],phases,t0 );
	  write_qqbar(qqbarto, qqbar, phases, "Dsbar",sink_type);

	  if(params.param.doVectorMesons)
	    for(int k(0);k<Nd-1;k++){
	      {
		ostringstream tag ;
		tag<<"Dst_"<<k;
		compute_qqbar(qqbar, (1<<k), qprop[0],qprop[2],phases,t0 );
		write_qqbar(qqbarto, qqbar, phases, tag.str(),sink_type);
		compute_qqbar(qqbar, (1<<k), qprop[2],qprop[0],phases,t0 );
		write_qqbar(qqbarto, qqbar, phases, "bar"+tag.str(),sink_type);
	      }
	      {
		ostringstream tag ;
		tag<<"Dsst_"<<k;
		compute_qqbar(qqbar, (1<<k), qprop[1],qprop[2],phases,t0 );
		write_qqbar(qqbarto, qqbar, phases, tag.str(),sink_type);
		compute_qqbar(qqbar, (1<<k), qprop[2],qprop[1],phases,t0 );
		write_qqbar(qqbarto, qqbar, phases, "bar"+tag.str(),sink_type);
	      }
	    }
	}
      }// end num props>2
    } 

    close(qqqto);
    close(qqbarto);
    
    pop(xml_out);  // qqqNucNuc_w
    

    snoop.stop();
    QDPIO::cout << InlineQQQNucNucEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << InlineQQQNucNucEnv::name << ": ran successfully" << endl;

    END_CODE();
  } 
  
}
