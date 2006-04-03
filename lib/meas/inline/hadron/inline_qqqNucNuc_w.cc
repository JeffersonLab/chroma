// $Id: inline_qqqNucNuc_w.cc,v 3.0 2006-04-03 04:59:02 edwards Exp $
/*! \file
 * \brief The QQQ and QQBAR object calculation
 *
 */

#include "inline_qqqNucNuc_w.h"
#include "meas/hadron/qqq_w.h"
#include "meas/hadron/qqbar_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/smear/ape_smear.h"
#include "meas/smear/sink_smear2.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "io/param_io.h"
#include "io/qprop_io.h"
#include "util/gauge/taproj.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/inline/io/default_gauge_field.h"

namespace Chroma 
{ 
  namespace InlineQQQNucNucEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineQQQNucNuc(InlineQQQNucNucParams(xml_in, path));
    }

    const std::string name = "QQQ_NUCNUC";
    const bool registered = TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
  };



  //! Reader for parameters
  void read(XMLReader& xml, const string& path, InlineQQQNucNucParams::Param_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    read(paramtop, "Pt_snk", param.Pt_snk);
    read(paramtop, "Sl_snk", param.Sl_snk);
    
    read(paramtop, "wvf_kind", param.wvf_kind);
    read(paramtop, "wvf_param", param.wvf_param);
    read(paramtop, "wvfIntPar", param.wvfIntPar);
    
    if (param.wvf_param.size() != param.wvfIntPar.size())
      {
	QDPIO::cerr << "wvf_param size inconsistent with wvfintpar size"<<endl;
	QDP_abort(1);
      }
    
    read(paramtop, "max_p2", param.max_p2);
  }


  //! Writer for parameters
  void write(XMLWriter& xml, const string& path, const InlineQQQNucNucParams::Param_t& param)
  {
    push(xml, path);

    int version = 1;
    write(xml, "version", version);

    write(xml, "Pt_snk", param.Pt_snk);
    write(xml, "Sl_snk", param.Sl_snk);
    
    write(xml, "wvf_kind", param.wvf_kind);
    write(xml, "wvf_param", param.wvf_param);
    write(xml, "wvfIntPar", param.wvfIntPar);
    
    if (param.wvf_param.size() != param.wvfIntPar.size())
      {
	QDPIO::cerr << "wvf_param size inconsistent with wvfintpar size"<<endl;
	QDP_abort(1);
      }
    
    write(xml, "max_p2", param.max_p2);

    pop(xml);
  }


  //! Propagator input
  void read(XMLReader& xml, const string& path, InlineQQQNucNucParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    input.gauge_id = InlineDefaultGaugeField::readGaugeId(inputtop, "gauge_id");
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
	qqq_file = "qqq.dat" ;

      // Possible alternate qqq output file
      if (paramtop.count("qqbar_file") != 0) 
	read(paramtop, "qqbar_file", qqbar_file);
      else // default qqq_file
	qqbar_file = "qqbar.dat" ;
      
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
    
    /*
     * Sanity checks
     */
    if (params.param.wvf_param.size() != 1){
      QDPIO::cerr << "only one wvf_param is accepted" << endl;
      QDP_abort(1);
    }

    int NumberOfqqq(0) ;
    switch( params.named_obj.prop_ids.size()){
    case 1: //only up down
      if (params.param.Pt_snk)
	NumberOfqqq++;
      if (params.param.Sl_snk) 
	NumberOfqqq++;
      break ;
    case 2: //only up down (0) and strange (1)
      // 0 ---> up down only: Proton// neutron (up-down degenerate)
      // 1 ---> up down strange: Lambda 
      // 2 3 repeat for smeared
      if (params.param.Pt_snk)
	NumberOfqqq+=2;
      if (params.param.Sl_snk) 
	NumberOfqqq+=2;
      break ;
    default:
      QDPIO::cerr << "OOOPS!! Don't know what to do with all theses propagators.... " << endl;
      throw;
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
    
    
    multi1d<ChromaProp_t> prop_header(params.named_obj.prop_ids.size());
    multi1d<PropSourceConst_t> source_header(params.named_obj.prop_ids.size());
    multi1d<LatticePropagator> qprop(params.named_obj.prop_ids.size());
    multi1d<Real> Mass(params.named_obj.prop_ids.size());

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
	    read(prop_xml, "/Propagator/ForwardProp", prop_header[loop]);
	    read(prop_xml, "/Propagator/PropSource", source_header[loop]);
	  }
	catch (const string& e) 
	  {
	    QDPIO::cerr << "Error extracting forward_prop header: " << e << endl;
	    throw;
	  }
	
	// Hunt around to find the mass and the boundary conditions
	// NOTE: this may be problematic in the future if actions are used with no
	// clear def. of a Mass
	std::istringstream  xml_s(prop_header[loop].fermact);
	XMLReader  fermacttop(xml_s);
	const string fermact_path = "/FermionAction";
	string fermact;
	
	QDPIO::cout << "Try action and mass" << endl;
	try
	  {
	    XMLReader top(fermacttop, fermact_path);
	    
	    read(top, "FermAct", fermact);
	    
	    // Yuk - need to hop some hoops. This should be isolated.
	    if (top.count("Mass") != 0) 
	      {
		read(top, "Mass", Mass[loop]);
	      }
	    else if (top.count("Kappa") != 0)
	      {
		Real Kappa;
		read(top, "Kappa", Kappa);
		Mass[loop] = kappaToMass(Kappa);    // Convert Kappa to Mass
	      }
	    else if (top.count("m_q") != 0) 
	      {
		read(top, "m_q", Mass[loop]);
	      }
	    else
	      {
		QDPIO::cerr << "Neither Mass nor Kappa found" << endl;
		throw "Neither Mass nor Kappa found";
	      }
	    bc[loop] = getFermActBoundary(prop_header[loop].fermact);
	  }
	catch (const string& e) 
	  {
	    QDPIO::cerr << "Error reading fermact or mass or boundary: " << e << endl;
	    QDP_abort(1);
	  }
	
	QDPIO::cout << "FermAct = " << fermact << endl;
	QDPIO::cout << "Mass = " << Mass[loop] << endl;
	QDPIO::cout << "boundary = "
		    << bc[loop][0]<<" "
		    << bc[loop][1]<<" "
		    << bc[loop][2]<<" "
		    << bc[loop][3]<< endl;
      }

    // Derived from input prop
    int j_decay = source_header[0].j_decay;
    //multi1d<int> t_source = source_header[0].t_source;
    //int t0     = t_source[j_decay];
    int t_source = source_header[0].t_source;
    int t0 = t_source ;
    int bc_spec = bc[0][j_decay] ;
    for (int loop(0); loop < params.named_obj.prop_ids.size(); ++loop){
      if(source_header[loop].j_decay!=j_decay){
	QDPIO::cerr << "Error!! j_decay must be the same for all propagators " << endl;
	throw;
      }
      if(bc[loop][j_decay]!=bc_spec){
	QDPIO::cerr << "Error!! bc must be the same for all propagators " << endl;
	throw;
      }
      for(int d(0);d<Nd;d++)
	if(source_header[loop].t_source!=t_source){
	  QDPIO::cerr << "Error!! t_source must be the same for all propagators " << endl;
	  throw;
	}
    }
  

    // phases with  momenta
    SftMom phases(params.param.max_p2, false, j_decay);
    
    
    push(xml_out,"Propagator_info") ;
    write(xml_out, "Masses", Mass);
    write(xml_out, "t_source", t_source);
    for (int loop=0; loop < params.named_obj.prop_ids.size(); ++loop){
      push(xml_out, "Propagator");
      write(xml_out, "ForwardProp",  prop_header[loop]);
      write(xml_out, "PropSource", source_header[loop]);
      multi1d<Double> qp_corr= sumMulti(localNorm2(qprop[loop]),phases.getSet());
      push(xml_out, "Qprop_correlator");
      write(xml_out, "qp_corr", qp_corr);
      pop(xml_out); //Qprop_correlator
      pop(xml_out); //Propagator
    }
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
    for (int loop=0; loop < params.named_obj.prop_ids.size(); ++loop){
      push(file_xml, "Propagator");
      write(file_xml, "ForwardProp",  prop_header[loop]);
      write(file_xml, "PropSource", source_header[loop]);
      multi1d<Double> qp_corr= sumMulti(localNorm2(qprop[loop]),phases.getSet());
      push(file_xml, "Qprop_correlator");
      write(file_xml, "qp_corr", qp_corr);
      pop(file_xml); //Qprop_correlator
      pop(file_xml); //Propagator
    }
    pop(file_xml); //Propagator_info
    write(file_xml, "MomNum", phases.numMom());
    pop(file_xml); //qqqNucNuc_w

    // Write the scalar data
    QDPFileWriter qqqto(file_xml, params.qqq_file, 
			QDPIO_SINGLEFILE, QDPIO_SERIAL, QDPIO_OPEN);
    QDPFileWriter qqbarto(file_xml, params.qqbar_file, 
			  QDPIO_SINGLEFILE, QDPIO_SERIAL, QDPIO_OPEN);
  

    multi2d<ThreeQuarks> qqq(phases.numMom(),phases.numSubsets()); 
    multi2d<DPropagator> qqbar(phases.numMom(),phases.numSubsets()); 
  
    if (params.param.Pt_snk){
      compute_qqq(qqq, qprop[0],qprop[0],qprop[0],phases,t0, bc_spec);
      write_qqq(qqqto, qqq, phases, "nucleon","POINT");

      compute_qqbar(qqbar, qprop[0],qprop[0],phases,t0 );
      write_qqbar(qqbarto, qqbar, phases, "pion","POINT");

      if(params.named_obj.prop_ids.size()==2){
	compute_qqq(qqq, qprop[0],qprop[0],qprop[1],phases,t0, bc_spec);
	write_qqq(qqqto, qqq, phases, "lambda","POINT");
	compute_qqq(qqq, qprop[1],qprop[0],qprop[0],phases,t0, bc_spec);
	write_qqq(qqqto, qqq, phases, "sigma","POINT");
	compute_qqq(qqq, qprop[0],qprop[1],qprop[1],phases,t0, bc_spec);
	write_qqq(qqqto, qqq, phases, "xi","POINT");

	compute_qqbar(qqbar, qprop[0],qprop[1],phases,t0 );
	write_qqbar(qqbarto, qqbar, phases, "kaon","POINT");
	compute_qqbar(qqbar, qprop[1],qprop[0],phases,t0 );
	write_qqbar(qqbarto, qqbar, phases, "kaonbar","POINT");
      }
      
    } 


    // Convolute the quark propagator with the sink smearing function.
    // Make a copy of the quark propagator and then overwrite it with
    // the convolution. 
    if (params.param.Sl_snk){
      multi1d<LatticePropagator> qprop_smr(params.named_obj.prop_ids.size());
      for (int loop=0; loop < params.named_obj.prop_ids.size(); ++loop){
	qprop_smr[loop] = qprop[loop];
	sink_smear2(u, qprop_smr[loop], 
		    params.param.wvf_kind, 
		    params.param.wvf_param[0],//we only use accept one wvf_param
		    params.param.wvfIntPar[0],//and one wvfIntPar. This can change
		    j_decay);
      }   
      compute_qqq(qqq,qprop_smr[0],qprop_smr[0],qprop_smr[0],phases,t0, bc_spec);
      write_qqq(qqqto, qqq, phases, "nucleon","SHELL");

      compute_qqbar(qqbar,qprop_smr[0],qprop_smr[0],phases,t0);
      write_qqbar(qqbarto, qqbar, phases, "pion","SHELL");

      if(params.named_obj.prop_ids.size()==2){
	compute_qqq(qqq,qprop_smr[0],qprop_smr[0],qprop_smr[1],phases,t0,bc_spec);
	write_qqq(qqqto, qqq, phases, "lambda","SHELL");
	compute_qqq(qqq,qprop_smr[1],qprop_smr[0],qprop_smr[0],phases,t0,bc_spec);
	write_qqq(qqqto, qqq, phases, "sigma","SHELL");
	compute_qqq(qqq,qprop_smr[0],qprop_smr[1],qprop_smr[1],phases,t0,bc_spec);
	write_qqq(qqqto, qqq, phases, "xi","SHELL");

	compute_qqbar(qqbar,qprop_smr[0],qprop_smr[1],phases,t0);
	write_qqbar(qqbarto, qqbar, phases, "kaon","SHELL");

	compute_qqbar(qqbar,qprop_smr[1],qprop_smr[0],phases,t0);
	write_qqbar(qqbarto, qqbar, phases, "kaonbar","SHELL");

      }
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
  
};
