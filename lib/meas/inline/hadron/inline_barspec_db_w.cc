/*! \file
 * \brief Inline construction of hadron spectrum
 *
 * Spectrum calculations
 */


#include <sstream>

#include "util/ferm/key_val_db.h"
#include "util/ferm/key_hadron_2pt_corr.h"

#include "inline_barspec_db_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "io/param_io.h"
#include "io/qprop_io.h"
#include "meas/hadron/mesons2_w.h"
#include "meas/hadron/barhqlq_w.h"
#include "meas/hadron/curcor2_w.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/smear/no_quark_displacement.h"

namespace Chroma 
{ 
  namespace InlineBarSpecEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMeas(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;

      const std::string name = "BARYON_DB_SPECTRUM";
    }

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



    void read(XMLReader& xml, const string& path, Params::SpinTerms_t& ter){
      XMLReader paramtop(xml, path);    
      read(paramtop, "spin"   ,ter.spin   );
      read(paramtop, "weight" ,ter.weight );
    }

    void write(XMLWriter& xml, const string& path, const Params::SpinTerms_t& ter){
      push(xml, path);    
      write(xml, "spin"   ,ter.spin   );
      write(xml, "weight" ,ter.weight );
      pop(xml);
    }

    void read(XMLReader& xml, const string& path, Params::SpinWF_t& spWF){
      XMLReader paramtop(xml, path);
      read(paramtop, "norm"     ,spWF.norm   );
      read(paramtop, "terms"    ,spWF.terms  );
    }


    void write(XMLWriter& xml,const string& path, const Params::SpinWF_t& spWF){
      push(xml, path);
      write(xml, "norm"     ,spWF.norm   );
      write(xml, "terms"    ,spWF.terms  );
      pop(xml);
    }

  
    void read(XMLReader& xml, const string& path, Params::Operators_t& op){
      XMLReader paramtop(xml, path);
      //QDPIO::cout<<"Reading Operators: "<<path<<endl ;
      //xml.print(std::cout);
      read(paramtop, "name"     ,op.name   );
      read(paramtop, "spinWF"   ,op.spinWF );
    }

    void write(XMLWriter& xml, const string& path, const Params::Operators_t& op){
      push(xml, path);
      write(xml, "name"     ,op.name   );
      write(xml, "spinWF"   ,op.spinWF );
      pop(xml);
    }

    void read(XMLReader& xml,const string& path, Params::State_t& s){
      XMLReader paramtop(xml, path);
      //QDPIO::cout<<"Reading Operators: "<<path<<endl ;
      //xml.print(std::cout);
      read(paramtop, "name"      ,s.name   );
      read(paramtop, "flavor"    ,s.flavor );
      read(paramtop, "db"        ,s.db     );
      read(paramtop, "spin"      ,s.spin   );
      read(paramtop, "Operators" ,s.ops    );
    }

    void write(XMLWriter& xml, const string& path, const Params::State_t& s){
      push(xml, path);
      write(xml, "name"      ,s.name   );
      write(xml, "flavor"    ,s.flavor );
      write(xml, "db"        ,s.db     );
      write(xml, "spin"      ,s.spin   );
      write(xml, "Operators" ,s.ops    );
      pop(xml);
    }



    //! Reader for parameters
    void read(XMLReader& xml, const string& path, Params::Param_t& param)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
	break;

      default:
	QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
	QDP_abort(1);
      }

      //read(paramtop, "MesonP", param.MesonP);
      //read(paramtop, "CurrentP", param.CurrentP);
      //read(paramtop, "BaryonP", param.BaryonP);

      read(paramtop, "time_rev", param.time_rev);
      read(paramtop, "mom2_max", param.mom2_max);
      read(paramtop, "avg_equiv_mom", param.avg_equiv_mom);
      read(paramtop, "ensemble", param.ensemble);
      read(paramtop, "States", param.states);
    }


    //! Writer for parameters
    void write(XMLWriter& xml, const string& path, const Params::Param_t& param)
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
      write(xml, "ensemble", param.ensemble);
      write(xml, "States", param.states);

      pop(xml);
    }

  

    //! Propagator input
    void read(XMLReader& xml, const string& path, Params::NamedObject_t::Props_t& input)
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
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const Params::NamedObject_t::Props_t& input)
    {
      push(xml, path);

      write(xml, "up_id", input.up_id);
      write(xml, "down_id", input.down_id);
      write(xml, "strange_id", input.strange_id);
      write(xml, "charm_id", input.charm_id);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id" , input.gauge_id);
      read(inputtop, "props"    , input.props);
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);
      write(xml, "gauge_id" , input.gauge_id);
      write(xml, "props"    , input.props);
      pop(xml);
    }


    // Param stuff
    Params::Params()
    { 
      frequency = 0; 
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
	QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << endl;
	QDP_abort(1);
      }
    }


    void
    Params::writeXML(XMLWriter& xml_out, const std::string& path) 
    {
      push(xml_out, path);
    
      write(xml_out, "Param", param);
      write(xml_out, "NamedObject", named_obj);
      write(xml_out, "xml_file", xml_file);

      pop(xml_out);
    }


    void epsilon_contract(LatticeComplex& res,
			  const multi2d<LatticeComplex>& l,
			  const multi2d<LatticeComplex>& m,
			  const multi2d<LatticeComplex>& r){
      //LatticeComplex res = 0.0;
    
      //generated by : make_epsilon_contr.pl | sort                

      res  = l(0,0)*m(1,1)*r(2,2);
      res -= l(0,0)*m(1,2)*r(2,1);
      res -= l(0,0)*m(2,1)*r(1,2);
      res += l(0,0)*m(2,2)*r(1,1);
      res -= l(0,1)*m(1,0)*r(2,2);
      res += l(0,1)*m(1,2)*r(2,0);
      res += l(0,1)*m(2,0)*r(1,2);
      res -= l(0,1)*m(2,2)*r(1,0);
      res += l(0,2)*m(1,0)*r(2,1);
      res -= l(0,2)*m(1,1)*r(2,0);
      res -= l(0,2)*m(2,0)*r(1,1);
      res += l(0,2)*m(2,1)*r(1,0);
      res -= l(1,0)*m(0,1)*r(2,2);
      res += l(1,0)*m(0,2)*r(2,1);
      res += l(1,0)*m(2,1)*r(0,2);
      res -= l(1,0)*m(2,2)*r(0,1);
      res += l(1,1)*m(0,0)*r(2,2);
      res -= l(1,1)*m(0,2)*r(2,0);
      res -= l(1,1)*m(2,0)*r(0,2);
      res += l(1,1)*m(2,2)*r(0,0);
      res -= l(1,2)*m(0,0)*r(2,1);
      res += l(1,2)*m(0,1)*r(2,0);
      res += l(1,2)*m(2,0)*r(0,1);
      res -= l(1,2)*m(2,1)*r(0,0);
      res += l(2,0)*m(0,1)*r(1,2);
      res -= l(2,0)*m(0,2)*r(1,1);
      res -= l(2,0)*m(1,1)*r(0,2);
      res += l(2,0)*m(1,2)*r(0,1);
      res -= l(2,1)*m(0,0)*r(1,2);
      res += l(2,1)*m(0,2)*r(1,0);
      res += l(2,1)*m(1,0)*r(0,2);
      res -= l(2,1)*m(1,2)*r(0,0);
      res += l(2,2)*m(0,0)*r(1,1);
      res -= l(2,2)*m(0,1)*r(1,0);
      res -= l(2,2)*m(1,0)*r(0,1);
      res += l(2,2)*m(1,1)*r(0,0);

      //return res ;
    }

    // Anonymous namespace
    namespace 
    {
      //! Useful structure holding sink props
      class SinkPropContainer_t{
      public:
	ForwardProp_t prop_header;
	string quark_propagator_id;
	Real Mass;

	bool exists ;

	multi1d<int> bc; 
    
	string source_type;
	string source_disp_type;
	string sink_type;
	string sink_disp_type;

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
		string xpath;
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
	      QDPIO::cerr << name << ": caught dynamic cast error" 
			  << endl;
	      QDP_abort(1);
	    }
	    catch (const string& e) 
	    {
	      QDPIO::cerr << name << ": error message: " << e 
			  << endl;
	      QDP_abort(1);
	    }
	
	
	    // Derived from input prop
	    // Hunt around to find the mass
	    // NOTE: this may be problematic in the future if actions are used with no
	    // clear def. of a Mass
	    QDPIO::cout << "Try action and mass" << endl;
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
	    catch (const string& e) 
	    {
	      QDPIO::cerr << name 
			  << ": caught exception. No BC found in these headers. Will assume dirichlet: " << e 
			  << endl;
	    }
	
	    QDPIO::cout << "FermAct = " << prop_header.prop_header.fermact.id << endl;
	    QDPIO::cout << "Mass = " << Mass << endl;
	  }
      
      };





      void sanity_check_props(const SinkPropContainer_t& p1, 
			      const SinkPropContainer_t& p2){
      
	int j_decay = p1.prop_header.source_header.j_decay;
	if (p2.prop_header.source_header.j_decay != 
	    p1.prop_header.source_header.j_decay){
	  QDPIO::cerr << "Error!! j_decay must be the same for all propagators " << endl;
	  QDP_abort(1);
	}
	if (p2.prop_header.source_header.t_source != 
	    p1.prop_header.source_header.t_source)
	{
	  QDPIO::cerr << "Error!! t_source must be the same for all propagators " << endl;
	  QDP_abort(1);
	}
	if (p1.source_type != p2.source_type)
	{
	  QDPIO::cerr << "Error!! source_type must be the same in a pair " << endl;
	  QDP_abort(1);
	}
	if (p1.sink_type != p2.sink_type)
	{
	  QDPIO::cerr << "Error!! source_type must be the same in a pair " << endl;
	  QDP_abort(1);
	}
      
	if (p2.bc[j_decay] != p1.bc[j_decay])
	{
	  QDPIO::cerr << "Error!! bc must be the same for all propagators " << endl;
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
	map<string,SinkPropContainer_t>  prop;
	map<string,BarSpec::RPropagator> rprop;

	//! Read all sinks
	AllSinkProps_t(const Params::NamedObject_t::Props_t& p){

	  QDPIO::cout<<"Attempt to parse forward propagator= "<<p.up_id<<endl;
	  prop["up"].readSinkProp(p.up_id);
	  rprop[p.up_id].ConvertProp(TheNamedObjMap::Instance().getData<LatticePropagator>(p.up_id));
	  QDPIO::cout<<"up quark  propagator successfully parsed" << endl;
	  j_decay = prop["up"].prop_header.source_header.j_decay;
	  t0      = prop["up"].prop_header.source_header.t_source;
	  t_srce  = prop["up"].prop_header.source_header.getTSrce();
	  bc_spec = prop["up"].bc[j_decay];

	  QDPIO::cout<<"Attempt to parse forward propagator= " <<p.down_id<<endl;
	  //Always need a down quark 
	  prop["down"].readSinkProp(p.down_id);
	  QDPIO::cout << "down quark propagator successfully parsed" << endl;
	  if(rprop.find(p.down_id) == rprop.end()){
	    QDPIO::cout<<__func__<<": Need to convert prop id: "
		       <<p.down_id<<endl;
	    rprop[p.down_id].ConvertProp(TheNamedObjMap::Instance().getData<LatticePropagator>(p.down_id));
	  }
	  QDPIO::cout<<"Attempt to parse forward propagator= "<<p.strange_id<<endl;
	  prop["strange"].readSinkProp(p.strange_id);
	  if(p.strange_id != "NULL"){
	    QDPIO::cout <<"strange quark propagator successfully parsed" << endl;
	    if(rprop.find(p.strange_id) == rprop.end()){
	      QDPIO::cout<<__func__<<": Need to convert prop id: "
			 <<p.strange_id<<endl;
	      rprop[p.strange_id].ConvertProp(TheNamedObjMap::Instance().getData<LatticePropagator>(p.strange_id));
	    }
	  }

	  QDPIO::cout<<"Attempt to parse forward propagator= "<<p.charm_id<<endl;
	  prop["charm"].readSinkProp(p.charm_id);
	  if(p.charm_id != "NULL"){
	    QDPIO::cout << "charm quark propagator successfully parsed" << endl;
	    if(rprop.find(p.charm_id) == rprop.end()){
	      QDPIO::cout<<__func__<<": Need to convert prop id: "
			 <<p.charm_id<<endl;
	      rprop[p.charm_id].ConvertProp(TheNamedObjMap::Instance().getData<LatticePropagator>(p.charm_id));
	    }
	  }
	
	 
	  sanity_check_props(prop["up"], prop["down"]) ;
	  if(prop["strange"].exists)
	    sanity_check_props(prop["up"], prop["strange"]) ;
	  if(prop["charm"].exists)
	    sanity_check_props(prop["up"], prop["charm"]) ;

	}


	string sink(const string& flavor){
	  return prop[flavor].sink_type ;
	}

	string source(const string& flavor){
	  return prop[flavor].source_type ;
	}

	const BarSpec::RPropagator& prop_ref(const string& flavor){

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
	  
	    QDPIO::cout << "propagator_id = " << p.quark_propagator_id << endl;
	    QDPIO::cout << "Source_type   = " << p.source_type << endl;
	    QDPIO::cout << "Sink_type     = " << p.sink_type << endl;
	  
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

    namespace BarSpec{
    
      vector<int> permutation(int k, const vector<int>& s){
	vector<int> p = s ;
	int ss ;
	int jj ;
	for(int j(2); j<=s.size(); j++){
	  jj=k%j ;
	  ss=p[jj]; p[jj]=p[j-1]; p[j-1]=ss;
	  k=k/j ;
	}
	return p ;
      }
    
      int permutation_sign(int k, int n){
	int jj ;
	int sign_(1) ;
	for(int j(2); j<=n; j++){
	  jj=k%j ;
	  if(jj!=(j-1))
	    sign_ *= (-1) ;
	  k=k/j ;
	}
	return sign_ ;
      }

      void contract(LatticeComplex& latC, 
		    const RPropagator& q1,
		    const RPropagator& q2,
		    const RPropagator& q3,
		    const SpinWF_t& snk,
		    const SpinWF_t& src){
	LatticeComplex cc ;
	latC = 0.0;
	for(int s(0);s<src.terms.size();s++)
	  for(int ss(0);ss<snk.terms.size();ss++){
	    epsilon_contract(cc,
			     q1.p(snk.terms[ss].spin[0],src.terms[s].spin[0]),
			     q2.p(snk.terms[ss].spin[1],src.terms[s].spin[1]),
			     q3.p(snk.terms[ss].spin[2],src.terms[s].spin[2]) );
	    latC += (snk.terms[ss].weight*src.terms[s].weight)*cc ;
	  }// loop over source sink wavefunction components
	latC *= (snk.norm*src.norm) ;
      }
    }//barspec name space


    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "barspec");
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
	QDPIO::cerr << name << ": caught dynamic cast error" 
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

      push(xml_out, "barspec");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << " BARSPEC: Spectroscopy for Wilson-like fermions" << endl;
      QDPIO::cout << endl << "     Gauge group: SU(" << Nc << ")" << endl;
      QDPIO::cout << "     volume: " << Layout::lattSize()[0];
      for (int i=1; i<Nd; ++i) {
	QDPIO::cout << " x " << Layout::lattSize()[i];
      }
      QDPIO::cout << endl;

      proginfo(xml_out);    // Print out basic program info

      // Write out the input
      params.writeXML(xml_out, "Input");

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
	const Params::NamedObject_t::Props_t named_obj = params.named_obj.props[p];

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

	pop(xml_out);

	push(xml_out, "Forward_Propagator_Properties"                     );
	map<string,SinkPropContainer_t>::iterator it;
	for(it=all_sinks.prop.begin();it!=all_sinks.prop.end();it++)
	  xml_print_prop_info(xml_out, it->second      , phases.getSet()  );
      
	pop(xml_out);

	int Nt = Layout::lattSize()[j_decay];

	LatticeComplex latC ;
	StopWatch tictoc;
	tictoc.reset();
	tictoc.start();
	for(int s(0);s<params.param.states.size();s++){
	  QDPIO::cout<<"Doing state "<<params.param.states[s].name<<endl;
	  QDPIO::cout<<"   Flavor structure: " ;
	  for(int k(0);k<params.param.states[s].flavor.size();k++)
	    QDPIO::cout<<"  "<<params.param.states[s].flavor[k] ;
	  QDPIO::cout<<endl;
	  multi1d<string> prop_id(params.param.states[s].flavor.size()) ;
	  for(int k(0);k<params.param.states[s].flavor.size();k++){
	    prop_id[k] = params.param.states[s].flavor[k] ;
	  }
	  // Open the file, and write the meta-data and the binary for this state
	  if (! qdp_db.fileExists(params.param.states[s].db)){
	    XMLBufferWriter file_xml;
	  
	    push(file_xml, "DBMetaData");
	    write(file_xml, "id", string("hadron2Pt"));
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
	  const BarSpec::RPropagator& q1 = all_sinks.prop_ref(prop_id[0]) ;
	  const BarSpec::RPropagator& q2 = all_sinks.prop_ref(prop_id[1]) ;
	  const BarSpec::RPropagator& q3 = all_sinks.prop_ref(prop_id[2]) ;

	  KeyHadron2PtCorr_t key  ;

	  ostringstream os ;
	  os<<all_sinks.prop[prop_id[0]].Mass<<" ";
	  os<<all_sinks.prop[prop_id[1]].Mass<<" ";
	  os<<all_sinks.prop[prop_id[2]].Mass ;
	  key.mass        = os.str();
	  key.ensemble    = params.param.ensemble;
	  key.num_vecs = Nc;  /*!< Number of color vectors */

	  key.src_smear    = all_sinks.source(prop_id[0])+" ";
	  key.src_smear   += all_sinks.source(prop_id[1])+" ";
	  key.src_smear   += all_sinks.source(prop_id[2]) ;	    
	  key.src_spin  =   params.param.states[s].spin ;
	  key.src_lorentz.resize(0);
	  //key.src_lorentz =  key.src_spin  ;

	  key.snk_smear    = all_sinks.sink(prop_id[0])+" ";
	  key.snk_smear   += all_sinks.sink(prop_id[1])+" ";
	  key.snk_smear   += all_sinks.sink(prop_id[2]) ;
	  key.snk_spin  =  params.param.states[s].spin ;
	  key.snk_lorentz.resize(0);
	  //key.snk_lorentz =  key.snk_spin  ;

	  //loop over momenta goes here
	  for(int oi(0);oi<params.param.states[s].ops.size();oi++){ //sink
	    BarSpec::SpinWF_t snk(params.param.states[s].ops[oi].spinWF);
	    snk.permutations(prop_id);
	    for(int oj(0);oj<params.param.states[s].ops.size();oj++){//source

	      BarSpec::SpinWF_t src(params.param.states[s].ops[oj].spinWF);

	      key.src_name    = params.param.states[s].ops[oj].name;
	      key.snk_name    = params.param.states[s].ops[oi].name;
	    
	      BarSpec::contract(latC,q1,q2,q3,snk,src);

	      multi2d<ComplexD> hsum;
	      hsum = phases.sft(latC) ;
	    
	      for(int mom(0);mom<phases.numMom();mom++){
		key.mom = phases.numToMom(mom);    /*<! Momentum  */
		SerialDBKey<KeyHadron2PtCorr_t> K;
		K.key() = key ;
		SerialDBData<multi1d<ComplexD> >  V ;
		V.data().resize(Nt);
		for(int t(0);t<Nt;t++){
		  int t_eff = (t - t0 + Nt) % Nt;
		  if ( bc_spec < 0 && (t_eff+t0) >= Nt)
		    V.data()[t_eff] = -hsum[mom][t];
		  else
		    V.data()[t_eff] =  hsum[mom][t];
		}//loop over time
		qdp_db.insert(K,V);
	      }// loop over momenta

	    }// loop  over source ops
	  }// loop over sink ops
	  qdp_db.close() ;
	}// loop over states
	tictoc.stop();
	QDPIO::cout << name << ": contraction time = "
		    << tictoc.getTimeInSeconds() 
		    << " secs" << endl;

	pop(xml_out); //elem
      }// loop over propagator groups
      pop(xml_out);  // barspec

      snoop.stop();
      QDPIO::cout << name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << name << ": ran successfully" << endl;

      END_CODE();
    } 

  }

}
