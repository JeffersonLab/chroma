// $Id: inline_grid_prop_w.cc,v 3.3 2008-12-12 04:13:34 kostas Exp $
/*! \file
 * \brief Compute the matrix element of   M^-1 * multi1d<LatticeColorVector>
 *
 * Propagator calculation on a colorvector
 */

#include "fermact.h"
#include "meas/inline/hadron/inline_grid_prop_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ferm/subset_vectors.h"
#include "util/ferm/map_obj.h"
#include "util/ferm/map_obj_memory.h"
#include "util/ferm/key_grid_prop.h"
#include "util/ferm/transf.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/inline/make_xml_file.h"

#include "meas/sources/diluteGrid_source_const.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlineGridPropEnv 
  {
    //! Propagator input
    void read(XMLReader& xml, const string& path, InlineGridPropEnv::Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "prop_id", input.prop_id);
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const InlineGridPropEnv::Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "prop_id", input.prop_id);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const string& path, InlineGridPropEnv::Params::Param_t::Sources_t& input)
    {
      XMLReader inputtop(xml, path);
      read(inputtop, "spatial_mask_size",input.spatial_mask_size);
      read(inputtop, "spatial_masks", input.spatial_masks);

      read(inputtop, "t_sources", input.t_sources);
      read(inputtop, "decay_dir", input.decay_dir);
      
      input.smear = false ;
      if(inputtop.count("Smearing") !=0 ) {
	input.smr = readXMLGroup(inputtop, "Smearing", "wvf_kind");
	input.link_smear = readXMLGroup(inputtop, "LinkSmearing", "LinkSmearingType");
	input.displace = readXMLGroup(inputtop, "Displacement","DisplacementType");
	input.smear = true ;
      }
      
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const InlineGridPropEnv::Params::Param_t::Sources_t& out)
    {
      push(xml, path);

      write(xml, "spatial_mask_size", out.spatial_mask_size);
      write(xml, "spatial_masks", out.spatial_masks);

      write(xml, "t_sources", out.t_sources);
      write(xml, "decay_dir", out.decay_dir);

      if(out.smear){
	 xml << out.smr.xml;
	 xml << out.displace.xml ;
	 xml << out.link_smear.xml ;
      }

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const string& path, InlineGridPropEnv::Params::Param_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "Propagator", input.prop) ;
      read(inputtop, "Sources"   , input.src)  ;
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const InlineGridPropEnv::Params::Param_t& input)
    {
      push(xml, path);

      write(xml, "Propagator", input.prop) ;
      write(xml, "Sources"   , input.src)  ;

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const string& path, InlineGridPropEnv::Params& input)
    {
      InlineGridPropEnv::Params tmp(xml, path);
      input = tmp;
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const InlineGridPropEnv::Params& input)
    {
      push(xml, path);
    
      write(xml, "Param", input.param);
      write(xml, "NamedObject", input.named_obj);

      pop(xml);
    }
  } // namespace InlineGridPropEnv 


  namespace InlineGridPropEnv 
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
    }
      
    const std::string name = "GRID_PROPAGATOR";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= WilsonTypeFermActsEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }


    //----------------------------------------------------------------------------
    // Param stuff
    Params::Params() { frequency = 0; }

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



    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "GridProp");
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
      multi1d<LatticeColorMatrix> u;
      XMLBufferWriter gauge_xml;
      try
      {
	u = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
	TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << name << ": caught dynamic cast error" << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << name << ": map call failed: " << e << endl;
	QDP_abort(1);
      }

      push(xml_out, "GridProp");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": propagator calculation" << endl;

      proginfo(xml_out);    // Print out basic program info

      // Write out the input
      write(xml_out, "Input", params);

      // Write out the config header
      write(xml_out, "Config_info", gauge_xml);

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      // Calculate some gauge invariant observables just for info.
      MesPlq(xml_out, "Observables", u);


      //
      // Create the output files
      //
      try
	{	
	// Create the object as a handle. 
	// This bit will and up changing to a Factory invocation
	Handle< MapObject<KeyGridProp_t, LatticeFermion> > new_map_obj_handle(new MapObjectMemory<KeyGridProp_t, LatticeFermion>() );

	// Create the entry
	TheNamedObjMap::Instance().create< Handle< MapObject<KeyGridProp_t,LatticeFermion> > >(params.named_obj.prop_id);

	// Insert
	TheNamedObjMap::Instance().getData< Handle<MapObject<KeyGridProp_t,LatticeFermion> > >(params.named_obj.prop_id) = new_map_obj_handle;

	}
      catch (std::bad_cast)
	{
	  QDPIO::cerr << name << ": caught dynamic cast error" << endl;
	  QDP_abort(1);
	}
      catch (const string& e) 
	{
	  QDPIO::cerr << name << ": error creating prop: " << e << endl;
	  QDP_abort(1);
	}

      // Cast should be valid now
      MapObject<KeyGridProp_t,LatticeFermion>& map_obj =
	*(TheNamedObjMap::Instance().getData< Handle<MapObject<KeyGridProp_t,LatticeFermion> > >(params.named_obj.prop_id));



      DiluteGridQuarkSourceConstEnv::Params srcParams ;
      srcParams.j_decay = params.param.src.decay_dir ;
      srcParams.spatial_mask_size = params.param.src.spatial_mask_size ;
      if(params.param.src.smear){
	srcParams.smear = true ;
	srcParams.smr = params.param.src.smr ;
	srcParams.displace = params.param.src.displace ;
	srcParams.link_smear = params.param.src.link_smear ;
	
      }

      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, params.param.src.decay_dir);

      // Total number of iterations
      int ncg_had = 0;
      try{
	StopWatch swatch;
	swatch.reset();
	QDPIO::cout << "Try the various factories" << endl;
	
	// Typedefs to save typing
	typedef LatticeFermion               T;
	typedef multi1d<LatticeColorMatrix>  P;
	typedef multi1d<LatticeColorMatrix>  Q;

	//
	// Initialize fermion action
	//
	std::istringstream  xml_s(params.param.prop.fermact.xml);
	XMLReader  fermacttop(xml_s);
	QDPIO::cout << "FermAct = " << params.param.prop.fermact.id << endl;

	// Generic Wilson-Type stuff
	Handle< FermionAction<T,P,Q> >
	  S_f(TheFermionActionFactory::Instance().createObject(params.param.prop.fermact.id,
							       fermacttop,
							       params.param.prop.fermact.path));

	Handle< FermState<T,P,Q> > state(S_f->createState(u));

	Handle< SystemSolver<LatticeFermion> > PP = S_f->qprop(state,
							       params.param.prop.invParam);
      
	QDPIO::cout << "Suitable factory found: compute all the quark props" << endl;
	swatch.start();

	for(int t0(0);t0<params.param.src.t_sources.size();t0++){//loop over timeslices
	  QDPIO::cout << name << ": Doing timeslice " ;
	  QDPIO::cout << params.param.src.t_sources[t0] << endl ;
	  srcParams.t_source = params.param.src.t_sources[t0] ;
	  for(int s(0);s<Ns;s++){// loop over spins
	    srcParams.spin = s;
	    for(int c(0);c<Nc;c++){// loop over colors
	      srcParams.color = c;
	      for(int g(0);g<params.param.src.spatial_masks.size();g++){//loop over grids
		srcParams.spatial_mask = params.param.src.spatial_masks[g];
		
		QDPIO::cout << name << ": Doing spin " << s ;
		QDPIO::cout <<" color " << c << " and  grid " << g<< " ( of " ;
		QDPIO::cout << params.param.src.spatial_masks.size() << ")"<<endl;
		
		//Constuct the source
		DiluteGridQuarkSourceConstEnv::SourceConst<LatticeFermion>  GridSrc(srcParams);
		LatticeFermion chi = GridSrc(u);// one can save time here by passed an presmeared field

		LatticeFermion quark_soln = zero;
		// Do the propagator inversion
		SystemSolverResults_t res = (*PP)(quark_soln, chi);
		ncg_had += res.n_count;

		KeyGridProp_t key;
		key.t_source  = srcParams.t_source ;
		key.color     = c ;
		key.spin      = s ;
		key.grid      = g ;
		map_obj.insert(key, quark_soln);
		
		//QDPIO::cout << params.param.src.spatial_masks[g]<<endl ;
	      }//grids
	    }//colors
	  }//spins
	}//t0

	swatch.stop();
	QDPIO::cout << "  Propagators computed: time= " 
		    << swatch.getTimeInSeconds() 
		    << " secs" << endl;
      }
      catch (const std::string& e){
	QDPIO::cout << name << ": caught exception around qprop: " << e << endl;
	QDP_abort(1);
      }

      push(xml_out,"Relaxation_Iterations");
      write(xml_out, "ncg_had", ncg_had);
      pop(xml_out);

      pop(xml_out);  // GridProp

      // Write the meta-data for this operator
      {
	XMLBufferWriter file_xml;

	push(file_xml, "GridProp");
	write(file_xml, "num_records", map_obj.size()); 
	write(file_xml, "Params", params.param);
	write(file_xml, "Config_info", gauge_xml);
	pop(file_xml);

	XMLBufferWriter record_xml;
	push(record_xml, "GridProp");
	write(record_xml, "num_records", map_obj.size()); 
	pop(record_xml);

	// Write the propagator xml info
	TheNamedObjMap::Instance().get(params.named_obj.prop_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.prop_id).setRecordXML(record_xml);
      }

      snoop.stop();
      QDPIO::cout << name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << name << ": ran successfully" << endl;

      END_CODE();
    } 

  }

} // namespace Chroma
