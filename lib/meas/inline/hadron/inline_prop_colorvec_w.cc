/*! \file
 * \brief Compute the matrix element of   M^-1 * multi1d<LatticeColorVector>
 *
 * Propagator calculation on a colorstd::vector
 */

#include "fermact.h"
#include "meas/inline/hadron/inline_prop_colorvec_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ferm/subset_vectors.h"
#include "util/ferm/map_obj/map_obj_aggregate_w.h"
#include "util/ferm/map_obj/map_obj_factory_w.h"
#include "util/ferm/key_prop_colorvec.h"
#include "util/ferm/transf.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlinePropColorVecEnv 
  {
    //! Propagator input
    void read(XMLReader& xml, const std::string& path, InlinePropColorVecEnv::Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "colorvec_id", input.colorvec_id);
      read(inputtop, "prop_id", input.prop_id);

      // User Specified MapObject tags
      input.prop_obj = readXMLGroup(inputtop, "PropMapObject", "MapObjType");
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const InlinePropColorVecEnv::Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "colorvec_id", input.colorvec_id);
      write(xml, "prop_id", input.prop_id);
      xml << input.prop_obj.xml;

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const std::string& path, InlinePropColorVecEnv::Params::Param_t::Contract_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "num_vecs", input.num_vecs);
      read(inputtop, "t_sources", input.t_sources);
      read(inputtop, "decay_dir", input.decay_dir);
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const InlinePropColorVecEnv::Params::Param_t::Contract_t& input)
    {
      push(xml, path);

      write(xml, "num_vecs", input.num_vecs);
      write(xml, "t_sources", input.t_sources);
      write(xml, "decay_dir", input.decay_dir);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const std::string& path, InlinePropColorVecEnv::Params::Param_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "Propagator", input.prop);
      read(inputtop, "Contractions", input.contract);
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const InlinePropColorVecEnv::Params::Param_t& input)
    {
      push(xml, path);

      write(xml, "Propagator", input.prop);
      write(xml, "Contractions", input.contract);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const std::string& path, InlinePropColorVecEnv::Params& input)
    {
      InlinePropColorVecEnv::Params tmp(xml, path);
      input = tmp;
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const InlinePropColorVecEnv::Params& input)
    {
      push(xml, path);
    
      write(xml, "Param", input.param);
      write(xml, "NamedObject", input.named_obj);

      pop(xml);
    }
  } // namespace InlinePropColorVecEnv 


  namespace InlinePropColorVecEnv 
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
      
    const std::string name = "PROP_COLORVEC";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= WilsonTypeFermActsEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	success &= MapObjectWilson4DEnv::registerAll();
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
	QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
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
	std::string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "PropColorVec");
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
	QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
	QDP_abort(1);
      }
      catch (const std::string& e) 
      {
	QDPIO::cerr << name << ": std::map call failed: " << e << std::endl;
	QDP_abort(1);
      }

      push(xml_out, "PropColorVec");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": propagator calculation" << std::endl;

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
      // Read in the source along with relevant information.
      // 
      XMLReader source_file_xml, source_record_xml;

      QDPIO::cout << "Snarf the source from a named buffer" << std::endl;
      try
      {
	TheNamedObjMap::Instance().getData< Handle< QDP::MapObject<int,EVPair<LatticeColorVector> > > >(params.named_obj.colorvec_id);

	// Snarf the source info. This is will throw if the colorvec_id is not there

	TheNamedObjMap::Instance().get(params.named_obj.colorvec_id).getFileXML(source_file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.colorvec_id).getRecordXML(source_record_xml);

	// Write out the source header
	write(xml_out, "Source_file_info", source_file_xml);
	write(xml_out, "Source_record_info", source_record_xml);
      }    
      catch (std::bad_cast) {
	QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
	QDP_abort(1);
      }
      catch (const std::string& e) {
	QDPIO::cerr << name << ": error extracting source_header: " << e << std::endl;
	QDP_abort(1);
      }
      catch( const char* e) {
	QDPIO::cerr << name <<": Caught some char* exception:" << std::endl;
	QDPIO::cerr << e << std::endl;
	QDPIO::cerr << "Rethrowing" << std::endl;
	throw;
      }

      // Cast should be valid now
      const QDP::MapObject<int,EVPair<LatticeColorVector> >& eigen_source = 
	*(TheNamedObjMap::Instance().getData< Handle< QDP::MapObject<int,EVPair<LatticeColorVector> > > >(params.named_obj.colorvec_id));

      QDPIO::cout << "Source successfully read and parsed" << std::endl;

      //
      // Create the output files
      //
      try
      {
	// Generate a metadata
	std::string file_str;
	if (1)
	{
	  XMLBufferWriter file_xml;

	  push(file_xml, "MODMetaData");
	  write(file_xml, "id", std::string("propColorVec"));
	  write(file_xml, "lattSize", QDP::Layout::lattSize());
	  write(file_xml, "decay_dir", params.param.contract.decay_dir);
	  write(file_xml, "num_vecs", params.param.contract.num_vecs);
	  // write(file_xml, "mass_label", params.param.contract.mass_label);  // no simple kind of mass label available. That is unfortunate...
	  proginfo(file_xml);    // Print out basic program info
	  write(file_xml, "Propagator", params.param.prop);
	  write(file_xml, "Contractions", params.param.contract);
	  write(file_xml, "Config_info", gauge_xml);
	  pop(file_xml);

	  file_str = file_xml.str();
	}
	  
	// Create the map object
	std::istringstream  xml_s(params.named_obj.prop_obj.xml);
	XMLReader MapObjReader(xml_s);
	
	// Create the entry
	TheNamedObjMap::Instance().create< Handle< QDP::MapObject<KeyPropColorVec_t,LatticeFermion> > >(params.named_obj.prop_id);
	TheNamedObjMap::Instance().getData< Handle< QDP::MapObject<KeyPropColorVec_t,LatticeFermion> > >(params.named_obj.prop_id) =
	  TheMapObjKeyPropColorVecFactory::Instance().createObject(params.named_obj.prop_obj.id,
								   MapObjReader,
								   params.named_obj.prop_obj.path,
								   file_str);

	TheNamedObjMap::Instance().getData< Handle<QDP::MapObject<KeyPropColorVec_t,LatticeFermion> > >(params.named_obj.prop_id)->insertUserdata(file_str);
      }
      catch (std::bad_cast)
      {
	QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
	QDP_abort(1);
      }
      catch (const std::string& e) 
      {
	QDPIO::cerr << name << ": error creating prop: " << e << std::endl;
	QDP_abort(1);
      }

      // Cast should be valid now
      QDP::MapObject<KeyPropColorVec_t,LatticeFermion>& prop_obj =
	*(TheNamedObjMap::Instance().getData< Handle<QDP::MapObject<KeyPropColorVec_t,LatticeFermion> > >(params.named_obj.prop_id));

      // Sanity check - write out the norm2 of the source in the Nd-1 direction
      // Use this for any possible verification
      {
	// Initialize the slow Fourier transform phases
	SftMom phases(0, true, Nd-1);

	EVPair<LatticeColorVector> tmpvec; eigen_source.get(0, tmpvec);
	multi1d<Double> source_corrs = sumMulti(localNorm2(tmpvec.eigenVector), phases.getSet());

	push(xml_out, "Source_correlators");
	write(xml_out, "source_corrs", source_corrs);
	pop(xml_out);
      }

      // Another sanity check
      if (params.param.contract.num_vecs > eigen_source.size())
      {
	QDPIO::cerr << __func__ << ": num_vecs= " << params.param.contract.num_vecs
		    << " is greater than the number of available colorvectors= "
		    << eigen_source.size() << std::endl;
	QDP_abort(1);
      }


      // Total number of iterations
      int ncg_had = 0;

      //
      // Try the factories
      //
      try
      {
	StopWatch swatch;
	swatch.reset();
	QDPIO::cout << "Try the various factories" << std::endl;

	// Typedefs to save typing
	typedef LatticeFermion               T;
	typedef multi1d<LatticeColorMatrix>  P;
	typedef multi1d<LatticeColorMatrix>  Q;

	//
	// Initialize fermion action
	//
	std::istringstream  xml_s(params.param.prop.fermact.xml);
	XMLReader  fermacttop(xml_s);
	QDPIO::cout << "FermAct = " << params.param.prop.fermact.id << std::endl;

	// Generic Wilson-Type stuff
	Handle< FermionAction<T,P,Q> >
	  S_f(TheFermionActionFactory::Instance().createObject(params.param.prop.fermact.id,
							       fermacttop,
							       params.param.prop.fermact.path));

	Handle< FermState<T,P,Q> > state(S_f->createState(u));

	Handle< SystemSolver<LatticeFermion> > PP = S_f->qprop(state,
							       params.param.prop.invParam);
      
	QDPIO::cout << "Suitable factory found: compute all the quark props" << std::endl;
	swatch.start();

	//
	// Loop over the source color and spin, creating the source
	// and calling the relevant propagator routines.
	//
	const int num_vecs            = params.param.contract.num_vecs;
	const int decay_dir           = params.param.contract.decay_dir;
	const multi1d<int>& t_sources = params.param.contract.t_sources;

	// Initialize the slow Fourier transform phases
	SftMom phases(0, true, decay_dir);


	// Loop over each operator 
	for(int tt=0; tt < t_sources.size(); ++tt)
	{
	  int t_source = t_sources[tt];
	  QDPIO::cout << "t_source = " << t_source << std::endl; 

	  // All the loops
	  for(int colorvec_source=0; colorvec_source < num_vecs; ++colorvec_source)
	  {
	    QDPIO::cout << "colorvec_source = " << colorvec_source << std::endl; 

	    // Pull out a time-slice of the color std::vector source
	    LatticeColorVector vec_srce = zero;
	    EVPair<LatticeColorVector> tmpvec;
	    eigen_source.get(colorvec_source, tmpvec);
	    vec_srce[phases.getSet()[t_source]] = tmpvec.eigenVector;
	
	    for(int spin_source=0; spin_source < Ns; ++spin_source)
	    {
	      QDPIO::cout << "spin_source = " << spin_source << std::endl; 

	      // Insert a ColorVector into spin index spin_source
	      // This only overwrites sections, so need to initialize first
	      LatticeFermion chi = zero;
	      CvToFerm(vec_srce, chi, spin_source);

	      LatticeFermion quark_soln = zero;

	      // Do the propagator inversion
	      SystemSolverResults_t res = (*PP)(quark_soln, chi);
	      ncg_had = res.n_count;

	      KeyPropColorVec_t key;
	      key.t_source     = t_source;
	      key.colorvec_src = colorvec_source;
	      key.spin_src     = spin_source;
		  
	      prop_obj.insert(key, quark_soln);
	    } // for spin_source
	  } // for colorvec_source
	} // for t_source

	swatch.stop();
	QDPIO::cout << "Propagators computed: time= " 
		    << swatch.getTimeInSeconds() 
		    << " secs" << std::endl;
      }
      catch (const std::string& e) 
      {
	QDPIO::cout << name << ": caught exception around qprop: " << e << std::endl;
	QDP_abort(1);
      }

      push(xml_out,"Relaxation_Iterations");
      write(xml_out, "ncg_had", ncg_had);
      pop(xml_out);

      pop(xml_out);  // prop_colorvec

      // Flush the object-std::map 
      prop_obj.flush();

      // Write the meta-data for this operator
      {
	XMLBufferWriter file_xml;

	push(file_xml, "PropColorVectors");
	write(file_xml, "num_records", prop_obj.size()); 
	write(file_xml, "Params", params.param);
	write(file_xml, "Config_info", gauge_xml);
	pop(file_xml);

	XMLBufferWriter record_xml;
	push(record_xml, "PropColorVector");
	write(record_xml, "num_records", prop_obj.size()); 
	pop(record_xml);

	// Write the propagator xml info
	TheNamedObjMap::Instance().get(params.named_obj.prop_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.prop_id).setRecordXML(record_xml);
      }

      snoop.stop();
      QDPIO::cout << name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << std::endl;

      QDPIO::cout << name << ": ran successfully" << std::endl;

      END_CODE();
    } 

  }

} // namespace Chroma
