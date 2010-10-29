/*! \file
 * \brief Compute the matrix element of  LatticeColorVector*M^-1*LatticeColorVector
 *
 * Propagator calculation on a colorvector
 */

#include "fermact.h"
#include "meas/inline/hadron/inline_prop_matelem_colorvec_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ferm/subset_vectors.h"
#include "qdp_map_obj.h"
#include "util/ferm/key_prop_colorvec.h"
#include "util/ferm/key_prop_matelem.h"
#include "util/ferm/key_val_db.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlinePropMatElemColorVecEnv 
  {
    //! Propagator input
    void read(XMLReader& xml, const string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "colorvec_id", input.colorvec_id);
      read(inputtop, "prop_id", input.prop_id);
      read(inputtop, "prop_op_file", input.prop_op_file);
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "colorvec_id", input.colorvec_id);
      write(xml, "prop_id", input.prop_id);
      write(xml, "prop_op_file", input.prop_op_file);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const string& path, Params::Param_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "num_vecs", input.num_vecs);
      read(inputtop, "t_sources", input.t_sources);
      read(inputtop, "decay_dir", input.decay_dir);
      read(inputtop, "mass_label", input.mass_label);
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const Params::Param_t& input)
    {
      push(xml, path);

      write(xml, "num_vecs", input.num_vecs);
      write(xml, "t_sources", input.t_sources);
      write(xml, "decay_dir", input.decay_dir);
      write(xml, "mass_label", input.mass_label);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const string& path, Params& input)
    {
      Params tmp(xml, path);
      input = tmp;
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const Params& input)
    {
      push(xml, path);
    
      write(xml, "Param", input.param);
      write(xml, "NamedObject", input.named_obj);

      pop(xml);
    }


    //---------------------------------------------------------------------------------
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
      
    const std::string name = "PROP_MATELEM_COLORVEC";

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

	push(xml_out, "PropMatElemColorVec");
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

      push(xml_out, "PropMatElemColorVec");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": propagator colorvector matrix element calculation" << endl;

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

      QDPIO::cout << "Snarf the source from a named buffer. Check for the prop map" << endl;
      try
      {	
	// NB We are just checking this is here.
	TheNamedObjMap::Instance().getData< Handle< QDP::MapObject<int,EVPair<LatticeColorVector> > > >(params.named_obj.colorvec_id);
	TheNamedObjMap::Instance().getData< Handle< QDP::MapObject<KeyPropColorVec_t,LatticeFermion> > >(params.named_obj.prop_id);

	// Snarf the source info. This is will throw if the colorvec_id is not there
	TheNamedObjMap::Instance().get(params.named_obj.colorvec_id).getFileXML(source_file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.colorvec_id).getRecordXML(source_record_xml);

	// Write out the source header
	write(xml_out, "Source_file_info", source_file_xml);
	write(xml_out, "Source_record_info", source_record_xml);
      }    
      catch (std::bad_cast)
      {
	QDPIO::cerr << name << ": caught dynamic cast error" << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << name << ": error extracting source_header or prop map: " << e << endl;
	QDP_abort(1);
      }

      // Cast should be valid now
      const QDP::MapObject<int,EVPair<LatticeColorVector> >& eigen_source = 
	*(TheNamedObjMap::Instance().getData< Handle< QDP::MapObject<int,EVPair<LatticeColorVector> > > >(params.named_obj.colorvec_id));

      // Cast should be valid now
      const QDP::MapObject<KeyPropColorVec_t,LatticeFermion>& prop_obj =
	*(TheNamedObjMap::Instance().getData< Handle< QDP::MapObject<KeyPropColorVec_t,LatticeFermion> > >(params.named_obj.prop_id));

      QDPIO::cout << "Source and prop map successfully found and parsed" << endl;

      // Sanity check - write out the norm2 of the source in the Nd-1 direction
      // Use this for any possible verification
      {
	// Initialize the slow Fourier transform phases
	SftMom phases(0, true, Nd-1);

	EVPair<LatticeColorVector> tmpvec; eigen_source.get(0,tmpvec);
	multi1d<Double> source_corrs = sumMulti(localNorm2(tmpvec.eigenVector), phases.getSet());

	push(xml_out, "Source_correlators");
	write(xml_out, "source_corrs", source_corrs);
	pop(xml_out);
      }

      // Another sanity check
      if (params.param.num_vecs > eigen_source.size())
      {
	QDPIO::cerr << __func__ << ": num_vecs= " << params.param.num_vecs
		    << " is greater than the number of available colorvectors= "
		    << eigen_source.size() << endl;
	QDP_abort(1);
      }


      //
      // DB storage
      //
      BinaryStoreDB< SerialDBKey<KeyPropElementalOperator_t>, SerialDBData<ValPropElementalOperator_t> > qdp_db;

      // Open the file, and write the meta-data and the binary for this operator
      if (! qdp_db.fileExists(params.named_obj.prop_op_file))
      {
	XMLBufferWriter file_xml;

	push(file_xml, "DBMetaData");
	write(file_xml, "id", string("propElemOp"));
	write(file_xml, "lattSize", QDP::Layout::lattSize());
	write(file_xml, "decay_dir", params.param.decay_dir);
	proginfo(file_xml);    // Print out basic program info
	write(file_xml, "Params", params.param);
	write(file_xml, "Config_info", gauge_xml);
	write(file_xml, "Weights", getEigenValues(eigen_source, params.param.num_vecs));
	pop(file_xml);

	std::string file_str(file_xml.str());
	qdp_db.setMaxUserInfoLen(file_str.size());

	qdp_db.open(params.named_obj.prop_op_file, O_RDWR | O_CREAT, 0664);

	qdp_db.insertUserdata(file_str);
      }
      else
      {
	qdp_db.open(params.named_obj.prop_op_file, O_RDWR, 0664);
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
	swatch.start();

	//
	// Loop over the source color and spin, creating the source
	// and calling the relevant propagator routines.
	//
	const int num_vecs            = params.param.num_vecs;
	const int decay_dir           = params.param.decay_dir;
	const multi1d<int>& t_sources = params.param.t_sources;

	// Initialize the slow Fourier transform phases
	SftMom phases(0, true, decay_dir);

	// Binary output
	push(xml_out, "ColorVecMatElems");

	// Loop over each time source
	for(int tt=0; tt < t_sources.size(); ++tt)
	{
	  int t_source = t_sources[tt];
	  QDPIO::cout << "t_source = " << t_source << endl; 

	  //
	  // All the loops
	  //
	  // NOTE: I pull the spin source and sink loops to the outside intentionally.
	  // The idea is to create a colorvector index (2d) array. These are not
	  // too big, but are big enough to make the IO efficient, and the DB efficient
	  // on reading. For N=32 and Lt=128, the mats are 2MB.
	  //
	  // Another NOTE: the colorvec_source is pulled all the way out as far as possible.
          // This means there is the minimal number of lookups of the fermion solution vector.
          // However, it also means that all spin_sink worth of a perambulator must be
          // held in memory.
	  //
	  for(int spin_source=0; spin_source < Ns; ++spin_source)
	  {
	    QDPIO::cout << "spin_source = " << spin_source << endl; 

	    //
	    // Initialize all the keys for this spin spin and all spin sinks on all t-slices
	    //
	    multi2d<KeyValPropElementalOperator_t> buf(phases.numSubsets(),Ns);
	    for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
	    {
	      for(int t=0; t < phases.numSubsets(); ++t)
	      {
		buf(t,spin_sink).key.key().t_slice      = t;
		buf(t,spin_sink).key.key().t_source     = t_source;
		buf(t,spin_sink).key.key().spin_src     = spin_source;
		buf(t,spin_sink).key.key().spin_snk     = spin_sink;
		buf(t,spin_sink).key.key().mass_label   = params.param.mass_label;
		buf(t,spin_sink).val.data().mat.resize(num_vecs,num_vecs);
	      }
	    }

	    for(int colorvec_source=0; colorvec_source < num_vecs; ++colorvec_source)
	    {
	      KeyPropColorVec_t key;
	      key.t_source     = t_source;
	      key.colorvec_src = colorvec_source;
	      key.spin_src     = spin_source;
		  
              multi1d<LatticeColorVector> vec_source(Ns);
              {
	        LatticeFermion tmp; prop_obj.get(key, tmp);
		for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
		{
                  vec_source[spin_sink] = peekSpin(tmp, spin_sink);
                }
              }

	      for(int colorvec_sink=0; colorvec_sink < num_vecs; ++colorvec_sink)
	      {
		EVPair<LatticeColorVector> vec_sink; eigen_source.get(colorvec_sink,vec_sink);

		for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
		{
		  multi1d<ComplexD> hsum(sumMulti(localInnerProduct(vec_sink.eigenVector, vec_source[spin_sink]), phases.getSet()));

		  for(int t=0; t < hsum.size(); ++t)
		  {
		    buf(t,spin_sink).val.data().mat(colorvec_sink,colorvec_source) = hsum[t];
		  }

		} // for spin_sink
	      } // for colorvec_sink
	    } // for colorvec_source
	      
	    for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
	    {
	      QDPIO::cout << "insert: spin_source= " << spin_source << " spin_sink= " << spin_sink << endl; 
	      for(int t=0; t < phases.numSubsets(); ++t)
	      {
		qdp_db.insert(buf(t,spin_sink).key, buf(t,spin_sink).val);
	      }
	    } // for spin_sink

	  } // for spin_source
	} // for t_source

	pop(xml_out);

	swatch.stop();
	QDPIO::cout << "Propagators computed: time= " 
		    << swatch.getTimeInSeconds() 
		    << " secs" << endl;
      }
      catch (const std::string& e) 
      {
	QDPIO::cout << name << ": caught exception around qprop: " << e << endl;
	QDP_abort(1);
      }

      push(xml_out,"Relaxation_Iterations");
      write(xml_out, "ncg_had", ncg_had);
      pop(xml_out);

      pop(xml_out);  // prop_matelem_colorvec

      snoop.stop();
      QDPIO::cout << name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << name << ": ran successfully" << endl;

      END_CODE();
    } 

  }

} // namespace Chroma
