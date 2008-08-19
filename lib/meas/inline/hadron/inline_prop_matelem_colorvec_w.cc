// $Id: inline_prop_matelem_colorvec_w.cc,v 1.9 2008-08-19 04:59:40 edwards Exp $
/*! \file
 * \brief Compute the matrix element of  LatticeColorVector*M^-1*LatticeColorVector
 *
 * Propagator calculation on a colorvector
 */

#include "fermact.h"
#include "meas/inline/hadron/inline_prop_matelem_colorvec_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ferm/eigeninfo.h"
#include "util/ferm/map_obj.h"
#include "util/ferm/key_prop_colorvec.h"
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
    void read(XMLReader& xml, const string& path, InlinePropMatElemColorVecEnv::Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "colorvec_id", input.colorvec_id);
      read(inputtop, "prop_id", input.prop_id);
      read(inputtop, "prop_op_file", input.prop_op_file);
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const InlinePropMatElemColorVecEnv::Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "colorvec_id", input.colorvec_id);
      write(xml, "prop_id", input.prop_id);
      write(xml, "prop_op_file", input.prop_op_file);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const string& path, InlinePropMatElemColorVecEnv::Params::Param_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "num_vecs", input.num_vecs);
      read(inputtop, "t_sources", input.t_sources);
      read(inputtop, "decay_dir", input.decay_dir);
      read(inputtop, "mass_label", input.mass_label);
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const InlinePropMatElemColorVecEnv::Params::Param_t& input)
    {
      push(xml, path);

      write(xml, "num_vecs", input.num_vecs);
      write(xml, "t_sources", input.t_sources);
      write(xml, "decay_dir", input.decay_dir);
      write(xml, "mass_label", input.mass_label);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const string& path, InlinePropMatElemColorVecEnv::Params& input)
    {
      InlinePropMatElemColorVecEnv::Params tmp(xml, path);
      input = tmp;
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const InlinePropMatElemColorVecEnv::Params& input)
    {
      push(xml, path);
    
      write(xml, "Param", input.param);
      write(xml, "NamedObject", input.named_obj);

      pop(xml);
    }
  } // namespace InlinePropMatElemColorVecEnv 


  namespace InlinePropMatElemColorVecEnv 
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
    //! Prop operator
    struct KeyPropElementalOperator_t
    {
      int                t_source;      /*!< Source time slice */
      int                colorvec_src;  /*!< Source colorvector index */
      int                colorvec_snk;  /*!< Sink colorvector index */
      int                spin_src;      /*!< Source spin index */
      int                spin_snk;      /*!< Sink spin index */
      std::string        mass_label;    /*!< A mass label */
    };


    //! Prop operator
    struct ValPropElementalOperator_t
    {
      multi1d<ComplexD>  corr;          /*!< Projected propagator */
    };


    //----------------------------------------------------------------------------
    //! PropElementalOperator reader
    void read(BinaryReader& bin, KeyPropElementalOperator_t& param)
    {
      read(bin, param.t_source);
      read(bin, param.colorvec_src);
      read(bin, param.colorvec_snk);
      read(bin, param.spin_src);
      read(bin, param.spin_snk);
      read(bin, param.mass_label, 32);
    }

    //! PropElementalOperator write
    void write(BinaryWriter& bin, const KeyPropElementalOperator_t& param)
    {
      write(bin, param.t_source);
      write(bin, param.colorvec_src);
      write(bin, param.colorvec_snk);
      write(bin, param.spin_src);
      write(bin, param.spin_snk);
      write(bin, param.mass_label);
    }

    //! PropElementalOperator reader
    void read(XMLReader& xml, const std::string& path, KeyPropElementalOperator_t& param)
    {
      XMLReader paramtop(xml, path);
    
      read(paramtop, "t_source", param.t_source);
      read(paramtop, "colorvec_src", param.colorvec_src);
      read(paramtop, "colorvec_snk", param.colorvec_snk);
      read(paramtop, "spin_src", param.spin_src);
      read(paramtop, "spin_snk", param.spin_snk);
      read(paramtop, "mass_label", param.mass_label);
    }

    //! PropElementalOperator writer
    void write(XMLWriter& xml, const std::string& path, const KeyPropElementalOperator_t& param)
    {
      push(xml, path);

      write(xml, "t_source", param.t_source);
      write(xml, "colorvec_src", param.colorvec_src);
      write(xml, "colorvec_snk", param.colorvec_snk);
      write(xml, "spin_src", param.spin_src);
      write(xml, "spin_snk", param.spin_snk);
      write(xml, "mass_label", param.mass_label);

      pop(xml);
    }


    //----------------------------------------------------------------------------
    //! PropElementalOperator reader
    void read(BinaryReader& bin, ValPropElementalOperator_t& param)
    {
      read(bin, param.corr);
    }

    //! PropElementalOperator write
    void write(BinaryWriter& bin, const ValPropElementalOperator_t& param)
    {
      write(bin, param.corr);
    }

    //! PropElementalOperator reader
    void read(XMLReader& xml, const std::string& path, ValPropElementalOperator_t& param)
    {
      XMLReader paramtop(xml, path);
    
      read(paramtop, "corr", param.corr);
    }

    //! PropElementalOperator writer
    void write(XMLWriter& xml, const std::string& path, const ValPropElementalOperator_t& param)
    {
      push(xml, path);

      write(xml, "corr", param.corr);

      pop(xml);
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
      // Read in the source along with relevant information.
      // 
      XMLReader source_file_xml, source_record_xml;

      QDPIO::cout << "Snarf the source from a named buffer. Check for the prop map" << endl;
      try
      {
	TheNamedObjMap::Instance().getData< EigenInfo<LatticeColorVector> >(params.named_obj.colorvec_id);
	TheNamedObjMap::Instance().getData< MapObject<KeyPropColorVec_t,LatticeFermion> >(params.named_obj.prop_id);

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
      const EigenInfo<LatticeColorVector>& eigen_source = 
	TheNamedObjMap::Instance().getData< EigenInfo<LatticeColorVector> >(params.named_obj.colorvec_id);

      // Cast should be valid now
      const MapObject<KeyPropColorVec_t,LatticeFermion>& map_obj =
	TheNamedObjMap::Instance().getData< MapObject<KeyPropColorVec_t,LatticeFermion> >(params.named_obj.prop_id);

      QDPIO::cout << "Source and prop map successfully found and parsed" << endl;

      // Sanity check - write out the norm2 of the source in the Nd-1 direction
      // Use this for any possible verification
      {
	// Initialize the slow Fourier transform phases
	SftMom phases(0, true, Nd-1);

	multi1d< multi1d<Double> > source_corrs(eigen_source.getEvalues().size());
	for(int m=0; m < source_corrs.size(); ++m)
	{
	  source_corrs[m] = sumMulti(localNorm2(eigen_source.getEvectors()[m]), phases.getSet());
	}

	push(xml_out, "Source_correlators");
	write(xml_out, "source_corrs", source_corrs);
	pop(xml_out);
      }

      // Another sanity check
      if (params.param.num_vecs > eigen_source.getEvalues().size())
      {
	QDPIO::cerr << __func__ << ": num_vecs= " << params.param.num_vecs
		    << " is greater than the number of available colorvectors= "
		    << eigen_source.getEvalues().size() << endl;
	QDP_abort(1);
      }


      // Total number of iterations
      int ncg_had = 0;

      // DB storage
      BinaryFxStoreDB< SerialDBKey<KeyPropElementalOperator_t>, SerialDBData<ValPropElementalOperator_t> > 
	qdp_db(params.named_obj.prop_op_file);

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

	// Loop over each operator 
	for(int tt=0; tt < t_sources.size(); ++tt)
	{
	  int t_source = t_sources[tt];
	  QDPIO::cout << "t_source = " << t_source << endl; 

	  // All the loops
	  for(int colorvec_source=0; colorvec_source < num_vecs; ++colorvec_source)
	  {
	    QDPIO::cout << "colorvec_source = " << colorvec_source << endl; 

	    // Pull out a time-slice of the color vector source
	    LatticeColorVector vec_srce = zero;
	    vec_srce[phases.getSet()[t_source]] = eigen_source.getEvectors()[colorvec_source];
	
	    for(int spin_source=0; spin_source < Ns; ++spin_source)
	    {
	      QDPIO::cout << "spin_source = " << spin_source << endl; 

	      KeyPropColorVec_t key;
	      key.t_source     = t_source;
	      key.colorvec_src = colorvec_source;
	      key.spin_src     = spin_source;
		  
	      LatticeFermion quark_soln(map_obj[key]);

	      for(int colorvec_sink=0; colorvec_sink < num_vecs; ++colorvec_sink)
	      {
		const LatticeColorVector& vec_sink = eigen_source.getEvectors()[colorvec_sink];

		for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
		{
		  SerialDBKey<KeyPropElementalOperator_t> key;
		  key.key().t_source     = t_source;
		  key.key().colorvec_src = colorvec_source;
		  key.key().colorvec_snk = colorvec_sink;
		  key.key().spin_src     = spin_source;
		  key.key().spin_snk     = spin_sink;
		  key.key().mass_label   = params.param.mass_label;

		  SerialDBData<ValPropElementalOperator_t> val;
		  val.data().corr        = sumMulti(localInnerProduct(vec_sink, peekSpin(quark_soln, spin_sink)), 
						    phases.getSet());
		  
		  qdp_db.insert(key, val);
		} // for spin_sink
	      } // for colorvec_sink
	    } // for spin_source
	  } // for colorvec_source
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

      // Write the meta-data and the binary for this operator
      {
	XMLBufferWriter file_xml;

	push(file_xml, "PropElementalOperators");
	write(file_xml, "Params", params.param);
	write(file_xml, "Config_info", gauge_xml);
	pop(file_xml);

	qdp_db.insertUserdata(file_xml.str());
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
