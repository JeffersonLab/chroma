// $Id: inline_genprop_matelem_colorvec_w.cc,v 1.6 2008-10-15 20:09:24 edwards Exp $
/*! \file
 * \brief Compute the matrix element of  LatticeColorVector*M^-1*Gamma*M^-1**LatticeColorVector
 *
 * Generalized propagator calculation on a colorvector
 */

#include "handle.h"
#include "meas/inline/hadron/inline_genprop_matelem_colorvec_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/displace.h"
#include "meas/glue/mesplq.h"
#include "util/ferm/map_obj.h"
#include "util/ferm/key_val_db.h"
#include "util/ferm/key_prop_colorvec.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  /*!
   * \ingroup hadron
   *
   * @{
   */
  namespace InlineGenPropMatElemColorVecEnv 
  { 
    // Reader for input parameters
    void read(XMLReader& xml, const string& path, Params::Param_t::DispGamma_t& param)
    {
      XMLReader paramtop(xml, path);
    
      read(paramtop, "gamma", param.gamma);
      read(paramtop, "displacement", param.displacement);
    }

    // Reader for input parameters
    void read(XMLReader& xml, const string& path, Params::Param_t& param)
    {
      XMLReader paramtop(xml, path);
    
      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
	/**************************************************************************/
	break;

      default :
	/**************************************************************************/

	QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
	QDP_abort(1);
      }

      read(paramtop, "t_source", param.t_source);
      read(paramtop, "t_sink", param.t_sink);
      read(paramtop, "mom2_max", param.mom2_max);
      read(paramtop, "mom_offset", param.mom_offset);
      read(paramtop, "displacement_length", param.displacement_length);
      read(paramtop, "DisplacementGammaList", param.disp_gamma_list);
      read(paramtop, "num_vecs", param.num_vecs);
      read(paramtop, "decay_dir", param.decay_dir);
      read(paramtop, "mass_label", param.mass_label);

      param.link_smearing  = readXMLGroup(paramtop, "LinkSmearing", "LinkSmearingType");
    }


    // Writer for input parameters
    void write(XMLWriter& xml, const string& path, const Params::Param_t::DispGamma_t& param)
    {
      push(xml, path);

      write(xml, "gamma", param.gamma);
      write(xml, "displacement", param.displacement);

      pop(xml);
    }

    // Writer for input parameters
    void write(XMLWriter& xml, const string& path, const Params::Param_t& param)
    {
      push(xml, path);

      int version = 1;

      write(xml, "version", version);
      write(xml, "t_source", param.t_source);
      write(xml, "t_sink", param.t_sink);
      write(xml, "mom_offset", param.mom_offset);
      write(xml, "mom2_max", param.mom2_max);
      write(xml, "displacement_length", param.displacement_length);
      write(xml, "DisplacementGammaList", param.disp_gamma_list);
      write(xml, "num_vecs", param.num_vecs);
      write(xml, "decay_dir", param.decay_dir);
      write(xml, "mass_label", param.mass_label);

      xml << param.link_smearing.xml;

      pop(xml);
    }

    //! Read named objects 
    void read(XMLReader& xml, const string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "source_prop_id", input.source_prop_id);
      read(inputtop, "sink_prop_id", input.sink_prop_id);
      read(inputtop, "genprop_op_file", input.genprop_op_file);
    }

    //! Write named objects
    void write(XMLWriter& xml, const string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "source_prop_id", input.source_prop_id);
      write(xml, "sink_prop_id", input.sink_prop_id);
      write(xml, "genprop_op_file", input.genprop_op_file);

      pop(xml);
    }

    // Writer for input parameters
    void write(XMLWriter& xml, const string& path, const Params& param)
    {
      param.writeXML(xml, path);
    }
  }


  namespace InlineGenPropMatElemColorVecEnv 
  { 
    // Anonymous namespace for registration
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

    const std::string name = "GENPROP_MATELEM_COLORVEC";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= LinkSmearingEnv::registerAll();
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }


    //! Anonymous namespace
    /*! Diagnostic stuff */
    namespace
    {
      StandardOutputStream& operator<<(StandardOutputStream& os, const multi1d<int>& d)
      {
	if (d.size() > 0)
	{
	  os << d[0];

	  for(int i=1; i < d.size(); ++i)
	    os << " " << d[i];
	}

	return os;
      }
    }


    //----------------------------------------------------------------------------
    // Param stuff
    Params::Params()
    { 
      frequency = 0; 
      param.mom2_max = 0;
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


    void
    Params::writeXML(XMLWriter& xml_out, const std::string& path) const
    {
      push(xml_out, path);
    
      // Parameters for source construction
      write(xml_out, "Param", param);

      // Write out the output propagator/source configuration info
      write(xml_out, "NamedObject", named_obj);

      pop(xml_out);
    }


    //----------------------------------------------------------------------------
    //! Generalized propagator operator
    struct KeyGenPropElementalOperator_t
    {
      int                t_slice;       /*!< Propagator time slice */
      int                t_source;      /*!< Source time slice */
      int                t_sink;        /*!< Source time slice */
      int                spin_l;        /*!< Source spin index */
      int                spin_r;        /*!< Sink spin index */
      int                gamma;         /*!< The gamma matrix number - [0,Ns^2) */
      multi1d<int>       displacement;  /*!< Displacement dirs of right colorvector */
      multi1d<int>       mom;           /*!< D-1 momentum of this operator */
      std::string        mass_label;    /*!< A mass label */
    };

    //! Generalized propagator operator
    struct ValGenPropElementalOperator_t
    {
      multi2d<ComplexD>  op;              /*!< Colorvector source and sink with momentum projection */
    };


    //----------------------------------------------------------------------------
    //! Holds key and value as temporaries
    struct KeyValGenPropElementalOperator_t
    {
      SerialDBKey<KeyGenPropElementalOperator_t>  key;
      SerialDBData<ValGenPropElementalOperator_t> val;
    };


    //----------------------------------------------------------------------------
    //! KeyGenPropElementalOperator reader
    void read(BinaryReader& bin, KeyGenPropElementalOperator_t& param)
    {
      read(bin, param.t_slice);
      read(bin, param.t_source);
      read(bin, param.t_sink);
      read(bin, param.spin_l);
      read(bin, param.spin_r);
      read(bin, param.gamma);
      read(bin, param.displacement);
      read(bin, param.mom);
      read(bin, param.mass_label, 32);
    }

    //! GenPropElementalOperator write
    void write(BinaryWriter& bin, const KeyGenPropElementalOperator_t& param)
    {
      write(bin, param.t_slice);
      write(bin, param.t_source);
      write(bin, param.t_sink);
      write(bin, param.spin_l);
      write(bin, param.spin_r);
      write(bin, param.gamma);
      write(bin, param.displacement);
      write(bin, param.mom);
      write(bin, param.mass_label);
    }


    //----------------------------------------------------------------------------
    //! PropElementalOperator reader
    void read(BinaryReader& bin, ValGenPropElementalOperator_t& param)
    {
      int n;
      read(bin, n);    // the size is always written, even if 0
      param.op.resize(n,n);
  
      for(int i=0; i < n; ++i)
      {
	for(int j=0; j < n; ++j)
	{
	  read(bin, param.op(i,j));
	}
      }
    }

    //! GenPropElementalOperator write
    void write(BinaryWriter& bin, const ValGenPropElementalOperator_t& param)
    {
      int n = param.op.size1();  // all sizes the same
      write(bin, n);
      for(int i=0; i < n; ++i)
      {
	for(int j=0; j < n; ++j)
	{
	  write(bin, param.op(i,j));
	}
      }
    }


    //----------------------------------------------------------------------------
    //! Normalize just one displacement array
    multi1d<int> normDisp(const multi1d<int>& orig)
    {
      START_CODE();

      multi1d<int> disp;
      multi1d<int> empty; 
      multi1d<int> no_disp(1); no_disp[0] = 0;

      // NOTE: a no-displacement is recorded as a zero-length array
      // Convert a length one array with no displacement into a no-displacement array
      if (orig.size() == 1)
      {
	if (orig == no_disp)
	  disp = empty;
	else
	  disp = orig;
      }
      else
      {
	disp = orig;
      }

      END_CODE();

      return disp;
    } // void normDisp


    //-------------------------------------------------------------------------------
    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "GenPropMatElemColorVec");
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


    // Function call
    void 
    InlineMeas::func(unsigned long update_no,
		     XMLWriter& xml_out) 
    {
      START_CODE();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      StopWatch swiss;
			
      push(xml_out, "GenPropMatElemColorVec");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": Generalized propagator color-vector matrix element" << endl;

      // Test and grab a reference to the gauge field
      XMLBufferWriter gauge_xml;
      XMLReader source_prop_file_xml, source_prop_record_xml;
      XMLReader sink_prop_file_xml, sink_prop_record_xml;
      try
      {
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
	TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);

	TheNamedObjMap::Instance().getData< MapObject<KeyPropColorVec_t,LatticeFermion> >(params.named_obj.source_prop_id);
	TheNamedObjMap::Instance().getData< MapObject<KeyPropColorVec_t,LatticeFermion> >(params.named_obj.sink_prop_id);

	// Snarf the prop info. This is will throw if the prop_id is not there
	TheNamedObjMap::Instance().get(params.named_obj.source_prop_id).getFileXML(source_prop_file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.source_prop_id).getRecordXML(source_prop_record_xml);

	// Snarf the prop info. This is will throw if the prop_id is not there
	TheNamedObjMap::Instance().get(params.named_obj.sink_prop_id).getFileXML(sink_prop_file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.sink_prop_id).getRecordXML(sink_prop_record_xml);
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
      // Cast should be valid now
      const multi1d<LatticeColorMatrix>& u = 
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

      const MapObject<KeyPropColorVec_t,LatticeFermion>& source_ferm_map =
	TheNamedObjMap::Instance().getData< MapObject<KeyPropColorVec_t,LatticeFermion> >(params.named_obj.source_prop_id);
      const MapObject<KeyPropColorVec_t,LatticeFermion>& sink_ferm_map =
	TheNamedObjMap::Instance().getData< MapObject<KeyPropColorVec_t,LatticeFermion> >(params.named_obj.sink_prop_id);

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      proginfo(xml_out);    // Print out basic program info

      // Write out the input
      params.writeXML(xml_out, "Input");

      // Write out the source header
      write(xml_out, "Source_prop_file_info", source_prop_file_xml);
      write(xml_out, "Source_prop_record_info", source_prop_record_xml);
      write(xml_out, "Sink_prop_file_info", sink_prop_file_xml);
      write(xml_out, "Sink_prop_record_info", sink_prop_record_xml);

      // Write out the config info
      write(xml_out, "Config_info", gauge_xml);

      //First calculate some gauge invariant observables just for info.
      //This is really cheap.
      MesPlq(xml_out, "Observables", u);

      //
      // Initialize the slow Fourier transform phases
      //
      multi1d<int> origin_offset(Nd);
      origin_offset = 0;
      SftMom phases(params.param.mom2_max, origin_offset, params.param.mom_offset, false, params.param.decay_dir);

      //
      // Smear the gauge field if needed
      //
      multi1d<LatticeColorMatrix> u_smr = u;

      try
      {
	std::istringstream  xml_l(params.param.link_smearing.xml);
	XMLReader  linktop(xml_l);
	QDPIO::cout << "Link smearing type = " << params.param.link_smearing.id << endl;
	
	
	Handle< LinkSmearing >
	  linkSmearing(TheLinkSmearingFactory::Instance().createObject(params.param.link_smearing.id,
								       linktop, 
								       params.param.link_smearing.path));

	(*linkSmearing)(u_smr);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": Caught Exception link smearing: " << e << endl;
	QDP_abort(1);
      }

      MesPlq(xml_out, "Smeared_Observables", u_smr);

      //
      // Generalized propagatos
      //
      QDPIO::cout << "Building generalized propagators" << endl;

      // DB storage
      BinaryFxStoreDB< SerialDBKey<KeyGenPropElementalOperator_t>, SerialDBData<ValGenPropElementalOperator_t> > 
        qdp_db(params.named_obj.genprop_op_file, 10*1024*1024, 64*1024);

      push(xml_out, "ElementalOps");

      // Loop over all time slices for the source. This is the same 
      // as the subsets for  phases

      const int num_vecs  = params.param.num_vecs;
      const int decay_dir = params.param.decay_dir;
      const int t_source  = params.param.t_source;
      const int t_sink    = params.param.t_sink;

      // Loop over each operator 
      for(int l=0; l < params.param.disp_gamma_list.size(); ++l)
      {
	StopWatch watch;

	QDPIO::cout << "Elemental operator: op = " << l << endl;

	// Make sure displacement is something sensible
	multi1d<int> disp = normDisp(params.param.disp_gamma_list[l].displacement);

	// Fold in the gamma_5 associated with hermiticity of the sink. 
	// Can multiply the desired Gamma on the left by gamma_5
	int gamma = params.param.disp_gamma_list[l].gamma;
	int gamma_tmp = (Ns*Ns-1) ^ gamma;

	QDPIO::cout << "gamma=" << gamma_tmp << "  displacement= " << disp << endl;

	// Build the operator
	swiss.reset();
	swiss.start();

	// Big loop over the momentum projection
	for(int mom_num = 0 ; mom_num < phases.numMom() ; ++mom_num) 
	{
	  // Loop over spins
	  for(int spin_r=0; spin_r < Ns; ++spin_r)
	  {
	    QDPIO::cout << "spin_r = " << spin_r << endl; 

	    for(int spin_l=0; spin_l < Ns; ++spin_l)
	    {
	      QDPIO::cout << "spin_l = " << spin_l << endl; 

	      // The keys for the spin and displacements for this particular elemental operator
	      // No displacement for left colorvector, only displace right colorvector
	      // Invert the time - make it an independent key
	      multi1d<KeyValGenPropElementalOperator_t> buf(phases.numSubsets());
	      for(int t=0; t < phases.numSubsets(); ++t)
	      {
		buf[t].key.key().t_slice       = t;
		buf[t].key.key().t_source      = t_source;
		buf[t].key.key().t_sink        = t_sink;
		buf[t].key.key().spin_r        = spin_r;
		buf[t].key.key().spin_l        = spin_l;
		buf[t].key.key().mass_label    = params.param.mass_label;
		buf[t].key.key().mom           = phases.numToMom(mom_num);
		buf[t].key.key().gamma         = gamma;
		buf[t].key.key().displacement  = disp; // only right colorvector

		buf[t].val.data().op.resize(num_vecs, num_vecs);
	      }

	      for(int j = 0; j < params.param.num_vecs; ++j)
	      {
		KeyPropColorVec_t key_r;
		key_r.t_source     = t_source;
		key_r.colorvec_src = j;
		key_r.spin_src     = spin_r;
		  
		// Displace the right vector and multiply by the momentum phase
		LatticeFermion shift_ferm = Gamma(gamma_tmp) * displace(u_smr, 
									source_ferm_map[key_r],
									params.param.displacement_length, 
									disp);

		for(int i = 0; i < params.param.num_vecs; ++i)
		{
		  KeyPropColorVec_t key_l;
		  key_l.t_source     = t_sink;
		  key_l.colorvec_src = i;
		  key_l.spin_src     = spin_l;
		  
		  watch.reset();
		  watch.start();

		  // Contract over color indices
		  // Do the relevant quark contraction
		  LatticeComplex lop = localInnerProduct(sink_ferm_map[key_l], shift_ferm);

		  // Slow fourier-transform
		  multi1d<ComplexD> op_sum = sumMulti(phases[mom_num] * lop, phases.getSet());

		  watch.stop();

		  for(int t=0; t < op_sum.size(); ++t)
		  {
		    buf[t].val.data().op(i,j) = op_sum[t];
		  }
		} // end for i
	      } // end for j

	      QDPIO::cout << "insert: mom= " << phases.numToMom(mom_num) << " displacement= " << disp << endl; 
	      for(int t=0; t < phases.numSubsets(); ++t)
	      {
		qdp_db.insert(buf[t].key, buf[t].val);
	      }

	    } // end for spin_l
	  } // end for spin_r
	} // mom_num

	swiss.stop();

	QDPIO::cout << "GenProp operator= " << l 
		    << "  time= "
		    << swiss.getTimeInSeconds() 
		    << " secs" << endl;

      } // for l

      pop(xml_out); // ElementalOps

      // Write the meta-data and the binary for this operator
      swiss.reset();
      swiss.start();
      {
	XMLBufferWriter file_xml;

	push(file_xml, "PropElementalOperators");
	write(file_xml, "lattSize", QDP::Layout::lattSize());
	write(file_xml, "decay_dir", params.param.decay_dir);
	write(file_xml, "Params", params.param);
	write(file_xml, "Op_Info", params.param.disp_gamma_list);
	write(file_xml, "Source_prop_record_info", source_prop_record_xml);
	write(file_xml, "Sink_prop_record_info", sink_prop_record_xml);
	write(file_xml, "Config_info", gauge_xml);
	pop(file_xml);

	qdp_db.insertUserdata(file_xml.str());
      }
      swiss.stop();

      QDPIO::cout << "Meson Operator written:"
		  << "  time= " << swiss.getTimeInSeconds() << " secs" << endl;

      // Close the namelist output file XMLDAT
      pop(xml_out);     // GenPropMatElemColorVector

      snoop.stop();
      QDPIO::cout << name << ": total time = " 
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << name << ": ran successfully" << endl;

      END_CODE();
    } // func
  } // namespace InlineGenPropMatElemColorVecEnv

  /*! @} */  // end of group hadron

} // namespace Chroma
