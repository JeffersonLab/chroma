// $Id: inline_genprop_matelem_colorvec_w.cc,v 1.1 2008-09-01 01:44:31 edwards Exp $
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
#include "meas/glue/mesplq.h"
#include "meas/smear/disp_colvec_map.h"
#include "util/ferm/eigeninfo.h"
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
    void read(XMLReader& xml, const string& path, InlineGenPropMatElemColorVecEnv::Params::Param_t& param)
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

      read(paramtop, "mom2_max", param.mom2_max);
      read(paramtop, "displacement_length", param.displacement_length);
      read(paramtop, "displacement_list", param.displacement_list);
      read(paramtop, "num_vecs", param.num_vecs);
      read(paramtop, "decay_dir", param.decay_dir);

      param.link_smearing  = readXMLGroup(paramtop, "LinkSmearing", "LinkSmearingType");
    }


    // Writer for input parameters
    void write(XMLWriter& xml, const string& path, const InlineGenPropMatElemColorVecEnv::Params::Param_t& param)
    {
      push(xml, path);

      int version = 1;

      write(xml, "version", version);
      write(xml, "mom2_max", param.mom2_max);
      write(xml, "displacement_length", param.displacement_length);
      write(xml, "displacement_list", param.displacement_list);
      write(xml, "num_vecs", param.num_vecs);
      write(xml, "decay_dir", param.decay_dir);
      xml << param.link_smearing.xml;

      pop(xml);
    }

    //! Read named objects 
    void read(XMLReader& xml, const string& path, InlineGenPropMatElemColorVecEnv::Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "colorvec_id", input.colorvec_id);
      read(inputtop, "prop_id", input.prop_id);
      read(inputtop, "genprop_op_file", input.genprop_op_file);
    }

    //! Write named objects
    void write(XMLWriter& xml, const string& path, const InlineGenPropMatElemColorVecEnv::Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "colorvec_id", input.colorvec_id);
      write(xml, "prop_id", input.prop_id);
      write(xml, "genprop_op_file", input.genprop_op_file);

      pop(xml);
    }

    // Writer for input parameters
    void write(XMLWriter& xml, const string& path, const InlineGenPropMatElemColorVecEnv::Params& param)
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

    const std::string name = "MESON_MATELEM_COLORVEC";

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
    //! Meson operator
    struct KeyGenPropElementalOperator_t
    {
      int                t_source;      /*!< Source time slice */
      int                t_sink;        /*!< Source time slice */
      int                colvec_l;      /*!< Left colorvector index */
      int                colvec_r;      /*!< Right colorvector index */
      int                spin_l;        /*!< Source spin index */
      int                spin_r;        /*!< Sink spin index */
      multi1d<int>       displacement;  /*!< Displacement dirs of right colorvector */
      multi1d<int>       mom;           /*!< D-1 momentum of this operator */
      std::string        mass_label;    /*!< A mass label */
    };

    //! Meson operator
    struct ValGenPropElementalOperator_t
    {
      multi1d<ComplexD>  op;           /*!< Momentum projected operator */
    };


    //----------------------------------------------------------------------------
    //! KeyGenPropElementalOperator reader
    void read(BinaryReader& bin, KeyGenPropElementalOperator_t& param)
    {
      read(bin, param.t_source);
      read(bin, param.t_sink);
      read(bin, param.colvec_l);
      read(bin, param.colvec_r);
      read(bin, param.spin_l);
      read(bin, param.spin_r);
      read(bin, param.displacement);
      read(bin, param.mom);
      read(bin, param.mass_label, 32);
    }

    //! GenPropElementalOperator write
    void write(BinaryWriter& bin, const KeyGenPropElementalOperator_t& param)
    {
      write(bin, param.t_source);
      write(bin, param.t_sink);
      write(bin, param.colvec_l);
      write(bin, param.colvec_r);
      write(bin, param.spin_l);
      write(bin, param.spin_r);
      write(bin, param.displacement);
      write(bin, param.mom);
      write(bin, param.mass_label);
    }

    //! GenPropElementalOperator reader
    void read(XMLReader& xml, const std::string& path, KeyGenPropElementalOperator_t& param)
    {
      XMLReader paramtop(xml, path);
    
      read(paramtop, "t_source", param.t_source);
      read(paramtop, "t_sink", param.t_sink);
      read(paramtop, "colvec_l", param.colvec_l);
      read(paramtop, "colvec_r", param.colvec_r);
      read(paramtop, "spin_l", param.spin_l);
      read(paramtop, "spin_r", param.spin_r);
      read(paramtop, "displacement", param.displacement);
      read(paramtop, "mom", param.mom);
      read(paramtop, "mass_label", param.mass_label);
    }

    //! GenPropElementalOperator writer
    void write(XMLWriter& xml, const std::string& path, const KeyGenPropElementalOperator_t& param)
    {
      push(xml, path);

      write(xml, "t_source", param.t_source);
      write(xml, "t_sink", param.t_sink);
      write(xml, "colvec_l", param.colvec_l);
      write(xml, "colvec_r", param.colvec_r);
      write(xml, "spin_l", param.spin_l);
      write(xml, "spin_r", param.spin_r);
      write(xml, "displacement", param.displacement);
      write(xml, "mom", param.mom);
      write(xml, "mass_label", param.mass_label);

      pop(xml);
    }


    //----------------------------------------------------------------------------
    //! GenPropElementalOperator reader
    void read(BinaryReader& bin, ValGenPropElementalOperator_t& param)
    {
      read(bin, param.op);
    }

    //! GenPropElementalOperator write
    void write(BinaryWriter& bin, const ValGenPropElementalOperator_t& param)
    {
      write(bin, param.op);
    }

    //! GenPropElementalOperator reader
    void read(XMLReader& xml, const std::string& path, ValGenPropElementalOperator_t& param)
    {
      XMLReader paramtop(xml, path);
    
      read(paramtop, "op", param.op);
    }

    //! GenPropElementalOperator writer
    void write(XMLWriter& xml, const std::string& path, const ValGenPropElementalOperator_t& param)
    {
      push(xml, path);

      write(xml, "op", param.op);

      pop(xml);
    }


    //----------------------------------------------------------------------------
    //! Make sure displacements are something sensible
    multi1d< multi1d<int> > normalizeDisplacements(const multi1d< multi1d<int> >& orig_list)
    {
      START_CODE();

      multi1d< multi1d<int> > displacement_list(orig_list.size());
      multi1d<int> empty(1); empty = 0;

      // Loop over displacements
      for(int n=0; n < orig_list.size(); ++n)
      {
	if (orig_list[n].size() == 0)
	{
	  displacement_list[n] = empty;
	}
	else
	{
	  displacement_list[n] = orig_list[n];
	}

	QDPIO::cout << "disp[" << n << "]= " << displacement_list[n] << endl;
      }

      END_CODE();

      return displacement_list;
    } // void normalizeDisplacements


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
			
      // Test and grab a reference to the gauge field
      XMLBufferWriter gauge_xml;
      XMLReader source_file_xml, source_record_xml;
      try
      {
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
	TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);

	TheNamedObjMap::Instance().getData< EigenInfo<LatticeColorVector> >(params.named_obj.colorvec_id).getEvectors();
	TheNamedObjMap::Instance().getData< MapObject<KeyPropColorVec_t,LatticeFermion> >(params.named_obj.prop_id);

	// Snarf the source info. This is will throw if the colorvec_id is not there
	TheNamedObjMap::Instance().get(params.named_obj.colorvec_id).getFileXML(source_file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.colorvec_id).getRecordXML(source_record_xml);

	// Write out the source header
	write(xml_out, "Source_file_info", source_file_xml);
	write(xml_out, "Source_record_info", source_record_xml);
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

      const EigenInfo<LatticeColorVector>& eigen_source = 
	TheNamedObjMap::Instance().getData< EigenInfo<LatticeColorVector> >(params.named_obj.colorvec_id);

      const MapObject<KeyPropColorVec_t,LatticeFermion>& map_obj =
	TheNamedObjMap::Instance().getData< MapObject<KeyPropColorVec_t,LatticeFermion> >(params.named_obj.prop_id);

      push(xml_out, "GenPropMatElemColorVec");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": Generalized propagator color-vector matrix element" << endl;

      proginfo(xml_out);    // Print out basic program info

      // Write out the input
      params.writeXML(xml_out, "Input");

      // Write out the config info
      write(xml_out, "Config_info", gauge_xml);

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      //First calculate some gauge invariant observables just for info.
      //This is really cheap.
      MesPlq(xml_out, "Observables", u);

      //
      // Initialize the slow Fourier transform phases
      //
      SftMom phases(params.param.mom2_max, false, params.param.decay_dir);

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
      // Make sure displacements are something sensible
      //
      QDPIO::cout << "Normalize displacement lengths" << endl;
      multi1d< multi1d<int> > displacement_list(normalizeDisplacements(params.param.displacement_list));

      for(int n=0; n < displacement_list.size(); ++n)
      {
	QDPIO::cout << "displa[" << n << "]= " << displacement_list[n] << endl;
      }

      // Keep track of no displacements and zero momentum
      multi1d<int> zero_displacement(1); zero_displacement = 0;
      multi1d<int> zero_mom(3); zero_mom = 0;

      //
      // The object holding the displaced color vector maps  
      //
      DispColorVectorMap smrd_disp_vecs(params.param.displacement_length,
					u_smr,
					eigen_source.getEvectors());

      //
      // Generalized propagatos
      //
      QDPIO::cout << "Building generalized propagators" << endl;

      // DB storage
      BinaryVarStoreDB< SerialDBKey<KeyGenPropElementalOperator_t>, SerialDBData<ValGenPropElementalOperator_t> > 
	qdp_db(params.named_obj.genprop_op_file);

      push(xml_out, "ElementalOps");

      // Loop over all time slices for the source. This is the same 
      // as the subsets for  phases

      const int num_vecs  = params.param.num_vecs;
      const int decay_dir = params.param.decay_dir;
      const int t_source  = params.param.t_source;
      const int t_sink    = params.param.t_sink;

      // Loop over each operator 
      for(int l=0; l < displacement_list.size(); ++l)
      {
	StopWatch watch;

	QDPIO::cout << "Elemental operator: op = " << l << endl;

	QDPIO::cout << "displacement = " << displacement_list[l] << endl;

	// Build the operator
	swiss.reset();
	swiss.start();

	// The keys for the spin and displacements for this particular elemental operator
	multi1d<KeyDispColorVector_t> keyDispColorVector(2);

	// No displacement for left colorvector, only displace right colorvector
	keyDispColorVector[0].displacement.resize(1);
	keyDispColorVector[0].displacement = 0;
	keyDispColorVector[1].displacement = displacement_list[l];

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

	    KeyPropColorVec_t key_source;
	    key.t_source     = t_source;
	    key.colorvec_src = colorvec_source;
	    key.spin_src     = spin_source;
	  
	    LatticeFermion quark_source(map_obj[key_source]);

	    for(int colorvec_sink=0; colorvec_sink < num_vecs; ++colorvec_sink)
	    {
	      KeyPropColorVec_t key_sink;
	      key.t_source     = t_sink;
	      key.colorvec_src = colorvec_sink;
	      key.spin_src     = spin_sink;
	  
	      LatticeFermion quark_sink(map_obj[key_sink]);

	      for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
	      {
		SerialDBKey<KeyPropElementalOperator_t> key;
		key.key().t_source     = t_source;
		key.key().colorvec_src = colorvec_source;
		key.key().colorvec_snk = colorvec_sink;
		key.key().spin_src     = spin_source;
		key.key().spin_snk     = spin_sink;
		key.key().mass_label   = params.param.mass_label;

		keyDispColorVector[0].colvec = colorvec_sink;
		keyDispColorVector[1].colvec = colorvec_source;

		watch.reset();
		watch.start();

		// Contract over color indices
		// Do the relevant quark contraction
		// Slow fourier-transform
		LatticeComplex lop = localInnerProduct(smrd_disp_vecs.getDispVector(keyDispColorVector[0]),
						       smrd_disp_vecs.getDispVector(keyDispColorVector[1]));

		multi2d<ComplexD> op_sum = phases.sft(lop);

		watch.stop();
		/*
		  QDPIO::cout << "Spatial Sums completed: time " << 
		  watch.getTimeInSeconds() << "secs" << endl;
		*/

		// Write the momentum projected fields
		for(int mom_num = 0 ; mom_num < phases.numMom() ; ++mom_num) 
		{
		  SerialDBKey<KeyGenPropElementalOperator_t> key;
		  key.key().colvec_l      = i;
		  key.key().colvec_r      = j;
		  key.key().mom           = phases.numToMom(mom_num);
		  key.key().displacement  = displacement_list[l]; // only right colorvector

		  SerialDBData<ValGenPropElementalOperator_t> val;
		  val.data().op           = op_sum[mom_num];
		  
		  // Insert into the DB
		  qdp_db.insert(key, val);

//	          write(xml_out, "elem", key.key());  // debugging
		} // mom_num
	      } // spin_sink
	    } // colorvec_sink
	  } // spin_source
	} // colorvec_source
	swiss.stop();

	QDPIO::cout << "Meson operator= " << l 
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

	push(file_xml, "GenPropElementalOperators");
	write(file_xml, "Params", params.param);
	write(file_xml, "Config_info", gauge_xml);
	write(file_xml, "Op_Info",displacement_list);
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
