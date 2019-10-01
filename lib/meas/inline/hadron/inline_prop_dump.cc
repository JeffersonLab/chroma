/*! \file
 * \brief Compute propagators from distillation
 *
 * Propagator calculation in distillation
 */

#include "qdp.h"
#include "fermact.h"
#include "meas/inline/hadron/inline_prop_dump.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "qdp_map_obj.h"
#include "qdp_map_obj_disk.h"
#include "qdp_map_obj_disk_multiple.h"
#include "qdp_map_obj_memory.h"
#include "qdp_disk_map_slice.h"
#include "util/ferm/subset_vectors.h"
#include "util/ferm/key_prop_distillation.h"
#include "util/ferm/key_timeslice_colorvec.h"
#include "util/ferm/key_prop_colorvec.h"
#include "util/ferm/key_prop_matelem.h"
#include "util/ferm/key_val_db.h"
#include "util/ferm/transf.h"
#include "util/ferm/spin_rep.h"
#include "util/ferm/diractodr.h"
#include "util/ferm/twoquark_contract_ops.h"
#include "util/ft/sftmom.h"
#include "util/ft/time_slice_set.h"
#include "util/info/proginfo.h"
#include "util/info/proginfo.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

#include "chroma_config.h"

#ifndef QDP_IS_QDPJIT_NO_NVPTX






namespace Chroma 
{ 

  //----------------------------------------------------------------------------
  namespace InlinePropDumpEnv 
  {
    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "prop_op_file", input.prop_op_file);
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "prop_op_file", input.prop_op_file);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::Param_t& input)
    {
      XMLReader inputtop(xml, path);

      int           t_source;      /*!< Array of time slice sources for props */
      std::string   mass_label;     /*!< Some kind of mass label */
      int           num_vecs;      /*!< In case of bad things happening in the solution vectors, do retries */

      read(inputtop, "num_vecs", input.num_vecs);
      read(inputtop, "t_source", input.t_source);
      read(inputtop, "Nt_forward", input.Nt_forward);
      read(inputtop, "mass_label", input.mass_label);
      read(inputtop, "decay_dir", input.decay_dir);
    }


    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t& input)
    {
      push(xml, path);

      write(xml, "num_vecs", input.num_vecs);
      write(xml, "t_source", input.t_source);
      write(xml, "mass_label", input.mass_label);
      write(xml, "decay_dir", input.decay_dir);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params& input)
    {
      Params tmp(xml, path);
      input = tmp;
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params& input)
    {
      push(xml, path);
    
      write(xml, "Param", input.param);
      write(xml, "NamedObject", input.named_obj);

      pop(xml);
    }
  } // namespace InlinePropDistillationEnv 


  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  namespace InlinePropDumpEnv 
  {
      //----------------------------------------------------------------------------
      //! Get active time-slices
      std::vector<bool> getActiveTSlices(int t_source, int Nt_forward, int Nt_backward)
      {
	// Initialize the active time slices
	const int decay_dir = Nd-1;
	const int Lt = Layout::lattSize()[decay_dir];

	std::vector<bool> active_t_slices(Lt);
	for(int t=0; t < Lt; ++t)
	{
	  active_t_slices[t] = false;
	}

	// Forward
	for(int dt=0; dt < Nt_forward; ++dt)
	{
	  int t = t_source + dt;
	  active_t_slices[t % Lt] = true;
	}

	// Backward
	for(int dt=0; dt < Nt_backward; ++dt)
	{
	  int t = t_source - dt;
	  while (t < 0) {t += Lt;} 

	  active_t_slices[t % Lt] = true;
	}

	return active_t_slices;
      }


    //----------------------------------------------------------------------------
    //! Get sink keys
    std::list<KeyPropElementalOperator_t> getSnkKeys(int t_source, int spin_source, int Nt_forward, int Nt_backward, const std::string mass)
    {
      std::list<KeyPropElementalOperator_t> keys;

      std::vector<bool> active_t_slices = getActiveTSlices(t_source, Nt_forward, Nt_backward);
	
      const int decay_dir = Nd-1;
      const int Lt = Layout::lattSize()[decay_dir];
	
      for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
	{
	  for(int t=0; t < Lt; ++t)
	    {
	      if (! active_t_slices[t]) {continue;}

	      KeyPropElementalOperator_t key;

	      key.t_source     = t_source;
	      key.t_slice      = t;
	      key.spin_src     = spin_source;
	      key.spin_snk     = spin_sink;
	      key.mass_label   = mass;

	      keys.push_back(key);
	    } // for t
	} // for spin_sink

      return keys;
    }
	
  } // end anonymous


	
  //----------------------------------------------------------------------------
  namespace InlinePropDumpEnv 
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
      
    const std::string name = "PROP_DUMP";

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
	  QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
	  QDP_abort(1);
	}
    }


    namespace 
    {
      std::ostream& operator<< (std::ostream& stream, const KeyPropElementalOperator_t& k) {
	stream << "t_slice = " << k.t_slice 
	       << "\nt_source = " << k.t_source 
	       << "\nspin_src = " << k.spin_src
	       << "\nspin_snk = " << k.spin_snk
	       << "\nmass_label = " << k.mass_label
	       << "\n";
	return stream;
      }
    }
  

    //----------------------------------------------------------------------------
    //----------------------------------------------------------------------------
    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
	{
	  std::string xml_file = makeXMLFileName(params.xml_file, update_no);

	  push(xml_out, "PropDistillation");
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

      proginfo(xml_out);    // Print out basic program info

      // Write out the input
      write(xml_out, "Input", params);


      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);


      // Will use TimeSliceSet-s a lot
      const int decay_dir = params.param.decay_dir;
      const int t_source  = params.param.t_source;
      const int Lt        = Layout::lattSize()[decay_dir];
	

      // A sanity check
      if (decay_dir != Nd-1)
	{
	  QDPIO::cerr << name << ": TimeSliceIO only supports decay_dir= " << Nd-1 << "\n";
	  QDP_abort(1);
	}


      //
      // DB storage
      //
      BinaryStoreDB< SerialDBKey<KeyPropElementalOperator_t>, SerialDBData<ValPropElementalOperator_t> > qdp_db;

      // Open the file, and write the meta-data and the binary for this operator
      if (! qdp_db.fileExists(params.named_obj.prop_op_file))
	{
	  QDP_error_exit("peram file not found");
	}
      else
	{
	  qdp_db.open(params.named_obj.prop_op_file, O_RDONLY, 0664);
	}

      QDPIO::cout << "Finished opening peram file" << std::endl;


	
      QDPIO::cout << "t_source = " << t_source << std::endl; 


      for(int spin_source=0; spin_source < Ns; ++spin_source)
	{
	  QDPIO::cout << "spin_source = " << spin_source << std::endl; 

	  // These are the common parts of a perambulator that are needed for this time source
	  std::list<KeyPropElementalOperator_t> snk_keys(getSnkKeys(t_source,
								    spin_source,
								    params.param.Nt_forward,
								    0,
								    params.param.mass_label));


	  // Initialize
	  for(std::list<KeyPropElementalOperator_t>::const_iterator key = snk_keys.begin();
	      key != snk_keys.end();
	      ++key)
	    {
	      QDPIO::cout << "t_slice    = " << key->t_slice
			  << ", t_source   = " << key->t_source
			  << ", spin_src   = " << key->spin_src
			  << ", spin_snk   = " << key->spin_snk
			  << ", mass_label = " << key->mass_label << "\n";

	      ValPropElementalOperator_t tmp;
	      SerialDBData<ValPropElementalOperator_t> serial(tmp);

	      if (qdp_db.get(*key, serial ) != 0)
		{
		  std::vector< SerialDBKey<KeyPropElementalOperator_t> > keys;
		  qdp_db.keys( keys );

		  for ( auto& k : keys )
		    std::cout << k.key() << "\n";

		  QDP_error_exit("key not found");
		}
	      
	      for (int n = 0 ; n < serial.data().mat.size1() ; ++n )
		for (int m = 0 ; m < serial.data().mat.size2() ; ++m )
		  QDPIO::cout << "(" << n << "," << m << ") = " << serial.data().mat(n,m) << "\n";
	    }


	} // for spin_src



      pop(xml_out);

      //pop(xml_out);  // prop_dist


      QDPIO::cout << name << ": ran successfully" << std::endl;

      END_CODE();
    }

  }

} // namespace Chroma

#endif
