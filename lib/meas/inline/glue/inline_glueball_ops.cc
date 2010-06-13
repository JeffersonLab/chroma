/*! \file
 * \brief Inline measurement of glueball operators
 */


#include "handle.h"
#include "meas/inline/glue/inline_glueball_ops.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/displace.h"
#include "meas/glue/mesplq.h"
#include "meas/glue/mesfield.h"
#include "util/ferm/subset_vectors.h"
#include "util/ferm/key_val_db.h"
#include "util/gauge/key_glue_matelem.h"
#include "util/gauge/taproj.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  /*!
   * \ingroup inlineglue
   *
   * @{
   */
  namespace InlineGlueballOpsEnv 
  { 
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

      read(paramtop, "mom2_max", param.mom2_max);
      read(paramtop, "displacement_length", param.displacement_length);
      read(paramtop, "displacement_list", param.displacement_list);
      read(paramtop, "decay_dir", param.decay_dir);

      param.link_smearing  = readXMLGroup(paramtop, "LinkSmearing", "LinkSmearingType");
    }


    // Writer for input parameters
    void write(XMLWriter& xml, const string& path, const Params::Param_t& param)
    {
      push(xml, path);

      int version = 1;

      write(xml, "version", version);
      write(xml, "mom2_max", param.mom2_max);
      write(xml, "displacement_length", param.displacement_length);
      write(xml, "displacement_list", param.displacement_list);
      write(xml, "decay_dir", param.decay_dir);
      xml << param.link_smearing.xml;

      pop(xml);
    }

    //! Read named objects 
    void read(XMLReader& xml, const string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "glue_op_file", input.glue_op_file);
    }

    //! Write named objects
    void write(XMLWriter& xml, const string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "glue_op_file", input.glue_op_file);

      pop(xml);
    }

    // Writer for input parameters
    void write(XMLWriter& xml, const string& path, const Params& param)
    {
      param.writeXML(xml, path);
    }
  }


  namespace InlineGlueballOpsEnv 
  { 
    // Anonymous namespace for registration
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMeas(Params(xml_in, path));
      }

      const std::string name = "GLUEBALL_OPS";

      //! Local registration flag
      bool registered = false;
    }

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

	push(xml_out, "GlueballOps");
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
      try
      {
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
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

      // Cast should be valid now
      const multi1d<LatticeColorMatrix>& u = 
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

      push(xml_out, "GlueballOps");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": Glueball operators" << endl;

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
      // Hack for the moment. Can only support 0-momentum. For non-zero momentum, 
      // we need leftNabla-s that have proper momentum.
      //
      if (params.param.mom2_max != 0)
      {
	QDPIO::cerr << name << ": only support zero momentum at the moment. Need generalizatin for left derivs\n";
	QDP_abort(1);
      }

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

      // Record the smeared observables
      MesPlq(xml_out, "Smeared_Observables", u_smr);


      //
      // Construct chroma-mag fields
      //
      multi1d<LatticeColorMatrix> B_mag(Nd-1); 

      try
      {
	if (Nd != 4)
	{
	  QDPIO::cerr << name << ": expects a 4d build\n";
	  QDP_abort(1);
	}

	// Computes anti-hermitian F fields
	multi1d<LatticeColorMatrix> f;
	mesField(f, u_smr);

	// The B fields
	multi1d<int> ind(3);
	ind[0] = 3;
	ind[1] = 1;
	ind[2] = 0;

	multi1d<int> sgnn(3);
	sgnn[0] = +1;
	sgnn[1] = -1;
	sgnn[2] = +1;
	
	// Already anti-hermitan. Remove trace and convert to hermitian to keep my sanity.
	for(int i=0; i < B_mag.size(); ++i)
	{
	  taproj(f[ind[i]]);
	  B_mag[i] = sgnn[i] * imag(f[ind[i]]);
	}
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << name << ": some exception while producing B fields: " << e << endl;
	QDP_abort(1);
      }
      

      //
      // DB storage
      //
      BinaryStoreDB< SerialDBKey<KeyGlueElementalOperator_t>, SerialDBData<ValGlueElementalOperator_t> > 
	qdp_db;

      // Open the file, and write the meta-data and the binary for this operator
      if (! qdp_db.fileExists(params.named_obj.glue_op_file))
      {
	XMLBufferWriter file_xml;

	// A hacked version of eigenvalues
	multi1d<SubsetVectorWeight_t> evals(1);
	evals[0].weights.resize(QDP::Layout::lattSize()[params.param.decay_dir]);
	evals[0].weights = Real(1);

	push(file_xml, "DBMetaData");
	write(file_xml, "id", string("glueElemOp"));
	write(file_xml, "lattSize", QDP::Layout::lattSize());
	write(file_xml, "decay_dir", params.param.decay_dir);
	proginfo(file_xml);    // Print out basic program info
	write(file_xml, "Params", params.param);
	write(file_xml, "Op_Info", params.param.displacement_list);
	write(file_xml, "Config_info", gauge_xml);
	write(file_xml, "Weights", evals);
	pop(file_xml);

	std::string file_str(file_xml.str());
	qdp_db.setMaxUserInfoLen(file_str.size());

	qdp_db.open(params.named_obj.glue_op_file, O_RDWR | O_CREAT, 0664);

	qdp_db.insertUserdata(file_str);
      }
      else
      {
	qdp_db.open(params.named_obj.glue_op_file, O_RDWR, 0664);
      }


      //
      // Glue operators
      //
      // The creation and annihilation operators are the same without the
      // spin matrices.
      //
      QDPIO::cout << "Building glue operators" << endl;

      push(xml_out, "ElementalOps");


      // Loop over all time slices for the source. This is the same 
      // as the subsets for  phases

      // Loop over each operator 
      for(int l=0; l < params.param.displacement_list.size(); ++l)
      {
	swiss.reset();
	swiss.start();

	// Build the operator
	QDPIO::cout << "Elemental operator: op = " << l << endl;

	// Make sure displacement is something sensible
	multi1d<int> disp = normDisp(params.param.displacement_list[l]);

	QDPIO::cout << "displacement = " << disp << endl;

	// Right mag field
	for(int j = 0 ; j < B_mag.size(); ++j)
	{
	  // Displace the right vector
	  LatticeColorMatrix shift_vec = rightNabla(u_smr, 
						    B_mag[j],
						    params.param.displacement_length, 
						    disp);

	  // Left mag field
	  for(int i = 0 ; i < B_mag.size(); ++i)
	  {
	    // Contract over color indices.
	    // NOTE: no need for daggers - B fields are hermitian
	    LatticeComplex lop = trace(B_mag[i] * shift_vec);


	    // Big loop over the momentum projection
	    for(int mom_num = 0 ; mom_num < phases.numMom() ; ++mom_num) 
	    {
	      // Slow fourier-transform
	      multi1d<ComplexD> op_sum = sumMulti(phases[mom_num] * lop, phases.getSet());

	      // The keys for the spin and displacements for this particular elemental operator
	      // No displacement for left B-field, only displace right B-field
	      // Invert the time - make it an independent key.
	      KeyValGlueElementalOperator_t buf;

	      for(int t=0; t < phases.numSubsets(); ++t)
	      {
		buf.key.key().t_slice       = t;
		buf.key.key().mom           = phases.numToMom(mom_num);
		buf.key.key().left          = i + 1;
		buf.key.key().right         = j + 1;
		buf.key.key().displacement  = disp;

		// Should not be a phase
		buf.val.data().op.resize(1);
		buf.val.data().op(0) = op_sum[t];

//		QDPIO::cout << "insert: mom= " << phases.numToMom(mom_num) << " displacement= " << disp << endl; 
		qdp_db.insert(buf.key, buf.val);

	      } // end for t
	    } // end for i
	  } // end for j

	} // mom_num

	swiss.stop();

	QDPIO::cout << "Glue operator= " << l 
		    << "  time= "
		    << swiss.getTimeInSeconds() 
		    << " secs" << endl;

      } // for l

      pop(xml_out); // ElementalOps

      // Close the namelist output file XMLDAT
      pop(xml_out);     // GlueballOpstor

      snoop.stop();
      QDPIO::cout << name << ": total time = " 
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << name << ": ran successfully" << endl;

      END_CODE();
    } // func
  } // namespace InlineGlueballOpsEnv

  /*! @} */  // end of group inlineglue

} // namespace Chroma


