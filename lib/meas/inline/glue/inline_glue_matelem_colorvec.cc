/*! \file
 * \brief Inline measurement of glueball operators via colorvector matrix elements
 */


#include "handle.h"
#include "meas/inline/glue/inline_glue_matelem_colorvec.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/displace.h"
#include "meas/glue/mesplq.h"
#include "util/ferm/subset_vectors.h"
#include "util/ferm/key_val_db.h"
#include "util/gauge/key_glue_matelem.h"
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
  namespace InlineGlueMatElemColorVecEnv 
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
      read(paramtop, "num_vecs", param.num_vecs);
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
      write(xml, "num_vecs", param.num_vecs);
      write(xml, "decay_dir", param.decay_dir);
     xml << param.link_smearing.xml;

      pop(xml);
    }

    //! Read named objects 
    void read(XMLReader& xml, const string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "colorvec_id", input.colorvec_id);
      read(inputtop, "glue_op_file", input.glue_op_file);
    }

    //! Write named objects
    void write(XMLWriter& xml, const string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "colorvec_id", input.colorvec_id);
      write(xml, "glue_op_file", input.glue_op_file);

      pop(xml);
    }

    // Writer for input parameters
    void write(XMLWriter& xml, const string& path, const Params& param)
    {
      param.writeXML(xml, path);
    }
  }


  namespace InlineGlueMatElemColorVecEnv 
  { 
    // Anonymous namespace for registration
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMeas(Params(xml_in, path));
      }

      const std::string name = "GLUEBALL_MATELEM_COLORVEC";

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

	push(xml_out, "GlueMatElemColorVec");
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

	// NB We are just checking this is here.
	TheNamedObjMap::Instance().getData< Handle< MapObject<int,EVPair<LatticeColorVector> > > >(params.named_obj.colorvec_id);
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

      const MapObject<int,EVPair<LatticeColorVector> >& eigen_source = 
	*(TheNamedObjMap::Instance().getData< Handle< MapObject<int,EVPair<LatticeColorVector> > > >(params.named_obj.colorvec_id));

      push(xml_out, "GlueMatElemColorVec");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": Glue color-vector matrix element" << endl;

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

      // Record the smeared observables
      MesPlq(xml_out, "Smeared_Observables", u_smr);

      //
      // DB storage
      //
      BinaryStoreDB< SerialDBKey<KeyGlueElementalOperator_t>, SerialDBData<ValGlueElementalOperator_t> > 
	qdp_db;

      // Open the file, and write the meta-data and the binary for this operator
      if (! qdp_db.fileExists(params.named_obj.glue_op_file))
      {
	XMLBufferWriter file_xml;

	push(file_xml, "DBMetaData");
	write(file_xml, "id", string("glueElemOp"));
	write(file_xml, "lattSize", QDP::Layout::lattSize());
	write(file_xml, "decay_dir", params.param.decay_dir);
	proginfo(file_xml);    // Print out basic program info
	write(file_xml, "Params", params.param);
	write(file_xml, "Op_Info", params.param.displacement_list);
	write(file_xml, "Config_info", gauge_xml);
	write(file_xml, "Weights", getEigenValues(eigen_source, params.param.num_vecs));
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
	StopWatch watch;

	QDPIO::cout << "Elemental operator: op = " << l << endl;

	// Make sure displacement is something sensible
	multi1d<int> disp = normDisp(params.param.displacement_list[l]);

	QDPIO::cout << "displacement = " << disp << endl;

	// The first version of this code used just rightNabla instead of leftRightNabla. 
	// Account for the change in normalization. Namely, for each derivative, multiply by -1/2.
	Real fnorm = 1.0;
	for(int mu=0; mu < disp.size(); ++mu)
	  fnorm *= -Real(0.5);

	// Build the operator
	swiss.reset();
	swiss.start();

	// Big loop over the momentum projection
	for(int mom_num = 0 ; mom_num < phases.numMom() ; ++mom_num) 
	{
	  // The keys for the spin and displacements for this particular elemental operator
	  // Note: the "disp" is actually left-right derivatives. Apply them to the right vector.
	  // Invert the time - make it an independent key
	  multi1d<KeyValGlueElementalOperator_t> buf(phases.numSubsets());
	  for(int t=0; t < phases.numSubsets(); ++t)
	  {
	    buf[t].key.key().t_slice       = t;
	    buf[t].key.key().mom           = phases.numToMom(mom_num);
	    buf[t].key.key().left          = 0;    // not used
	    buf[t].key.key().right         = 0;    // not used
	    buf[t].key.key().displacement  = disp; // only right colorvector
	    buf[t].val.data().op.resize(params.param.num_vecs);
	  }

	  for(int j = 0 ; j < params.param.num_vecs; ++j)
	  {
	    watch.reset();
	    watch.start();
	    
	    // Apply a left-right derivative onto the right vector. The momentum is required.
	    // Multiply by the momentum phase to project onto a definite momentum.
	    EVPair<LatticeColorVector> tmpvec; eigen_source.get(j,tmpvec);
	    LatticeColorVector lvec(tmpvec.eigenVector);

	    LatticeColorVector shift_vec = leftRightNabla(u_smr, 
							  lvec,
							  params.param.displacement_length, 
							  disp,
							  phases.numToMom(mom_num));

	    // Contract over color indices
	    // Do the relevant quark contraction
	    LatticeComplex lop = localInnerProduct(lvec, shift_vec);

	    // Slow fourier-transform
	    multi1d<ComplexD> op_sum = sumMulti(phases[mom_num] * lop, phases.getSet());

	    watch.stop();

	    for(int t=0; t < op_sum.size(); ++t)
	    {
	      buf[t].val.data().op(j) = fnorm * op_sum[t];
	    }
	  } // end for j

//	  QDPIO::cout << "insert: mom= " << phases.numToMom(mom_num) << " displacement= " << disp << endl; 
	  for(int t=0; t < phases.numSubsets(); ++t)
	  {
	    qdp_db.insert(buf[t].key, buf[t].val);
	  }

	} // mom_num

	swiss.stop();

	QDPIO::cout << "Glue operator= " << l 
		    << "  time= "
		    << swiss.getTimeInSeconds() 
		    << " secs" << endl;

      } // for l

      pop(xml_out); // ElementalOps

      // Close the namelist output file XMLDAT
      pop(xml_out);     // GlueMatElemColorVector

      snoop.stop();
      QDPIO::cout << name << ": total time = " 
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << name << ": ran successfully" << endl;

      END_CODE();
    } // func
  } // namespace InlineGlueMatElemColorVecEnv

  /*! @} */  // end of group inlineglue

} // namespace Chroma


