// $Id: inline_baryon_matelem_colorvec_w.cc,v 3.2 2008-08-06 15:20:52 edwards Exp $
/*! \file
 * \brief Inline measurement of baryon operators via colorvector matrix elements
 */

#include "handle.h"
#include "meas/inline/hadron/inline_baryon_matelem_colorvec_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/smear/disp_colvec_map.h"
#include "util/ferm/eigeninfo.h"
#include "util/ferm/key_val_db.h"
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
  namespace InlineBaryonMatElemColorVecEnv 
  { 
    // Reader for input parameters
    void read(XMLReader& xml, const string& path, InlineBaryonMatElemColorVecEnv::Params::Param_t::Displacement_t& param)
    {
      XMLReader paramtop(xml, path);

      read(paramtop, "left", param.left);
      read(paramtop, "middle", param.middle);
      read(paramtop, "right", param.right);
    }


    // Writer for input parameters
    void write(XMLWriter& xml, const string& path, const InlineBaryonMatElemColorVecEnv::Params::Param_t::Displacement_t& param)
    {
      push(xml, path);

      write(xml, "left", param.left);
      write(xml, "middle", param.middle);
      write(xml, "right", param.right);

      pop(xml);
    }


    // Reader for input parameters
    void read(XMLReader& xml, const string& path, InlineBaryonMatElemColorVecEnv::Params::Param_t& param)
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
    void write(XMLWriter& xml, const string& path, const InlineBaryonMatElemColorVecEnv::Params::Param_t& param)
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
    void read(XMLReader& xml, const string& path, InlineBaryonMatElemColorVecEnv::Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "colorvec_id", input.colorvec_id);
      read(inputtop, "baryon_op_file", input.baryon_op_file);
    }

    //! Write named objects
    void write(XMLWriter& xml, const string& path, const InlineBaryonMatElemColorVecEnv::Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "colorvec_id", input.colorvec_id);
      write(xml, "baryon_op_file", input.baryon_op_file);

      pop(xml);
    }

    // Writer for input parameters
    void write(XMLWriter& xml, const string& path, const InlineBaryonMatElemColorVecEnv::Params& param)
    {
      param.writeXML(xml, path);
    }
  }


  namespace InlineBaryonMatElemColorVecEnv 
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

    const std::string name = "BARYON_MATELEM_COLORVEC";

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

      StandardOutputStream& operator<<(StandardOutputStream& os, const Params::Param_t::Displacement_t& d)
      {
	os << "left= " << d.left 
	   << " middle= " << d.middle
	   << " right= " << d.right;

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
    //! Baryon operator
    struct KeyBaryonElementalOperator_t
    {
      int                colvec_l;     /*!< Left colorvector index */
      int                colvec_m;     /*!< Middle colorvector index */
      int                colvec_r;     /*!< Right colorvector index */
      multi1d<int>       left;         /*!< Displacement dirs of left colorvector */
      multi1d<int>       middle;       /*!< Displacement dirs of middle colorvector */
      multi1d<int>       right;        /*!< Displacement dirs of right colorvector */
      multi1d<int>       mom;          /*!< D-1 momentum of this operator */
    };

    //! Baryon operator
    struct ValBaryonElementalOperator_t
    {
      multi1d<ComplexD>  op;           /*!< Momentum projected operator */
    };


    //----------------------------------------------------------------------------
    //! BaryonElementalOperator reader
    void read(BinaryReader& bin, KeyBaryonElementalOperator_t& param)
    {
      read(bin, param.colvec_l);
      read(bin, param.colvec_m);
      read(bin, param.colvec_r);
      read(bin, param.left);
      read(bin, param.middle);
      read(bin, param.right);
      read(bin, param.mom);
    }

    //! BaryonElementalOperator write
    void write(BinaryWriter& bin, const KeyBaryonElementalOperator_t& param)
    {
      write(bin, param.colvec_l);
      write(bin, param.colvec_m);
      write(bin, param.colvec_r);
      write(bin, param.left);
      write(bin, param.middle);
      write(bin, param.right);
      write(bin, param.mom);
    }

    //! BaryonElementalOperator reader
    void read(XMLReader& xml, const std::string& path, KeyBaryonElementalOperator_t& param)
    {
      XMLReader paramtop(xml, path);
    
      read(paramtop, "colvec_l", param.colvec_l);
      read(paramtop, "colvec_m", param.colvec_m);
      read(paramtop, "colvec_r", param.colvec_r);
      read(paramtop, "left", param.left);
      read(paramtop, "middle", param.middle);
      read(paramtop, "right", param.right);
      read(paramtop, "mom", param.mom);
    }

    //! BaryonElementalOperator writer
    void write(XMLWriter& xml, const std::string& path, const KeyBaryonElementalOperator_t& param)
    {
      push(xml, path);

      write(xml, "colvec_l", param.colvec_l);
      write(xml, "colvec_m", param.colvec_m);
      write(xml, "colvec_r", param.colvec_r);
      write(xml, "left", param.left);
      write(xml, "middle", param.middle);
      write(xml, "right", param.right);
      write(xml, "mom", param.mom);

      pop(xml);
    }


    //----------------------------------------------------------------------------
    //! BaryonElementalOperator reader
    void read(BinaryReader& bin, ValBaryonElementalOperator_t& param)
    {
      read(bin, param.op);
    }

    //! BaryonElementalOperator write
    void write(BinaryWriter& bin, const ValBaryonElementalOperator_t& param)
    {
      write(bin, param.op);
    }

    //! BaryonElementalOperator reader
    void read(XMLReader& xml, const std::string& path, ValBaryonElementalOperator_t& param)
    {
      XMLReader paramtop(xml, path);
    
      read(paramtop, "op", param.op);
    }

    //! BaryonElementalOperator writer
    void write(XMLWriter& xml, const std::string& path, const ValBaryonElementalOperator_t& param)
    {
      push(xml, path);

      write(xml, "op", param.op);

      pop(xml);
    }


    //----------------------------------------------------------------------------
    //! Make sure displacements are something sensible
    multi1d<Params::Param_t::Displacement_t> 
    normalizeDisplacements(const multi1d<Params::Param_t::Displacement_t>& orig_list)
    {
      START_CODE();

      multi1d<Params::Param_t::Displacement_t> displacement_list(orig_list.size());

      multi1d<int> empty(1); empty = 0;

      // Loop over displacements
      for(int n=0; n < orig_list.size(); ++n)
      {
	const Params::Param_t::Displacement_t& o = orig_list[n];
	Params::Param_t::Displacement_t&       d = displacement_list[n];

	if (o.left.size() == 0)
	  d.left = empty;
	else
	  d.left = o.left;

	if (o.middle.size() == 0)
	  d.middle = empty;
	else
	  d.middle = o.middle;

	if (o.right.size() == 0)
	  d.right = empty;
	else
	  d.right = o.right;

	QDPIO::cout << "disp[" << n << "]="
		    << "  left= " << d.left
		    << "  middle= " << d.middle
		    << "  right= " << d.right
		    << endl;
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

	push(xml_out, "BaryonMatElemColorVec");
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

	TheNamedObjMap::Instance().getData< EigenInfo<LatticeColorVector> >(params.named_obj.colorvec_id).getEvectors();
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
      const multi1d<LatticeColorMatrix>& u = 
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

      const EigenInfo<LatticeColorVector>& eigen_source = 
	TheNamedObjMap::Instance().getData< EigenInfo<LatticeColorVector> >(params.named_obj.colorvec_id);

      push(xml_out, "BaryonMatElemColorVec");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": Baryon color-vector matrix element" << endl;

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
      multi1d<Params::Param_t::Displacement_t> displacement_list(normalizeDisplacements(params.param.displacement_list));

      //
      // The object holding the displaced color vector maps  
      //
      DispColorVectorMap smrd_disp_vecs(params.param.displacement_length,
					u_smr,
					eigen_source.getEvectors());

      //
      // Baryon operators
      //
      // The creation and annihilation operators are the same without the
      // spin matrices.
      //
      QDPIO::cout << "Building baryon operators" << endl;

      // DB storage
      BinaryVarStoreDB< SerialDBKey<KeyBaryonElementalOperator_t>, SerialDBData<ValBaryonElementalOperator_t> > 
	qdp_db(params.named_obj.baryon_op_file);

      push(xml_out, "ElementalOps");

      // Loop over all time slices for the source. This is the same 
      // as the subsets for  phases

      // Loop over each operator 
      for(int l=0; l < displacement_list.size(); ++l)
      {
	StopWatch watch;

	QDPIO::cout << "Elemental operator: op = " << l << endl;

	QDPIO::cout << "displacement: " << displacement_list[l] << endl;

	// Build the operator
	swiss.reset();
	swiss.start();

	// The keys for the spin and displacements for this particular elemental operator
	multi1d<KeyDispColorVector_t> keyDispColorVector(3);

	// No displacement for left colorvector, only displace right colorvector
	keyDispColorVector[0].displacement = displacement_list[l].left;
	keyDispColorVector[1].displacement = displacement_list[l].middle;
	keyDispColorVector[2].displacement = displacement_list[l].right;

	for(int i = 0 ; i <  params.param.num_vecs; ++i)
	{
	  for(int j = 0 ; j < params.param.num_vecs; ++j)
	  {
	    for(int k = 0 ; k < params.param.num_vecs; ++k)
	    {
	      keyDispColorVector[0].colvec = i;
	      keyDispColorVector[1].colvec = j;
	      keyDispColorVector[2].colvec = k;

	      watch.reset();
	      watch.start();

	      // Contract over color indices
	      // Do the relevant quark contraction
	      // Slow fourier-transform
	      LatticeComplex lop = colorContract(smrd_disp_vecs.getDispVector(keyDispColorVector[0]),
						 smrd_disp_vecs.getDispVector(keyDispColorVector[1]),
						 smrd_disp_vecs.getDispVector(keyDispColorVector[2]));

	      multi2d<ComplexD> op_sum = phases.sft(lop);

	      watch.stop();
	      /*
		QDPIO::cout << "Spatial Sums completed: time " << 
		watch.getTimeInSeconds() << "secs" << endl;
	      */		

	      // Write the momentum projected fields
	      for(int mom_num = 0 ; mom_num < phases.numMom() ; ++mom_num) 
	      {
		SerialDBKey<KeyBaryonElementalOperator_t> key;
		key.key().colvec_l     = i;
		key.key().colvec_m     = j;
		key.key().colvec_r     = k;
		key.key().mom          = phases.numToMom(mom_num);
		key.key().left         = displacement_list[l].left;
		key.key().middle       = displacement_list[l].middle;
		key.key().right        = displacement_list[l].right;

		SerialDBData<ValBaryonElementalOperator_t> val;
		val.data().op          = op_sum[mom_num];

		qdp_db.insert(key, val);

//	        write(xml_out, "elem", mop);  // debugging
	      }
	    } // end for k
	  } // end for j
	} // end for i
	swiss.stop();

	QDPIO::cout << "Baryon operator= " << l 
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

	push(file_xml, "BaryonElementalOperators");
	write(file_xml, "Params", params.param);
	write(file_xml, "Config_info", gauge_xml);
	write(file_xml, "Op_Info",displacement_list);
	pop(file_xml);

	qdp_db.insertUserdata(file_xml.str());
      }
      swiss.stop();

      QDPIO::cout << "Baryon Operator written:"
		  << "  time= " << swiss.getTimeInSeconds() << " secs" << endl;

      // Close the namelist output file XMLDAT
      pop(xml_out);     // BaryonMatElemColorVector

      snoop.stop();
      QDPIO::cout << name << ": total time = " 
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << name << ": ran successfully" << endl;

      END_CODE();
    } // func
  } // namespace InlineBaryonMatElemColorVecEnv

  /*! @} */  // end of group hadron

} // namespace Chroma
