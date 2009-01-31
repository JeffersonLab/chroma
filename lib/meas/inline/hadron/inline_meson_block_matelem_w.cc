// $Id: inline_meson_block_matelem_w.cc,v 3.2 2009-01-31 15:29:07 kostas Exp $
/*! \file
 * \brief Inline measurement of meson operators via colorvector matrix elements
 */

#include "handle.h"
#include "meas/inline/hadron/inline_meson_block_matelem_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/displace.h"
#include "meas/glue/mesplq.h"
#include "util/ferm/block_subset.h"
#include "util/ferm/block_couplings.h"
#include "util/ferm/subset_vectors.h"
#include "util/ferm/key_val_db.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

#define COLORVEC_MATELEM_TYPE_ZERO       0
#define COLORVEC_MATELEM_TYPE_ONE        1
#define COLORVEC_MATELEM_TYPE_MONE       -1
#define COLORVEC_MATELEM_TYPE_GENERIC    10

namespace Chroma 
{ 
  /*!
   * \ingroup hadron
   *
   * @{
   */
  namespace InlineMesonBlockMatElemEnv 
  { 
    // Reader for input parameters
    void read(XMLReader& xml, const string& path, InlineMesonBlockMatElemEnv::Params::Param_t& param)
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
      read(paramtop, "orthog_basis", param.orthog_basis);
      read(paramtop, "block", param.block);

      param.link_smearing  = readXMLGroup(paramtop, "LinkSmearing", "LinkSmearingType");
    }


    // Writer for input parameters
    void write(XMLWriter& xml, const string& path, const InlineMesonBlockMatElemEnv::Params::Param_t& param)
    {
      push(xml, path);

      int version = 1;

      write(xml, "version", version);
      write(xml, "mom2_max", param.mom2_max);
      write(xml, "displacement_length", param.displacement_length);
      write(xml, "displacement_list", param.displacement_list);
      write(xml, "num_vecs", param.num_vecs);
      write(xml, "decay_dir", param.decay_dir);
      write(xml, "orthog_basis", param.orthog_basis);
      write(xml, "block", param.block);
     xml << param.link_smearing.xml;

      pop(xml);
    }

    //! Read named objects 
    void read(XMLReader& xml, const string& path, InlineMesonBlockMatElemEnv::Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "colorvec_id", input.colorvec_id);
      read(inputtop, "meson_op_file", input.meson_op_file);
    }

    //! Write named objects
    void write(XMLWriter& xml, const string& path, const InlineMesonBlockMatElemEnv::Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "colorvec_id", input.colorvec_id);
      write(xml, "meson_op_file", input.meson_op_file);

      pop(xml);
    }

    // Writer for input parameters
    void write(XMLWriter& xml, const string& path, const InlineMesonBlockMatElemEnv::Params& param)
    {
      param.writeXML(xml, path);
    }
  }


  namespace InlineMesonBlockMatElemEnv 
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

    const std::string name = "MESON_BLOCK_MATELEM";

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
    struct KeyMesonElementalOperator_t
    {
      int                t_slice;      /*!< Meson operator time slice */
      multi1d<int>       displacement; /*!< Displacement dirs of right colorvector */
      multi1d<int>       mom;          /*!< D-1 momentum of this operator */
      int b_left ;  // the block of the anti-quark
      int b_right ; // the block for the quark
    };

    //! Meson operator
    struct ValMesonElementalOperator_t
    {
      int                type_of_data; /*!< Flag indicating type of data (maybe trivial) */
      multi2d<ComplexD>  op;          /*!< Colorvector source and sink with momentum projection */
    };


    //----------------------------------------------------------------------------
    //! Holds key and value as temporaries
    struct KeyValMesonElementalOperator_t
    {
      SerialDBKey<KeyMesonElementalOperator_t>  key;
      SerialDBData<ValMesonElementalOperator_t> val;
    };


    //----------------------------------------------------------------------------
    //! KeyMesonElementalOperator reader
    void read(BinaryReader& bin, KeyMesonElementalOperator_t& param)
    {
      read(bin, param.t_slice);
      read(bin, param.displacement);
      read(bin, param.mom);
      read(bin, param.b_left);
      read(bin, param.b_right);
    }

    //! MesonElementalOperator write
    void write(BinaryWriter& bin, const KeyMesonElementalOperator_t& param)
    {
      write(bin, param.t_slice);
      write(bin, param.displacement);
      write(bin, param.mom);
      write(bin, param.b_left);
      write(bin, param.b_right);
    }

    //! MesonElementalOperator reader
    void read(XMLReader& xml, const std::string& path, KeyMesonElementalOperator_t& param)
    {
      XMLReader paramtop(xml, path);
    
      read(paramtop, "t_slice", param.t_slice);
      read(paramtop, "displacement", param.displacement);
      read(paramtop, "mom", param.mom);
      read(paramtop, "left_block", param.b_left);
      read(paramtop, "right_block", param.b_right);
    }

    //! MesonElementalOperator writer
    void write(XMLWriter& xml, const std::string& path, const KeyMesonElementalOperator_t& param)
    {
      push(xml, path);

      write(xml, "t_slice", param.t_slice);
      write(xml, "displacement", param.displacement);
      write(xml, "mom", param.mom);
      write(xml, "left_block", param.b_left);
      write(xml, "right_block", param.b_right);
      pop(xml);
    }


    //----------------------------------------------------------------------------
    //! MesonElementalOperator reader
    void read(BinaryReader& bin, ValMesonElementalOperator_t& param)
    {
      read(bin, param.type_of_data);

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

    //! MesonElementalOperator write
    void write(BinaryWriter& bin, const ValMesonElementalOperator_t& param)
    {
      write(bin, param.type_of_data);

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

	push(xml_out, "MesonMatElemColorVec");
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

	TheNamedObjMap::Instance().getData< SubsetVectors<LatticeColorVector> >(params.named_obj.colorvec_id).getEvectors();
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

      const SubsetVectors<LatticeColorVector>& eigen_source = 
	TheNamedObjMap::Instance().getData< SubsetVectors<LatticeColorVector> >(params.named_obj.colorvec_id);

      push(xml_out, "MesonMatElemColorVec");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": Meson color-vector matrix element" << endl;

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

      // Keep track of no displacements and zero momentum
      multi1d<int> no_displacement;
      multi1d<int> zero_mom(3); zero_mom = 0;

      //
      // Meson operators
      //
      // The creation and annihilation operators are the same without the
      // spin matrices.
      //
      QDPIO::cout << "Building meson operators" << endl;

      // DB storage
      BinaryFxStoreDB< SerialDBKey<KeyMesonElementalOperator_t>, SerialDBData<ValMesonElementalOperator_t> > 
	qdp_db(params.named_obj.meson_op_file, 10*1024*1024, 64*1024);

      push(xml_out, "ElementalOps");

      // Loop over all time slices for the source. This is the same 
      // as the subsets for  phases

      //Make the block Set                             
      Set blocks ;
      blocks.make(BlockFunc(params.param.decay_dir, params.param.block));

      // Loop over each operator 
      for(int l=0; l < params.param.displacement_list.size(); ++l)
      {
	StopWatch watch;

	QDPIO::cout << "Elemental operator: op = " << l << endl;

	// Make sure displacement is something sensible
	multi1d<int> disp = normDisp(params.param.displacement_list[l]);

	QDPIO::cout << "displacement = " << disp << endl;

	// Build the operator
	swiss.reset();
	swiss.start();
	for(int b(0);b<blocks.numSubsets();b++){
	  int blk_left = b ;
	  vector<int> blk_couplings   ; 
	  // call a routine that calculate the couplings
	  blk_couplings = block_couplings(b, blocks, disp, params.param.displacement_length);
	  for(int blk_c(0);blk_c< blk_couplings.size();blk_c++){
	    int blk_right = blk_couplings[blk_c] ;
	    // Big loop over the momentum projection
	    for(int mom_num = 0 ; mom_num < phases.numMom() ; ++mom_num) 
	      {
		// The keys for the spin and displacements for this particular elemental operator
		// No displacement for left colorvector, only displace right colorvector
		// Invert the time - make it an independent key
		multi1d<KeyValMesonElementalOperator_t> buf(phases.numSubsets());
		for(int t=0; t < phases.numSubsets(); ++t)
		  {
		    buf[t].key.key().t_slice       = t;
		    buf[t].key.key().mom           = phases.numToMom(mom_num);
		    buf[t].key.key().displacement  = disp; // only right colorvector
		    buf[t].key.key().b_left        = blk_left ;
		    buf[t].key.key().b_right       = blk_right ;
		    buf[t].val.data().op.resize(params.param.num_vecs,params.param.num_vecs);
		    
		    if ( params.param.orthog_basis && 
			 (phases.numToMom(mom_num)) == zero_mom && 
			 (disp == no_displacement) )
		      {
			buf[t].val.data().type_of_data = COLORVEC_MATELEM_TYPE_ONE;
		      }
		    else
		      {
			buf[t].val.data().type_of_data = COLORVEC_MATELEM_TYPE_GENERIC;
		      }
		  }

		for(int j = 0 ; j < params.param.num_vecs; ++j)
		  {
		    // Displace the right vector and multiply by the momentum phase
		    LatticeColorVector tt = zero ;
		    tt[blocks[blk_right]] = eigen_source.getEvectors()[j] ;
		    LatticeColorVector shift_vec = phases[mom_num] * 
		      displace(u_smr, tt, params.param.displacement_length, disp);

		    for(int i = 0 ; i <  params.param.num_vecs; ++i)
		      {
			watch.reset();
			watch.start();
			
			tt=zero ;
			tt[blocks[blk_left]] = eigen_source.getEvectors()[i] ;

			// Contract over color indices
			// Do the relevant quark contraction
			LatticeComplex lop = localInnerProduct(tt, shift_vec);
			
			// Slow fourier-transform
			multi1d<ComplexD> op_sum = sumMulti(lop, phases.getSet());
			
			watch.stop();
			
			for(int t=0; t < op_sum.size(); ++t)
			  {
			    buf[t].val.data().op(i,j) = op_sum[t];
			  }

			//  write(xml_out, "elem", key.key());  // debugging
		      } // end for j
		  } // end for i

		QDPIO::cout << "insert: mom= " << phases.numToMom(mom_num) << " displacement= " << disp << endl; 
		for(int t=0; t < phases.numSubsets(); ++t)
		  {
		    qdp_db.insert(buf[t].key, buf[t].val);
		  }
		
	      } // mom_num

	  }// block couplings loop
	}// loop over blocks
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

	push(file_xml, "MesonElementalOperators");
	write(file_xml, "lattSize", QDP::Layout::lattSize());
	write(file_xml, "decay_dir", params.param.decay_dir);
	write(file_xml, "Weights", eigen_source.getEvalues());
	write(file_xml, "Params", params.param);
	write(file_xml, "Config_info", gauge_xml);
	write(file_xml, "Op_Info",params.param.displacement_list);
	pop(file_xml);

	qdp_db.insertUserdata(file_xml.str());
      }
      swiss.stop();

      QDPIO::cout << "Meson Operator written:"
		  << "  time= " << swiss.getTimeInSeconds() << " secs" << endl;

      // Close the namelist output file XMLDAT
      pop(xml_out);     // MesonMatElemColorVector

      snoop.stop();
      QDPIO::cout << name << ": total time = " 
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << name << ": ran successfully" << endl;

      END_CODE();
    } // func
  } // namespace InlineMesonBlockMatElemEnv

  /*! @} */  // end of group hadron

} // namespace Chroma
