// $Id: inline_meson_matelem_colorvec_w.cc,v 1.2 2008-06-21 05:14:49 edwards Exp $
/*! \file
 * \brief Inline measurement of meson operators via colorvector matrix elements
 */

#include "handle.h"
#include "meas/inline/hadron/inline_meson_matelem_colorvec_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/smear/displacement.h"
#include "util/ferm/eigeninfo.h"
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
  namespace InlineMesonMatElemColorVecEnv 
  { 
  // Reader for input parameters
  void read(XMLReader& xml, const string& path, InlineMesonMatElemColorVecEnv::Params::Param_t& param)
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
    read(paramtop, "num_vecs", param.num_vecs);
    read(paramtop, "decay_dir", param.decay_dir);

    param.link_smearing  = readXMLGroup(paramtop, "LinkSmearing", "LinkSmearingType");
  }


  // Writer for input parameters
  void write(XMLWriter& xml, const string& path, const InlineMesonMatElemColorVecEnv::Params::Param_t& param)
  {
    push(xml, path);

    int version = 1;

    write(xml, "version", version);
    write(xml, "mom2_max", param.mom2_max);
    write(xml, "displacement_length", param.displacement_length);
    write(xml, "num_vecs", param.num_vecs);
    write(xml, "decay_dir", param.decay_dir);
    xml << param.link_smearing.xml;

    pop(xml);
  }

  //! Read named objects 
  void read(XMLReader& xml, const string& path, InlineMesonMatElemColorVecEnv::Params::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "colorvec_id", input.colorvec_id);
    read(inputtop, "displacement_file", input.displacement_file);
    read(inputtop, "meson_op_file", input.meson_op_file);
  }

  //! Write named objects
  void write(XMLWriter& xml, const string& path, const InlineMesonMatElemColorVecEnv::Params::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "colorvec_id", input.colorvec_id);
    write(xml, "displacement_file", input.displacement_file);
    write(xml, "meson_op_file", input.meson_op_file);

    pop(xml);
  }

  // Writer for input parameters
  void write(XMLWriter& xml, const string& path, const InlineMesonMatElemColorVecEnv::Params& param)
  {
    param.writeXML(xml, path);
  }
  }


  namespace InlineMesonMatElemColorVecEnv 
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



    //--------------------------------------------------------------
    //! 2-quark operator structure
    struct TwoQuarkOps_t
    {
      struct Displacement_t
      {
	multi1d<int> displacement;   /*!< Orig plus/minus 1-based directional displacements */
      };

      multi1d<Displacement_t> ops; /*!< 2-quark ops within a file */
    };

    //! Read two quark op
    void read(XMLReader& xml, const string& path, 
	      TwoQuarkOps_t::Displacement_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "displacement", input.displacement);
    }

    //! Read two quark ops
    void read(XMLReader& xml, const string& path, 
	      TwoQuarkOps_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "Operators", input.ops);
    }
	
    //! Write two quark op
    void write(XMLWriter& xml, const string& path, 
	       const TwoQuarkOps_t::Displacement_t& input)
    {
      push(xml, path);

      write(xml, "Displacement", input.displacement);

      pop(xml);
    }

    //! Write two quark op 
    void write(XMLWriter& xml, const string& path, 
	       const TwoQuarkOps_t& input)
    {
      push(xml, path);

      write(xml, "Operators", input.ops);

      pop(xml);
    }
	


    //----------------------------------------------------------------------------
    //! The key for smeared and displaced color vectors
    struct KeySmearedDispColorVector_t
    {
      int  colvec;                  /*!< Colorvector index */
      multi1d<int> displacement;    /*!< Orig plus/minus 1-based directional displacements */
    };


    //! Support for the keys of smeared and displaced color vectors
    bool operator<(const KeySmearedDispColorVector_t& a, const KeySmearedDispColorVector_t& b)
    {
      multi1d<int> lgaa(1);
      lgaa[0] = a.colvec;
      multi1d<int> lga = concat(lgaa, a.displacement);

      multi1d<int> lgbb(1);
      lgbb[0] = b.colvec;
      multi1d<int> lgb = concat(lgbb, b.displacement);

      return (lga < lgb);
    }


    //! The value of the map
    struct SmearedDispColorVector_t
    {
      LatticeColorVector vec;
    };


    //----------------------------------------------------------------------------
    //! The smeared and displaced objects
    class SmearedDispObjects
    {
    public:
      //! Constructor from smeared map 
      SmearedDispObjects(int disp_length,
			 const std::string& colorvec_id,
			 const multi1d<LatticeColorMatrix> & u_smr);

      //! Destructor
      ~SmearedDispObjects() {}

      //! Accessor
      const LatticeColorVector& getDispVector(const KeySmearedDispColorVector_t& key);

    protected:
      //! Displace an object
      const LatticeColorVector& displaceObject(const KeySmearedDispColorVector_t& key);
			
    private:
      //! Lattice color vectors
      const EigenInfo<LatticeColorVector>& eigen_source;

      //! Named object id of eigenvectors
      std::string colorvec_id;
		
      //! Gauge field 
      const multi1d<LatticeColorMatrix>& u;
			
      //! Displacement length
      int displacement_length;
			
      //!Maps of smeared displaced color vectors 
      map<KeySmearedDispColorVector_t, SmearedDispColorVector_t> disp_src_map;
    };

	
    // Constructor from smeared map 
    SmearedDispObjects::SmearedDispObjects(int disp_length,
					   const std::string& colorvec_id,
					   const multi1d<LatticeColorMatrix>& u_smr) :
      displacement_length(disp_length), 
      eigen_source(TheNamedObjMap::Instance().getData< EigenInfo<LatticeColorVector> >(colorvec_id)), 
      u(u_smr)
    {
    }


    //! Accessor
    const LatticeColorVector&
    SmearedDispObjects::getDispVector(const KeySmearedDispColorVector_t& key)
    {
      //Check if any displacement is needed
      if (displacement_length == 0) 
      {
	return eigen_source.getEvectors()[key.colvec];
      }
      else
      {
	return displaceObject(key);
      }
    }


    //! Accessor
    const LatticeColorVector&
    SmearedDispObjects::displaceObject(const KeySmearedDispColorVector_t& key)
    {
      StopWatch snoop;

      // If no entry, then create a displaced version of the quark
      if (disp_src_map.find(key) == disp_src_map.end())
      {
	//	      cout << __func__ 
	//		   << ": n=" << n
	//		   << " l=" << l
	//		   << " i=" << i 
	//		   << " disp=" << term.quark[i].displacement
	//		   << " len=" << term.quark[i].disp_len
	//		   << " dir=" << term.quark[i].disp_dir
	//		   << endl;

	// Insert an empty entry and then modify it. This saves on
	// copying the data around
	{
	  SmearedDispColorVector_t disp_empty;

	  snoop.reset();
	  snoop.start();

	  disp_src_map.insert(std::make_pair(key, disp_empty));

	  snoop.stop();

	  QDPIO::cout<<"Inserted key in map: time = "<< snoop.getTimeInSeconds() << "secs"<<endl;

	  // Sanity check - the entry better be there
	  if (disp_src_map.find(key) == disp_src_map.end())
	  {
	    QDPIO::cerr << __func__ 
			<< ": internal error - could not insert empty key in map"
			<< endl;
	    QDP_abort(1);
	  }		      
	}

	// Modify the previous empty entry
	SmearedDispColorVector_t& disp_q = disp_src_map.find(key)->second;
	disp_q.vec = eigen_source.getEvectors()[key.colvec];

	snoop.reset();
	snoop.start();

	for(int i=0; i < key.displacement.size(); ++i)
	{
	  if (key.displacement[i] > 0)
	  {
	    int disp_dir = key.displacement[i] - 1;
	    int disp_len = displacement_length;
	    displacement(u, disp_q.vec, disp_len, disp_dir);
	  }
	  else if (key.displacement[i] < 0)
	  {
	    int disp_dir = -key.displacement[i] - 1;
	    int disp_len = -displacement_length;
	    displacement(u, disp_q.vec, disp_len, disp_dir);
	  }
	}

	snoop.stop();

	QDPIO::cout << "Displaced Vector:  Disp = "
		    << key.displacement <<" Time = "<<snoop.getTimeInSeconds() <<" sec"<<endl;

      } // if find in map

      snoop.reset();
      snoop.start();

      // The key now must exist in the map, so return the vector
      SmearedDispColorVector_t& disp_q = disp_src_map.find(key)->second;

      snoop.stop(); 

      QDPIO::cout << "Retrieved entry from map : time = "<< snoop.getTimeInSeconds() << "secs "<<endl;

      return disp_q.vec;
    }

    //----------------------------------------------------------------------------
    //! Meson operator
    struct MesonElementalOperator_t
    {
      int                colvec_l;     /*!< Left colorvector index */
      int                colvec_r;     /*!< Right colorvector index */
      multi1d<int>       displacement; /*!< Displacement dirs of right colorvector */
      multi1d<int>       mom;          /*!< D-1 momentum of this operator */
      multi1d<DComplex>  op;           /*!< Momentum projected operator */
    };


    //----------------------------------------------------------------------------
    //! MesonElementalOperator writer
    void read(BinaryReader& bin, MesonElementalOperator_t& param)
    {
      read(bin, param.colvec_l);
      read(bin, param.colvec_r);
      read(bin, param.displacement);
      read(bin, param.mom);
      read(bin, param.op);
    }

    //! MesonElementalOperator writer
    void write(BinaryWriter& bin, const MesonElementalOperator_t& param)
    {
      write(bin, param.colvec_l);
      write(bin, param.colvec_r);
      write(bin, param.displacement);
      write(bin, param.mom);
      write(bin, param.op);
    }

    //! MesonElementalOperator writer
    void read(XMLReader& xml, const std::string& path, MesonElementalOperator_t& param)
    {
      XMLReader paramtop(xml, path);
    
      read(paramtop, "colvec_l", param.colvec_l);
      read(paramtop, "colvec_r", param.colvec_r);
      read(paramtop, "displacement", param.displacement);
      read(paramtop, "mom", param.mom);
      read(paramtop, "op", param.op);
    }

    //! MesonElementalOperator writer
    void write(XMLWriter& xml, const std::string& path, const MesonElementalOperator_t& param)
    {
      push(xml, path);

      write(xml, "colvec_l", param.colvec_l);
      write(xml, "colvec_r", param.colvec_r);
      write(xml, "displacement", param.displacement);
      write(xml, "mom", param.mom);
      write(xml, "op", param.op);

      pop(xml);
    }


    //----------------------------------------------------------------------------
    //! Read 2-quark operators file, assign correct displacement length
    void readOps(TwoQuarkOps_t& oplist, 
		 const std::string& displacement_file)
    {
      START_CODE();

      TextFileReader reader(displacement_file);

      int num_ops;
      reader >> num_ops;
      oplist.ops.resize(num_ops);

      //Loop over ops within a file
      for(int n=0; n < oplist.ops.size(); ++n)
      {
	TwoQuarkOps_t::Displacement_t& qq = oplist.ops[n];

	// Read 1-based displacement only for the right quark
	int ndisp;
	reader >> ndisp;
	multi1d<int> displacement(ndisp);
	for(int i=0; i < ndisp; ++i)
	  reader >> displacement[i];

	qq.displacement = displacement;
      } //n

      reader.close();

      END_CODE();
    } //void readOps


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

	push(xml_out, "stoch_meson");
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
      const multi1d<LatticeColorMatrix>& u = 
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

      push(xml_out, "MesonMatElemColorVec");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": Stochastic Meson Operator" << endl;

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
      // Read operator coefficients
      //
      QDPIO::cout << "Reading 2-quark operators" << endl;
      TwoQuarkOps_t qq_oplist; 

      readOps(qq_oplist, params.named_obj.displacement_file);

      //
      // The object holding the displaced color vector maps  
      //
      SmearedDispObjects smrd_disp_vecs(params.param.displacement_length,
					params.named_obj.colorvec_id, 
					u_smr);

      //
      // Meson operators
      //
      // The creation and annihilation operators are the same without the
      // spin matrices.
      //
      QDPIO::cout << "Building meson operators" << endl;

      BinaryBufferWriter src_record_bin;

      write(src_record_bin, params.param.mom2_max);
      write(src_record_bin, params.param.decay_dir);
      write(src_record_bin, params.param.num_vecs);
      write(src_record_bin, qq_oplist.ops.size());

      push(xml_out, "ElementalOps");

      // Loop over all time slices for the source. This is the same 
      // as the subsets for  phases

      // Loop over each operator 
      int total_num_elem = 0;
      for(int l=0; l < qq_oplist.ops.size(); ++l)
      {
	StopWatch watch;

	QDPIO::cout << "Elemental operator: op = " << l << endl;

	// Build the operator
	swiss.reset();
	swiss.start();

	// The keys for the spin and displacements for this particular elemental operator
	multi1d<KeySmearedDispColorVector_t> keySmearedDispColorVector(2);

	// No displacement for left colorvector, only displace right colorvector
	keySmearedDispColorVector[0].displacement.resize(1);
	keySmearedDispColorVector[0].displacement = 0;
	keySmearedDispColorVector[1].displacement = qq_oplist.ops[l].displacement;

	for(int i = 0 ; i <  params.param.num_vecs; ++i)
	{
	  for(int j = 0 ; j < params.param.num_vecs; ++j)
	  {
	    keySmearedDispColorVector[0].colvec = i;
	    keySmearedDispColorVector[1].colvec = j;

	    watch.reset();
	    watch.start();

	    // Contract over color indices
	    // Do the relevant quark contraction
	    // Slow fourier-transform
	    LatticeComplex lop = localInnerProduct(smrd_disp_vecs.getDispVector(keySmearedDispColorVector[0]),
						   smrd_disp_vecs.getDispVector(keySmearedDispColorVector[1]));

	    multi2d<DComplex> op_sum = phases.sft(lop);

	    watch.stop();
	    /*
	      QDPIO::cout << "Spatial Sums completed: time " << 
	      watch.getTimeInSeconds() << "secs" << endl;
	    */		

	    // Write the momentum projected fields
	    for(int mom_num = 0 ; mom_num < phases.numMom() ; ++mom_num) 
	    {
	      MesonElementalOperator_t mop;
	      mop.colvec_l     = i;
	      mop.colvec_r     = j;
	      mop.mom          = phases.numToMom(mom_num);
	      mop.displacement = qq_oplist.ops[l].displacement; // only right colorvector
	      mop.op           = op_sum[mom_num];

	      write(src_record_bin, mop);
	      write(xml_out, "elem", mop);  // debugging
	      ++total_num_elem;
	    }
	  } // end for j
	} // end for i
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
	XMLBufferWriter src_record_xml, file_xml;

	push(file_xml, "MesonElementalOperators");
	write(file_xml, "Params", params.param);
	write(file_xml, "Config_info", gauge_xml);
	write(file_xml, "Op_Info",qq_oplist);
	pop(file_xml);

	QDPFileWriter qdp_file(file_xml, params.named_obj.meson_op_file,
			       QDPIO_SINGLEFILE, QDPIO_SERIAL, QDPIO_OPEN);

	push(src_record_xml, "MesonElementalOperator");
	write(src_record_xml, "mom2_max", params.param.mom2_max);
	write(src_record_xml, "decay_dir", params.param.decay_dir);
	write(src_record_xml, "num_vecs", params.param.num_vecs);
	write(src_record_xml, "num_disp", qq_oplist.ops.size());
	write(src_record_xml, "total_num_elem", total_num_elem);
	pop(src_record_xml);

	write(qdp_file, src_record_xml, src_record_bin);
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
  } // namespace InlineMesonMatElemColorVecEnv

  /*! @} */  // end of group hadron

} // namespace Chroma
