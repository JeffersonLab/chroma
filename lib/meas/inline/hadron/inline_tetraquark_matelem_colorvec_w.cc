#include "util/ferm/colorvec_mod_reader.h"
#include "util/ferm/superb_contractions.h"
/*! \file
 * \brief Inline measurement of tetraquark operators via colorvector matrix elements
 */

#include "chromabase.h"

#include "qdp_map_obj_memory.h"
#include "qdp_map_obj_disk.h"
#include "qdp_map_obj_disk_multiple.h"
#include "handle.h"

#include "meas/inline/hadron/inline_tetraquark_matelem_colorvec_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ferm/key_val_db.h"
#include "util/ferm/key_timeslice_colorvec.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

using std::string;
using std::endl;

namespace Chroma 
{ 
  /*!
   * \ingroup hadron
   *
   * @{
   */
  namespace InlineTetraQuarkMatElemColorVecEnv 
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

      read(paramtop, "mom2_min", param.mom2_min);
      read(paramtop, "mom2_max", param.mom2_max);
      read(paramtop, "mom_list", param.mom_list);
      read(paramtop, "num_vecs", param.num_vecs);
      read(paramtop, "decay_dir", param.decay_dir);
      #ifdef BUILD_SB
      param.use_superb_format=true;
#else
      param.use_superb_format=false;
#endif
      if(paramtop.count("use_superb_format"))read(paramtop,"use_superb_format",param.use_superb_format);
#ifndef BUILD_SB
      if(param.use_superb_format)throw std::runtime_error("SUPERB output requires a BUILD_SB Chroma build");
#endif
    }


    // Writer for input parameters
    void write(XMLWriter& xml, const string& path, const Params::Param_t& param)
    {
      push(xml, path);

      int version = 3;

      write(xml, "version", version);
      write(xml, "mom2_min", param.mom2_min);
      write(xml, "mom2_max", param.mom2_max);
      write(xml, "mom_list", param.mom_list);
      write(xml, "num_vecs", param.num_vecs);
      write(xml, "decay_dir", param.decay_dir);
      write(xml,"use_superb_format",param.use_superb_format);

      pop(xml);
    }

    //! Read named objects 
    void read(XMLReader& xml, const string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "colorvec_file", input.colorvec_file);
      read(inputtop, "tetraquark_op_file", input.tetraquark_op_file);
    }

    //! Write named objects
    void write(XMLWriter& xml, const string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "colorvec_file", input.colorvec_file);
      write(xml, "tetraquark_op_file", input.tetraquark_op_file);

      pop(xml);
    }

    // Writer for input parameters
    void write(XMLWriter& xml, const string& path, const Params& param)
    {
      param.writeXML(xml, path);
    }
  }


  namespace InlineTetraQuarkMatElemColorVecEnv 
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

    const std::string name = "TETRAQUARK_MATELEM_COLORVEC";

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
      param.mom2_min = 0;
      param.mom2_max = 0;
      param.mom_list.resize(0);
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
    //! Key for TetraQuark operator
    struct KeyTetraQuarkElementalOperator_t
    {
      int              t_slice;             /*!< TetraQuark operator time slice */ 
      std::string      color;               /*!< Type of color construction */
      multi1d<int>     mom;                 /*!< D-1 momentum of this operator */
    };

    //! Value for TetraQuark operator
    struct ValTetraQuarkElementalOperator_t
    {
      multi4d<ComplexD>  op;           /*!< Four quarks diquark-antidiquark with momentum projection */
    };


    //----------------------------------------------------------------------------
    //! Used for error output
    StandardOutputStream& operator<<(StandardOutputStream& os, const KeyTetraQuarkElementalOperator_t& d)
    {
      os << "KeyTetraQuarkElementalOperator_t:"
	 << " t_slice= " << d.t_slice
	 << " color= " << d.color
	 << " mom= " << d.mom
	 << std::endl;

      return os;
    }

    
    //----------------------------------------------------------------------------
    // Read a key
    void read(XMLReader& xml, const std::string& path, KeyTetraQuarkElementalOperator_t& param)
    {
      XMLReader paramtop(xml, path);

      read(paramtop, "t_slice", param.t_slice);
      read(paramtop, "color", param.color);
      read(paramtop, "mom", param.mom);
    }

    // Write a key
    void write(XMLWriter& xml, const std::string& path, const KeyTetraQuarkElementalOperator_t& param)
    {
      push(xml, path);

      write(xml, "t_slice", param.t_slice);
      write(xml, "color", param.color);
      write(xml, "mom", param.mom);

      pop(xml);
    }


    //----------------------------------------------------------------------------
    //! KeyTetraQuarkElementalOperator reader
    void read(BinaryReader& bin, KeyTetraQuarkElementalOperator_t& param)
    {
      readDesc(bin, param.color);
      read(bin, param.t_slice);
      read(bin, param.mom);
    }

    //! TetraQuarkElementalOperator write
    void write(BinaryWriter& bin, const KeyTetraQuarkElementalOperator_t& param)
    {
      writeDesc(bin, param.color);
      write(bin, param.t_slice);
      write(bin, param.mom);
    }


    //----------------------------------------------------------------------------
    //! TetraQuarkElementalOperator reader
    void read(BinaryReader& bin, ValTetraQuarkElementalOperator_t& param)
    {
      read(bin, param.op);
    }

    //! TetraQuarkElementalOperator write
    void write(BinaryWriter& bin, const ValTetraQuarkElementalOperator_t& param)
    {
      write(bin, param.op);
    }


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

	push(xml_out, "TetraQuarkMatElemColorVec");
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
			
      //
      // Read in the source along with relevant information.
      // 
      QDPIO::cout << "Snarf the source from a map object disk file" << endl;

      ColorvecMODReader eigen_source;
      eigen_source.setDebug(0);

      try
      {
	// Open
	QDPIO::cout << "Open file= " << params.named_obj.colorvec_file << endl;
	eigen_source.open({params.named_obj.colorvec_file});

	// Snarf the source info. 
	string user_str;
	QDPIO::cout << "Get user data" << endl;
	eigen_source.getUserdata(user_str);
	//	QDPIO::cout << "User data= " << user_str << endl;

	// Write it
	QDPIO::cout << "Write to an xml file" << endl;
	// An XML declaration cannot be embedded inside Source_info.
        if(user_str.find("<?xml")==0)user_str.erase(0,user_str.find("?>")+2);
        XMLBufferWriter xml_buf(user_str);
	write(xml_out, "Source_info", xml_buf);
      }    
      catch (std::bad_cast) {
	QDPIO::cerr << name << ": caught dynamic cast error" << endl;
	QDP_abort(1);
      }
      catch (const string& e) {
	QDPIO::cerr << name << ": error extracting source_header: " << e << endl;
	QDP_abort(1);
      }
      catch( const char* e) {
	QDPIO::cerr << name <<": Caught some char* exception:" << endl;
	QDPIO::cerr << e << endl;
	QDPIO::cerr << "Rethrowing" << endl;
	throw;
      }

      QDPIO::cout << "Source successfully read and parsed" << endl;

      // Sanity check
      if (params.param.num_vecs > eigen_source.size())
      {
	QDPIO::cerr << name << ": number of available eigenvectors is too small\n";
	QDP_abort(1);
      }

      push(xml_out, "TetraQuarkMatElemColorVec");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": TetraQuark color-vector matrix element" << endl;

      proginfo(xml_out);    // Print out basic program info

      // Write out the input
      params.writeXML(xml_out, "Input");

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      //
      // Initialize the slow Tetraier transform phases
      //
      SftMom phases(params.param.mom2_max, false, params.param.decay_dir);
      
      //
      // If a list of momenta has been specified only need phases corresponding to these
      //
      if (params.param.mom_list.size() > 0)
      {
	int num_mom = params.param.mom_list.size();
	int mom_size = params.param.mom_list[0].size();
	multi2d<int> moms(num_mom,mom_size);
	for(int i = 0; i < num_mom; i++)
	{
	  moms[i] = params.param.mom_list[i];
	}
	SftMom temp_phases(moms, params.param.decay_dir);
	phases = temp_phases;
	params.param.mom2_min = 0;
      }

      //
      // DB storage
      //
      BinaryStoreDB< SerialDBKey<KeyTetraQuarkElementalOperator_t>, SerialDBData<ValTetraQuarkElementalOperator_t> > qdp_db;

      // Open the file, and write the meta-data and the binary for this operator
      if (!params.param.use_superb_format && ! qdp_db.fileExists(params.named_obj.tetraquark_op_file))
      {
	XMLBufferWriter file_xml;

	push(file_xml, "DBMetaData");
	write(file_xml, "id", string("tetraQuarkElemOp"));
	write(file_xml, "lattSize", QDP::Layout::lattSize());
	write(file_xml, "decay_dir", params.param.decay_dir);
	proginfo(file_xml);    // Print out basic program info
	write(file_xml, "Params", params.param);
	pop(file_xml);

	std::string file_str(file_xml.str());
	qdp_db.setMaxUserInfoLen(file_str.size());

	qdp_db.open(params.named_obj.tetraquark_op_file, O_RDWR | O_CREAT, 0664);

	qdp_db.insertUserdata(file_str);
      }
      else if (!params.param.use_superb_format)
      {
	qdp_db.open(params.named_obj.tetraquark_op_file, O_RDWR, 0664);
      }


      //
      // TetraQuark operators
      //
      // The creation and annihilation operators are the same without the
      // spin matrices.
      //
      QDPIO::cout << "Building tetraquark operators" << endl;

      push(xml_out, "ElementalOps");


      // Convenience copy - this could be extracted at any point later, but easier to hold
      const int Nt = Layout::lattSize()[params.param.decay_dir];   // time extent

      //
      // The SU(3) irrep outputs
      //
      typedef std::vector< std::pair<std::string, int> > Irreps_t;
      Irreps_t irreps;

      irreps.push_back(std::make_pair("3barx3->1,1", -1));
      irreps.push_back(std::make_pair("6x6bar->1,1", +1));


#ifdef BUILD_SB
      // Tensor convention: i is fastest, then j,k,l,time,momentum,color.
      // Keep the two color constructions explicit in metadata.
      std::vector<std::string> colors;
      for(const auto& ir:irreps)colors.push_back(ir.first);
      std::vector<multi1d<int>> stored_moms;
      for(int m=0;m<phases.numMom();++m)stored_moms.push_back(phases.numToMom(m));
      SB::StorageTensor<7,SB::ComplexD> tensor_store;
      if(params.param.use_superb_format){
        XMLBufferWriter meta;push(meta,"DBMetaData");
        write(meta,"id",std::string("tetraQuarkElemOp"));write(meta,"storage_version",1);
        write(meta,"lattSize",Layout::lattSize());write(meta,"decay_dir",params.param.decay_dir);
        proginfo(meta);write(meta,"tensorOrder",std::string("ijkltmc"));
        write(meta,"colors",colors);write(meta,"moms",stored_moms);write(meta,"Params",params.param);pop(meta);
        tensor_store=SB::StorageTensor<7,SB::ComplexD>(params.named_obj.tetraquark_op_file,SB::broadcast(meta.str()),"ijkltmc",
          SB::Coor<7>{{params.param.num_vecs,params.param.num_vecs,params.param.num_vecs,params.param.num_vecs,Nt,phases.numMom(),int(colors.size())}},
          SB::Sparse,SB::checksum_type::BlockChecksum,SB::SharedFSFile);
      }

#endif
      //
      // Loop over all time slices for the source. This is the same 
      // as the subsets for  phases
      // Build the operator
      //
      for(int t_slice = 0; t_slice < Nt; ++t_slice)
      {
	swiss.reset();
	swiss.start();

	//
	// Read in the eigenvectors on this time-slice
	//
	QDPIO::cout << name << ": read colorvec" << std::endl;
	multi2d<LatticeComplex> vecs(params.param.num_vecs, Nc);

	StopWatch sniss;
	sniss.reset(); sniss.start();

	for(int i=0; i < params.param.num_vecs; ++i)
	{
	  // Key
	  KeyTimeSliceColorVec_t time_key;
	  time_key.t_slice  = t_slice;
	  time_key.colorvec = i;

	  // Get the value
	  LatticeColorVector v;
	  TimeSliceIO<LatticeColorVector> slice(v,t_slice);
	  if (eigen_source.get(time_key, slice) != 0)
	  {
	    QDPIO::cerr << name << ": colorvec time-slice key not found\n";
	    QDP_abort(1);
	  }

	  // Crack open the color indices
	  for(int c=0; c < Nc; ++c)
	    vecs(i,c) = peekColor(v, c);
	}

	sniss.stop(); 
	QDPIO::cout << " Time to read colorvecs: t= " << t_slice 
		    << " time = " << sniss.getTimeInSeconds() << " secs " <<endl;


	
	//
	// Make the diquark. No conjugation yet.
	//
	multi2d< multi2d<LatticeComplex> > diquark(params.param.num_vecs, params.param.num_vecs);

	// Multiply the left and right indices
	for(int i = 0 ; i <  params.param.num_vecs; ++i)
	{
	  for(int j = 0 ; j < params.param.num_vecs; ++j)
	  {
	    diquark(i,j).resize(Nc,Nc);
		
	    for(int a = 0 ; a < Nc; ++a)
	    {
	      for(int b = 0 ; b < Nc; ++b)
	      {
		diquark(i,j)(a,b) = vecs(i,a) * vecs(j,b);
	      } // end for b
	    } // end for a
	  } // end for j
	} // end for i
	

	//
	// Make the two terms in a tetraquark. Contract over the color indices
	//
	multi4d<LatticeComplex> tetra1(params.param.num_vecs,params.param.num_vecs,params.param.num_vecs,params.param.num_vecs);
	multi4d<LatticeComplex> tetra2(params.param.num_vecs,params.param.num_vecs,params.param.num_vecs,params.param.num_vecs);

	tetra1 = zero;
	tetra2 = zero;

	for(int i = 0 ; i <  params.param.num_vecs; ++i)
	{
	  for(int j = 0 ; j < params.param.num_vecs; ++j)
	  {
	    for(int k = 0 ; k <  params.param.num_vecs; ++k)
	    {
	      for(int l = 0 ; l < params.param.num_vecs; ++l)
	      {
		for(int a = 0 ; a < Nc; ++a)
		{
		  for(int b = 0 ; b < Nc; ++b)
		  {
		    tetra1(i,j,k,l) += diquark(i,j)(a,b) * conj(diquark(k,l)(a,b));
		    tetra2(i,j,k,l) += diquark(i,j)(a,b) * conj(diquark(k,l)(b,a));
		  } // end for b
		} // end for a
	      } // end for k
	    } // end for l
	  } // end for j
	} // end for i


	//
	// Loop over the irrep structures
	//
	for(Irreps_t::const_iterator ir = irreps.begin(); ir != irreps.end(); ++ir)
	{
	  // Combine the two terms for this irrep. It's either symmetrization or antisymmetrization
	  Real sym(ir->second);

	  multi4d<LatticeComplex> tetra(params.param.num_vecs,params.param.num_vecs,params.param.num_vecs,params.param.num_vecs);
	  tetra = zero;

	  for(int i = 0 ; i <  params.param.num_vecs; ++i)
	  {
	    for(int j = 0 ; j < params.param.num_vecs; ++j)
	    {
	      for(int k = 0 ; k <  params.param.num_vecs; ++k)
	      {
		for(int l = 0 ; l < params.param.num_vecs; ++l)
		{
		  tetra(i,j,k,l) += tetra1(i,j,k,l) + sym * tetra2(i,j,k,l);
		} // end for l
	      } // end for k
	    } // end for j
	  } // end for i


	  // Big loop over the momentum projection
	  for(int mom_num = 0 ; mom_num < phases.numMom() ; ++mom_num) 
	  {
	    if ( norm2(phases.numToMom(mom_num)) < params.param.mom2_min ) continue;

	    // The keys for this particular elemental operator
	    SerialDBKey<KeyTetraQuarkElementalOperator_t> key;
	    key.key().t_slice = t_slice;
	    key.key().mom     = phases.numToMom(mom_num);
	    key.key().color   = ir->first;

	    // Tetraquark output object
	    SerialDBData<ValTetraQuarkElementalOperator_t> val;
	    val.data().op.resize(params.param.num_vecs,params.param.num_vecs,params.param.num_vecs,params.param.num_vecs);

	    for(int i = 0 ; i <  params.param.num_vecs; ++i)
	    {
	      for(int j = 0 ; j < params.param.num_vecs; ++j)
	      {
		for(int k = 0 ; k <  params.param.num_vecs; ++k)
		{
		  for(int l = 0 ; l < params.param.num_vecs; ++l)
		  {
		    // Do the slow fourier-transform.
		    // Accumulate with the appropriate weight.
		    val.data().op(i,j,k,l) = sum(phases[mom_num] * tetra(i,j,k,l));
		  } // end for l
		} // end for k
	      } // end for i
	    } // end for ff

	    QDPIO::cout << "insert: mom= " << phases.numToMom(mom_num) << endl; 
	    
#ifdef BUILD_SB
if(params.param.use_superb_format){
              const int n=params.param.num_vecs;
              SB::Tensor<4,SB::ComplexD> block("ijkl",SB::Coor<4>{{n,n,n,n}},SB::OnHost,SB::OnMaster);
              if(Layout::primaryNode())for(int l=0;l<n;++l)for(int k=0;k<n;++k)for(int j=0;j<n;++j)for(int i=0;i<n;++i){
                const auto z=val.data().op(i,j,k,l);
                block.data()[i+n*(j+n*(k+n*l))]=SB::ComplexD(toDouble(real(z)),toDouble(imag(z)));
              }
              tensor_store.kvslice_from_size({{'t',t_slice},{'m',mom_num},{'c',int(ir-irreps.begin())}},{{'t',1},{'m',1},{'c',1}}).copyFrom(block);
            }else
#endif
            qdp_db.insert(key, val);

	  } // end for mom_num
	} // end for ir

	swiss.stop();

	QDPIO::cout << "TetraQuark operator: t_slice= " << t_slice
		    << "  time= "
		    << swiss.getTimeInSeconds() 
		    << " secs" << endl;

      } // end for t_slice

      pop(xml_out); // ElementalOps

      // Close the namelist output file XMLDAT
      pop(xml_out);     // TetraQuarkMatElemColorVector

      snoop.stop();
      QDPIO::cout << name << ": total time = " 
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << name << ": ran successfully" << endl;

      END_CODE();
    } // func
  } // namespace InlineTetraQuarkMatElemColorVecEnv

  /*! @} */  // end of group hadron

} // namespace Chroma


