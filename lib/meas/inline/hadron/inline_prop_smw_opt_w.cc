/*! \file
 * \brief Inline measurement that construct unsmeared hadron nodes using distillation
 */

#include "meas/inline/hadron/inline_prop_smw_opt_w.h"

#include "qdp_map_obj_memory.h"
#include "qdp_map_obj_disk.h"
#include "qdp_map_obj_disk_multiple.h"
#include "qdp_disk_map_slice.h"
#include "handle.h"

#include "meas/glue/mesplq.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"
#include "util/ferm/key_timeslice_colorvec.h"
#include "util/ferm/key_val_db.h"
#include "util/ferm/key_prop_matelem.h"
#include "util/info/proginfo.h"
#include "meas/smear/displace.h"
#include "util/ft/sftmom.h"
#include "util/ft/time_slice_set.h"
#include "io/xml_group_reader.h"
#include "meas/inline/make_xml_file.h"

#include "util/ferm/key_val_db.h"
#include "util/ferm/transf.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

#include <tensor/tensor.h>
#include <tensor/linalg.h>
#include <tensor/vector.h>
#include <complex>

namespace Chroma 
{ 
  /*!
   * \ingroup hadron
   *
   * @{
   */
  namespace InlinePropSMWOptEnv
  { 
    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "colorvec_files", input.colorvec_files);
      read(inputtop, "prop_op_file", input.prop_op_file);
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "colorvec_files", input.colorvec_files);
      write(xml, "prop_op_file", input.prop_op_file);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::Param_t::Contract_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "num_vecs", input.num_vecs);
      read(inputtop, "decay_dir", input.decay_dir);
      read(inputtop, "mass_label", input.mass_label);
      read(inputtop, "num_tries", input.num_tries);
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t::Contract_t& input)
    {
      push(xml, path);

      write(xml, "num_vecs", input.num_vecs);
      write(xml, "decay_dir", input.decay_dir);
      write(xml, "mass_label", input.mass_label);
      write(xml, "num_tries", input.num_tries);

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const std::string& path, Params::Param_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "Propagator", input.prop);
      read(inputtop, "Contractions", input.contract);
    }

    //! Propagator output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t& input)
    {
      push(xml, path);

      write(xml, "Propagator", input.prop);
      write(xml, "Contractions", input.contract);

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
  } // end namespace



  //-------------------------------------------------------------------------------
  //-------------------------------------------------------------------------------
  namespace InlinePropSMWOptEnv
  { 
    //----------------------------------------------------------------------------
    // Convenience type
    typedef QDP::MapObjectDisk<KeyTimeSliceColorVec_t, TimeSliceIO<LatticeColorVectorF> > MOD_t;

    // Convenience type
    typedef QDP::MapObjectDiskMultiple<KeyTimeSliceColorVec_t, TimeSliceIO<LatticeColorVectorF> > MODS_t;

    // Convenience type
    typedef QDP::MapObjectMemory<KeyTimeSliceColorVec_t, SubLatticeColorVectorF> SUB_MOD_t;

  } // end namespace


  //-------------------------------------------------------------------------------
  namespace InlinePropSMWOptEnv
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

    const std::string name = "PROP_SMW_OPT";

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

      StandardOutputStream& operator<<(StandardOutputStream& os, const std::vector<int>& d)
      {
	if (d.size() > 0)
	{
	  os << d[0];

	  for(int i=1; i < d.size(); ++i)
	    os << " " << d[i];
	}

	return os;
      }

      //----------------------------------------------------------------------------
      // Lazy conversion
      ComplexD quickConvert(const std::complex<double>& src)
      {
	return cmplx(Real(src.real()),Real(src.imag()));
      }


      // Lazy conversion
      std::complex<double> quickConvert(const ComplexD& src)
      {
	return std::complex<double>(toDouble(real(src)),toDouble(imag(src)));
      }

      //----------------------------------------------------------------------------
      //! Convenience
      inline int packSpinDist(int s, int d) {return s + Ns*d;}
    }


    //-------------------------------------------------------------------------------
    //-------------------------------------------------------------------------------
    //! Invert off of each source, and do all the checking
    LatticeFermion doInversion(const SystemSolver<LatticeFermion>& PP, const LatticeColorVector& vec_srce, int spin_source, int num_tries)
    {
      //
      // Loop over each spin source and invert. 
      //
      LatticeFermion chi = zero;
      CvToFerm(vec_srce, chi, spin_source);

      LatticeFermion quark_soln = zero;
      int ncg_had;

      // Do the propagator inversion
      // Check if bad things are happening
      bool badP = true;
      for(int nn = 1; nn <= num_tries; ++nn)
      {	
	// Reset
	quark_soln = zero;
	badP = false;
	      
#if 0
	// Solve for the solution std::vector
	SystemSolverResults_t res = PP(quark_soln, chi);
#else
	SystemSolverResults_t res;
        res.n_count = 1;	
        res.resid = 1.0e-7;
#endif
	ncg_had = res.n_count;

	// Some sanity checks
	if (toDouble(res.resid) > 1.0e-3)
	{
	  QDPIO::cerr << name << ": have a resid > 1.0e-3. That is unacceptable" << std::endl;
	  QDP_abort(1);
	}

	// Check for finite values - neither NaN nor Inf
	if (isfinite(quark_soln))
	{
	  // Okay
	  break;
	}
	else
	{
	  QDPIO::cerr << name << ": WARNING - found something not finite, may retry\n";
	  badP = true;
	}
      }

      // Sanity check
      if (badP)
      {
	QDPIO::cerr << name << ": this is bad - did not get a finite solution vector after num_tries= " << num_tries << std::endl;
	QDP_abort(1);
      }

      return quark_soln;
    }
    

    //-------------------------------------------------------------------------------
    //! Get full propagator
    void getSoln(LatticeColorVectorSpinMatrix& soln, const SystemSolver<LatticeFermion>& PP, const LatticeColorVector& vec_srce)
    {
      // Get the source vector
      StopWatch swatch;
      swatch.reset();
      swatch.start();

      // Loop over each spin source
      for(int spin_ind=0; spin_ind < Ns; ++spin_ind)
      {	
	QDPIO::cout << __func__ << ": do spin= " << spin_ind << std::endl; 
	    
	StopWatch snarss1;
	snarss1.reset();
	snarss1.start();

	//
	// Loop over each spin source and invert. 
	// Use the same color vector source. No spin dilution will be used.
	//
	LatticeFermion tmp = doInversion(PP, vec_srce, spin_ind, 1);

	// shove a fermion into a colorvec-spinmatrix 
	FermToProp(tmp, soln, spin_ind);

	snarss1.stop();
	QDPIO::cout << __func__ << ": time to compute prop for spin_ind= " << spin_ind
		    << "  time = " << snarss1.getTimeInSeconds() 
		    << " secs" << std::endl;
      } // for spin_src

      swatch.stop(); 
      QDPIO::cout << __func__ << ": FINISHED: total time = " << swatch.getTimeInSeconds() << " secs" << std::endl;
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


    //-------------------------------------------------------------------------------
    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	std::string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "HadronNode");
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


    //-------------------------------------------------------------------------------
    // Function call
    void 
    InlineMeas::func(unsigned long update_no,
		     XMLWriter& xml_out) 
    {
      START_CODE();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      QDPIO::cout << name << ": construct unsmeared hadron nodes via distillation" << std::endl;
      
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
	QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
	QDP_abort(1);
      }
      catch (const std::string& e) 
      {
	QDPIO::cerr << name << ": std::map call failed: " << e << std::endl;
	QDP_abort(1);
      }

      push(xml_out, "PropSMW");
      write(xml_out, "update_no", update_no);

      // Write out the input
      write(xml_out, "Input", params);

      // Write out the config header
      write(xml_out, "Config_info", gauge_xml);

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      // Calculate some gauge invariant observables just for info.
      MesPlq(xml_out, "Observables", u);

      // Will use TimeSliceSet-s a lot
      const int decay_dir = params.param.contract.decay_dir;
      const int Lt        = Layout::lattSize()[decay_dir];
      const int sLt       = Layout::subgridLattSize()[decay_dir];
      const int num_vecs  = params.param.contract.num_vecs;

      // A sanity check
      if (decay_dir != Nd-1)
      {
	QDPIO::cerr << name << ": TimeSliceIO only supports decay_dir= " << Nd-1 << "\n";
	QDP_abort(1);
      }

      // Reset
      if (params.param.contract.num_tries <= 0)
      {
	params.param.contract.num_tries = 1;
      }


      //
      // Read in the source along with relevant information.
      // 
      QDPIO::cout << "Snarf the source from a std::map object disk file" << std::endl;

      MODS_t eigen_source;
      eigen_source.setDebug(0);

      std::string eigen_meta_data;   // holds the eigenvalues

      try
      {
	// Open
	QDPIO::cout << "Open file= " << params.named_obj.colorvec_files[0] << std::endl;
	eigen_source.open(params.named_obj.colorvec_files);

	// Snarf the source info. 
	QDPIO::cout << "Get user data" << std::endl;
	eigen_source.getUserdata(eigen_meta_data);
	//	QDPIO::cout << "User data= " << eigen_meta_data << std::endl;

	// Write it
	//	QDPIO::cout << "Write to an xml file" << std::endl;
	//	XMLBufferWriter xml_buf(eigen_meta_data);
	//	write(xml_out, "Source_info", xml_buf);
      }    
      catch (std::bad_cast) {
	QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
	QDP_abort(1);
      }
      catch (const std::string& e) {
	QDPIO::cerr << name << ": error extracting source_header: " << e << std::endl;
	QDP_abort(1);
      }
      catch( const char* e) {
	QDPIO::cerr << name <<": Caught some char* exception:" << std::endl;
	QDPIO::cerr << e << std::endl;
	QDPIO::cerr << "Rethrowing" << std::endl;
	throw;
      }

      QDPIO::cout << "Source successfully read and parsed" << std::endl;

#if 0
      // Sanity check
      if (num_vecs > eigen_source.size())
      {
	QDPIO::cerr << name << ": number of available eigenvectors is too small\n";
	QDP_abort(1);
      }
      QDPIO::cout << "Number of vecs available is large enough" << std::endl;
#endif

      //
      // Initialize the slow Fourier transform phases
      //
      QDPIO::cout << name << ": initialize phases" << std::endl;
      SftMom phases(0, true, decay_dir);

    
      //
      // DB storage
      //
      BinaryStoreDB< SerialDBKey<KeyPropElementalOperator_t>, SerialDBData<ValPropElementalOperator_t> > qdp_db;

      // Open the file, and write the meta-data and the binary for this operator
      if (! qdp_db.fileExists(params.named_obj.prop_op_file))
      {
	XMLBufferWriter file_xml;

	push(file_xml, "DBMetaData");
	write(file_xml, "id", std::string("propElemOp"));
	write(file_xml, "lattSize", QDP::Layout::lattSize());
	write(file_xml, "decay_dir", params.param.contract.decay_dir);
	proginfo(file_xml);    // Print out basic program info
	pop(file_xml);

	std::string file_str(file_xml.str());
	qdp_db.setMaxUserInfoLen(file_str.size());

	qdp_db.open(params.named_obj.prop_op_file, O_RDWR | O_CREAT, 0664);

	qdp_db.insertUserdata(file_str);
      }
      else
      {
	qdp_db.open(params.named_obj.prop_op_file, O_RDWR, 0664);
      }

      QDPIO::cout << "Finished opening distillation file" << std::endl;


      //
      // Hack up the gauge field
      //
      multi1d<LatticeColorMatrix> u_hack = u;
      u_hack[decay_dir] = zero;



      //
      // Try the factories
      //
      try
      {
	StopWatch swatch;
	swatch.reset();
	QDPIO::cout << "Try the various factories" << std::endl;

	// Typedefs to save typing
	typedef LatticeFermion               T;
	typedef multi1d<LatticeColorMatrix>  P;
	typedef multi1d<LatticeColorMatrix>  Q;

	//
	// Initialize fermion action
	//
	std::istringstream  xml_s(params.param.prop.fermact.xml);
	XMLReader  fermacttop(xml_s);
	QDPIO::cout << "FermAct = " << params.param.prop.fermact.id << std::endl;

	// Generic Wilson-Type stuff
	Handle< FermionAction<T,P,Q> >
	  S_f(TheFermionActionFactory::Instance().createObject(params.param.prop.fermact.id,
							       fermacttop,
							       params.param.prop.fermact.path));

	Handle< FermState<T,P,Q> > state(S_f->createState(u_hack));

	Handle< SystemSolver<LatticeFermion> > PP = S_f->qprop(state,
							       params.param.prop.invParam);
      
	QDPIO::cout << "Suitable factory found: compute all the quark props" << std::endl;
	swatch.start();

	//
	// The perambulator tensor for all time-slices
	// 
	std::vector<tensor::CTensor> spatial_peram(Lt);
	for(int t = 0; t < Lt; ++t)
	  spatial_peram[t] = tensor::CTensor(num_vecs*Ns, num_vecs*Ns);


	// Load up all time-slices into one vector
	// Because the vectors have been split out in individual time-slices, must reload
	// NOTE: we do have the lime version of the eigs which have all time-slices packed in
	std::vector<LatticeColorVector> eig_vecs(num_vecs);

	for(int colorvec_src=0; colorvec_src < num_vecs; ++colorvec_src)
	{
	  LatticeColorVectorF vec_srce = zero;

	  for(int t_slice=0; t_slice < Lt; ++t_slice)
	  {
	    QDPIO::cout << "load t_slice = " << t_slice << "  colorvec_src= " << colorvec_src << std::endl; 

	    KeyTimeSliceColorVec_t src_key(t_slice, colorvec_src);
	    TimeSliceIO<LatticeColorVectorF> time_slice_io(vec_srce, t_slice);
	    eigen_source.get(src_key, time_slice_io);
	  }

	  eig_vecs[colorvec_src] = vec_srce;
	}
	
	QDPIO::cout << "Finished loading up eigs" << std::endl; 


	//
	// Make a t-slice perambulator
	//
	std::vector<tensor::CTensor> time_plus_peram(Lt);
	std::vector<tensor::CTensor> time_minus_peram(Lt);
	for(int t = 0; t < Lt; ++t)
	{
	  time_plus_peram[t]  = tensor::CTensor(num_vecs*Ns, num_vecs*Ns);
	  time_minus_peram[t] = tensor::CTensor(num_vecs*Ns, num_vecs*Ns);
	}

	if (1)
	{
	  SpinMatrix one(1.0);
	  SpinMatrix g4p = 0.5 * (one + Gamma(8)*one);
	  SpinMatrix g4m = 0.5 * (one - Gamma(8)*one);

	  for(int colorvec_src=0; colorvec_src < num_vecs; ++colorvec_src)
	  {
	    QDPIO::cout << "colorvec_src = " << colorvec_src << std::endl; 

	    // The temporal perambulator part
	    for(int colorvec_snk=0; colorvec_snk < num_vecs; ++colorvec_snk)
	    {
	      multi1d<ComplexD> fredf = sumMulti(localInnerProduct(eig_vecs[colorvec_snk],
								  u[Nd-1]*shift(eig_vecs[colorvec_src], FORWARD, Nd-1)),
						phases.getSet());

	      multi1d<ComplexD> fredb = sumMulti(localInnerProduct(eig_vecs[colorvec_snk],
								   shift(adj(u[Nd-1])*eig_vecs[colorvec_src], BACKWARD, Nd-1)),
						phases.getSet());

	      for(int spin_snk=0; spin_snk < Ns; ++spin_snk)
	      {
		for(int spin_src=0; spin_src < Ns; ++spin_src)
		{
		  for(int t=0; t < Lt; ++t)
		  {
		    time_plus_peram[t].at(packSpinDist(spin_snk,colorvec_snk), packSpinDist(spin_src,colorvec_src))  = quickConvert(-fredf[t] * peekSpin(g4m, spin_snk, spin_src));
		    time_minus_peram[t].at(packSpinDist(spin_snk,colorvec_snk), packSpinDist(spin_src,colorvec_src)) = quickConvert( fredb[t] * peekSpin(g4p, spin_snk, spin_src));
		  } // for t
		} // for spin_src
	      } // for spin_snk
	    } // for colorvec_snk
	  } // for colorvec_src
	  
	  QDPIO::cout << "Finished constructing temporal perams" << std::endl;
	}


	//
	// Loop over each distillation source and generate the full propagator
	//
	for(int colorvec_src=0; colorvec_src < num_vecs; ++colorvec_src)
	{
	  QDPIO::cout << "colorvec_src = " << colorvec_src << std::endl; 

	  // Generate the propagator
	  LatticeColorVectorSpinMatrix quark_soln;
	  getSoln(quark_soln, *PP, eig_vecs[colorvec_src]);

	  //
	  // The spatial perambulator part
	  //
	  // LOTS of inefficiencies here. Should pack spin into the inner product. Of course, there's the fancy tensor stuff.
	  // Another day... 
	  //
	  for(int colorvec_snk=0; colorvec_snk < num_vecs; ++colorvec_snk)
	  {
	    for(int spin_snk=0; spin_snk < Ns; ++spin_snk)
	    {
	      for(int spin_src=0; spin_src < Ns; ++spin_src)
	      {
		multi1d<ComplexD> fred = sumMulti(localInnerProduct(eig_vecs[colorvec_snk],
								    peekSpin(quark_soln, spin_snk, spin_src)),
						  phases.getSet());
		  
		for(int t=0; t < Lt; ++t)
		{
		  spatial_peram[t].at(packSpinDist(spin_snk,colorvec_snk), packSpinDist(spin_src,colorvec_src)) = quickConvert(fred[t]);
		} // for t
	      } // for spin_src
	    } // for spin_snk
	  } // for colorvec_snk
	} // for colorvec_src

	QDPIO::cout << "Finished constructing spatial perams" << std::endl;


	//
	// Do some fun and wacky stuff. Impress your friends. 
	//
	if (1)
	{
	  QDPIO::cout << "It is play time!" << std::endl;

	  // See tensor  http://juanjosegarciaripoll.github.io/tensor
	  // Fun stuff in linalg.h
	  // This is an identity matrix
	  tensor::CTensor one = tensor::CTensor::eye(num_vecs*Ns);

	  for(int t=0; t < Lt; ++t)
	  {
	    QDPIO::cout << "fun t= " << t << std::endl;
	    // Totally rad
	    // Extract eigenvectors. NOTE: not checking if this makes sense
	    tensor::CTensor R;
	    tensor::RTensor s = linalg::eig_sym(spatial_peram[t], &R);
	    tensor::RTensor dS = diag(s);
	    
	    // Just goofing. I think this is foo = inv(peram)
	    // Note: "solve" says my matrix is singular, so cannot invert. Need SVD.
	    // This seems to fail with unable to fold() tensors
//	    tensor::CTensor foo = linalg::solve_with_svd(spatial_peram[t], one, 1.0e-3);
	  }

	  QDPIO::cout << "Play time is over" << std::endl;
	}

	//
	// Write out something. This will change...
	//
	for(int t=0; t < Lt; ++t)
	{
	  // Repack the tensor into the format used for perambulators
	  for(int spin_snk=0; spin_snk < Ns; ++spin_snk)
	  {
	    for(int spin_src=0; spin_src < Ns; ++spin_src)
	    {
	      // Key stuff
	      KeyPropElementalOperator_t key;

	      key.t_slice      = t;
	      key.t_source     = t;
	      key.spin_src     = spin_src;
	      key.spin_snk     = spin_snk;
	      key.mass_label   = params.param.contract.mass_label;

	      // Have to unpack tensor to make the spin keys explict
	      ValPropElementalOperator_t val;
	      val.mat.resize(num_vecs,num_vecs);

	      for(int colorvec_snk=0; colorvec_snk < num_vecs; ++colorvec_snk)
	      {
		for(int colorvec_src=0; colorvec_src < num_vecs; ++colorvec_src)
		{
		  val.mat(colorvec_snk,colorvec_src) = quickConvert(spatial_peram[t].at(packSpinDist(spin_snk,colorvec_snk), packSpinDist(spin_src,colorvec_src)));
		}
	      }

	      QDPIO::cout << "insert: t_slice= " << t << " spin_src= " << spin_src << " spin_snk= " << spin_snk << std::endl; 
	      for(int t=0; t < phases.numSubsets(); ++t)
	      {
		qdp_db.insert(key, val);
	      }
	    } // for spin_sink
	  } // for spin_source
	} // for t

	// Finished writing
	swatch.stop();
	QDPIO::cout << "Propagators computed: time= " 
		    << swatch.getTimeInSeconds() 
		    << " secs" << std::endl;
      }
      catch (const std::string& e) 
      {
	QDPIO::cout << name << ": caught exception around qprop: " << e << std::endl;
	QDP_abort(1);
      }

      // Close db
      qdp_db.close();

      // Close the xml output file
      pop(xml_out);     // UnsmearedHadronNode

      snoop.stop();
      QDPIO::cout << name << ": total time = " 
		  << snoop.getTimeInSeconds() 
		  << " secs" << std::endl;

      QDPIO::cout << name << ": ran successfully" << std::endl;

      END_CODE();
    } // func
  } // namespace InlineUnsmearedHadronNodeDistillationEnv

  /*! @} */  // end of group hadron

} // namespace Chroma

