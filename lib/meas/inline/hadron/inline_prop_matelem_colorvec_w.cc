// $Id: inline_prop_matelem_colorvec_w.cc,v 1.5 2008-06-29 03:06:57 edwards Exp $
/*! \file
 * \brief Compute the matrix element of  LatticeColorVector*M^-1*LatticeColorVector
 *
 * Propagator calculation on a colorvector
 */

#include "fermact.h"
#include "meas/inline/hadron/inline_prop_matelem_colorvec_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ferm/eigeninfo.h"
#include "util/ferm/transf.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlinePropMatElemColorVecEnv 
  {
  //! Propagator input
  void read(XMLReader& xml, const string& path, InlinePropMatElemColorVecEnv::Params::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "colorvec_id", input.colorvec_id);
    read(inputtop, "prop_op_file", input.prop_op_file);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlinePropMatElemColorVecEnv::Params::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "colorvec_id", input.colorvec_id);
    write(xml, "prop_op_file", input.prop_op_file);

    pop(xml);
  }


  //! Propagator input
  void read(XMLReader& xml, const string& path, InlinePropMatElemColorVecEnv::Params::Param_t::Contract_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "num_vecs", input.num_vecs);
    read(inputtop, "t_sources", input.t_sources);
    read(inputtop, "decay_dir", input.decay_dir);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlinePropMatElemColorVecEnv::Params::Param_t::Contract_t& input)
  {
    push(xml, path);

    write(xml, "num_vecs", input.num_vecs);
    write(xml, "t_sources", input.t_sources);
    write(xml, "decay_dir", input.decay_dir);

    pop(xml);
  }


  //! Propagator input
  void read(XMLReader& xml, const string& path, InlinePropMatElemColorVecEnv::Params::Param_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "Propagator", input.prop);
    read(inputtop, "Contractions", input.contract);
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlinePropMatElemColorVecEnv::Params::Param_t& input)
  {
    push(xml, path);

    write(xml, "Propagator", input.prop);
    write(xml, "Contractions", input.contract);

    pop(xml);
  }


  //! Propagator input
  void read(XMLReader& xml, const string& path, InlinePropMatElemColorVecEnv::Params& input)
  {
    InlinePropMatElemColorVecEnv::Params tmp(xml, path);
    input = tmp;
  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, const InlinePropMatElemColorVecEnv::Params& input)
  {
    push(xml, path);
    
    write(xml, "Param", input.param);
    write(xml, "NamedObject", input.named_obj);

    pop(xml);
  }
  } // namespace InlinePropMatElemColorVecEnv 


  namespace InlinePropMatElemColorVecEnv 
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
      
    const std::string name = "PROP_MATELEM_COLORVEC";

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
    //! Prop operator
    struct PropElementalOperator_t
    {
      int                t_source;      /*!< Source time slice */
      int                colorvec_src;  /*!< Source colorvector index */
      int                colorvec_snk;  /*!< Sink colorvector index */
      int                spin_src;      /*!< Source spin index */
      int                spin_snk;      /*!< Sink spin index */
      multi1d<ComplexD>  corr;          /*!< Projected propagator */
    };


    //----------------------------------------------------------------------------
    //! PropElementalOperator reader
    void read(BinaryReader& bin, PropElementalOperator_t& param)
    {
      read(bin, param.t_source);
      read(bin, param.colorvec_src);
      read(bin, param.colorvec_snk);
      read(bin, param.spin_src);
      read(bin, param.spin_snk);
      read(bin, param.corr);
    }

    //! PropElementalOperator write
    void write(BinaryWriter& bin, const PropElementalOperator_t& param)
    {
      write(bin, param.t_source);
      write(bin, param.colorvec_src);
      write(bin, param.colorvec_snk);
      write(bin, param.spin_src);
      write(bin, param.spin_snk);
      write(bin, param.corr);
    }

    //! PropElementalOperator reader
    void read(XMLReader& xml, const std::string& path, PropElementalOperator_t& param)
    {
      XMLReader paramtop(xml, path);
    
      read(paramtop, "t_source", param.t_source);
      read(paramtop, "colorvec_src", param.colorvec_src);
      read(paramtop, "colorvec_snk", param.colorvec_snk);
      read(paramtop, "spin_src", param.spin_src);
      read(paramtop, "spin_snk", param.spin_snk);
      read(paramtop, "corr", param.corr);
    }

    //! PropElementalOperator writer
    void write(XMLWriter& xml, const std::string& path, const PropElementalOperator_t& param)
    {
      push(xml, path);

      write(xml, "t_source", param.t_source);
      write(xml, "colorvec_src", param.colorvec_src);
      write(xml, "colorvec_snk", param.colorvec_snk);
      write(xml, "spin_src", param.spin_src);
      write(xml, "spin_snk", param.spin_snk);
      write(xml, "corr", param.corr);

      pop(xml);
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
	QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << endl;
	QDP_abort(1);
      }
    }



    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "PropMatElemColorVec");
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
	QDPIO::cerr << name << ": caught dynamic cast error" << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << name << ": map call failed: " << e << endl;
	QDP_abort(1);
      }

      push(xml_out, "PropMatElemColorVec");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": propagator calculation" << endl;

      proginfo(xml_out);    // Print out basic program info

      // Write out the input
      write(xml_out, "Input", params);

      // Write out the config header
      write(xml_out, "Config_info", gauge_xml);

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      // Calculate some gauge invariant observables just for info.
      MesPlq(xml_out, "Observables", u);

      //
      // Read in the source along with relevant information.
      // 
      XMLReader source_file_xml, source_record_xml;

      QDPIO::cout << "Snarf the source from a named buffer" << endl;
      try
      {
	TheNamedObjMap::Instance().getData< EigenInfo<LatticeColorVector> >(params.named_obj.colorvec_id);

	// Snarf the source info. This is will throw if the colorvec_id is not there
	TheNamedObjMap::Instance().get(params.named_obj.colorvec_id).getFileXML(source_file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.colorvec_id).getRecordXML(source_record_xml);

	// Write out the source header
	write(xml_out, "Source_file_info", source_file_xml);
	write(xml_out, "Source_record_info", source_record_xml);
      }    
      catch (std::bad_cast)
      {
	QDPIO::cerr << name << ": caught dynamic cast error" << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << name << ": error extracting source_header: " << e << endl;
	QDP_abort(1);
      }
      const EigenInfo<LatticeColorVector>& eigen_source = 
	TheNamedObjMap::Instance().getData< EigenInfo<LatticeColorVector> >(params.named_obj.colorvec_id);

      QDPIO::cout << "Source successfully read and parsed" << endl;

      // Sanity check - write out the norm2 of the source in the Nd-1 direction
      // Use this for any possible verification
      {
	// Initialize the slow Fourier transform phases
	SftMom phases(0, true, Nd-1);

	multi1d< multi1d<Double> > source_corrs(eigen_source.getEvalues().size());
	for(int m=0; m < source_corrs.size(); ++m)
	{
	  source_corrs[m] = sumMulti(localNorm2(eigen_source.getEvectors()[m]), phases.getSet());
	}

	push(xml_out, "Source_correlators");
	write(xml_out, "source_corrs", source_corrs);
	pop(xml_out);
      }

      // Another sanity check
      if (params.param.contract.num_vecs > eigen_source.getEvalues().size())
      {
	QDPIO::cerr << __func__ << ": num_vecs= " << params.param.contract.num_vecs
		    << " is greater than the number of available colorvectors= "
		    << eigen_source.getEvalues().size() << endl;
	QDP_abort(1);
      }


      // Total number of iterations
      int ncg_had = 0;

      // Number of props written
      int total_num_elem = 0;

      // Binary output
      BinaryBufferWriter src_record_bin;

      //
      // Try the factories
      //
      try
      {
	StopWatch swatch;
	swatch.reset();
	QDPIO::cout << "Try the various factories" << endl;

	// Typedefs to save typing
	typedef LatticeFermion               T;
	typedef multi1d<LatticeColorMatrix>  P;
	typedef multi1d<LatticeColorMatrix>  Q;

	//
	// Initialize fermion action
	//
	std::istringstream  xml_s(params.param.prop.fermact.xml);
	XMLReader  fermacttop(xml_s);
	QDPIO::cout << "FermAct = " << params.param.prop.fermact.id << endl;

	// Generic Wilson-Type stuff
	Handle< FermionAction<T,P,Q> >
	  S_f(TheFermionActionFactory::Instance().createObject(params.param.prop.fermact.id,
							       fermacttop,
							       params.param.prop.fermact.path));

	Handle< FermState<T,P,Q> > state(S_f->createState(u));

	Handle< SystemSolver<LatticeFermion> > PP = S_f->qprop(state,
							       params.param.prop.invParam);
      
	QDPIO::cout << "Suitable factory found: compute all the quark props" << endl;
	swatch.start();

	//
	// Loop over the source color and spin, creating the source
	// and calling the relevant propagator routines.
	//
	const int num_vecs            = params.param.contract.num_vecs;
	const int decay_dir           = params.param.contract.decay_dir;
	const multi1d<int>& t_sources = params.param.contract.t_sources;

	// Initialize the slow Fourier transform phases
	SftMom phases(0, true, decay_dir);

	// Binary output
	int num_elem_to_write = num_vecs * num_vecs * t_sources.size() * Ns * Ns;

	write(src_record_bin, decay_dir);
	write(src_record_bin, num_vecs);
	write(src_record_bin, t_sources);
	write(src_record_bin, num_elem_to_write);

	push(xml_out, "ColorVecMatElems");

	// Loop over each operator 
	for(int tt=0; tt < t_sources.size(); ++tt)
	{
	  int t_source = t_sources[tt];
	  QDPIO::cout << "t_source = " << t_source << endl; 

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

	      // Insert a ColorVector into spin index spin_source
	      // This only overwrites sections, so need to initialize first
	      LatticeFermion chi = zero;
	      CvToFerm(vec_srce, chi, spin_source);

	      LatticeFermion quark_soln = zero;

	      // Do the propagator inversion
	      SystemSolverResults_t res = (*PP)(quark_soln, chi);
	      ncg_had = res.n_count;

	      for(int colorvec_sink=0; colorvec_sink < num_vecs; ++colorvec_sink)
	      {
		const LatticeColorVector& vec_sink = eigen_source.getEvectors()[colorvec_sink];

		for(int spin_sink=0; spin_sink < Ns; ++spin_sink)
		{
		  PropElementalOperator_t op;
		  op.t_source     = t_source;
		  op.colorvec_src = colorvec_source;
		  op.colorvec_snk = colorvec_sink;
		  op.spin_src     = spin_source;
		  op.spin_snk     = spin_sink;
		  op.corr         = sumMulti(localInnerProduct(vec_sink, peekSpin(quark_soln, spin_sink)), 
					     phases.getSet());
		  
		  write(src_record_bin, op);    // binary output
		  write(xml_out, "elem", op);   // xml output (debugging)
		  ++total_num_elem;
		} // for spin_sink
	      } // for colorvec_sink
	    } // for spin_source
	  } // for colorvec_source
	} // for t_source

	pop(xml_out);

	// Sanity check
	if (total_num_elem != num_elem_to_write)
	{
	  QDPIO::cerr << name << ": inconsistent number of elemental ops written = "
		      << total_num_elem
		    << endl;
	  QDP_abort(1);
	}

	swatch.stop();
	QDPIO::cout << "Propagators computed: time= " 
		    << swatch.getTimeInSeconds() 
		    << " secs" << endl;
      }
      catch (const std::string& e) 
      {
	QDPIO::cout << name << ": caught exception around qprop: " << e << endl;
	QDP_abort(1);
      }

      push(xml_out,"Relaxation_Iterations");
      write(xml_out, "ncg_had", ncg_had);
      pop(xml_out);

      pop(xml_out);  // prop_matelem_colorvec

      // Write the meta-data and the binary for this operator
      {
	XMLBufferWriter src_record_xml, file_xml;

	push(file_xml, "PropElementalOperators");
	write(file_xml, "Params", params.param);
	write(file_xml, "Config_info", gauge_xml);
	pop(file_xml);

	QDPFileWriter qdp_file(file_xml, params.named_obj.prop_op_file,
			       QDPIO_SINGLEFILE, QDPIO_SERIAL, QDPIO_OPEN);

	push(src_record_xml, "PropElementalOperator");
	write(src_record_xml, "decay_dir", params.param.contract.decay_dir);
	write(src_record_xml, "num_vecs", params.param.contract.num_vecs);
	write(src_record_xml, "t_sources", params.param.contract.t_sources);
	write(src_record_xml, "total_num_elem", total_num_elem);
	pop(src_record_xml);

	write(qdp_file, src_record_xml, src_record_bin);
      }

      snoop.stop();
      QDPIO::cout << name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << name << ": ran successfully" << endl;

      END_CODE();
    } 

  }

} // namespace Chroma
