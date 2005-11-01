// $Id: inline_make_source_w.cc,v 2.3 2005-11-01 04:10:01 edwards Exp $
/*! \file
 * \brief Inline construction of make_source
 *
 * Construct source for propagator calculations
 */

#include "meas/inline/hadron/inline_make_source_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/smear/ape_smear.h"
#include "meas/smear/gaus_smear.h"
#include "meas/smear/displacement.h"
#include "meas/smear/laplacian.h"
#include "meas/sources/srcfil.h"
#include "meas/sources/walfil_w.h"
#include "meas/sources/p_src_w.h"
#include "meas/sources/d_src_w.h"
#include "util/ferm/transf.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlineMakeSourceEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineMakeSource(InlineMakeSourceParams(xml_in, path));
    }

    const std::string name = "MAKE_SOURCE";
    const bool registered = TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
  };


  //! MakeSource input
  void read(XMLReader& xml, const string& path, InlineMakeSourceParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "source_id", input.source_id);
  }

  //! MakeSource output
  void write(XMLWriter& xml, const string& path, const InlineMakeSourceParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "source_id", input.source_id);

    pop(xml);
  }


  // Param stuff
  InlineMakeSourceParams::InlineMakeSourceParams() { frequency = 0; }

  InlineMakeSourceParams::InlineMakeSourceParams(XMLReader& xml_in, const std::string& path) 
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
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  void
  InlineMakeSourceParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    // Parameters for source construction
    Chroma::write(xml_out, "Param", param);

    // Write out the buffer ids
    Chroma::write(xml_out, "NamedObject", named_obj);

    pop(xml_out);
  }


  void 
  InlineMakeSource::operator()(const multi1d<LatticeColorMatrix>& u,
			       XMLBufferWriter& gauge_xml,
			       unsigned long update_no,
			       XMLWriter& xml_out) 
  {
    START_CODE();

    StopWatch snoop;
    snoop.reset();
    snoop.start();

    // Record the initial state of the RNG (needed for noisy sources)
    QDP::Seed ran_seed ;
    QDP::RNG::savern(ran_seed) ;

    push(xml_out, "make_source");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlineMakeSourceEnv::name << ": propagator source constructor" << endl;

    switch(params.param.source_type)
    {
    case SRC_TYPE_POINT_SOURCE:
      QDPIO::cout << "Point source" << endl;
      break;
    case SRC_TYPE_WALL_SOURCE:
      QDPIO::cout << "Wall source" << endl;
      break;
    case SRC_TYPE_SHELL_SOURCE:
      if (params.param.sourceSmearParam.wvf_kind != WVF_KIND_GAUGE_INV_GAUSSIAN)
      {
	QDPIO::cout << "Unsupported source smearing type" << endl;
	QDP_abort(1);
      }
      QDPIO::cout << "Smeared source wvf_param= " 
		  << params.param.sourceSmearParam.wvf_param 
		  << ": wvfIntPar= "
		  << params.param.sourceSmearParam.wvfIntPar << endl
		  << "Power of Laplacian operator= " << params.param.laplace_power << endl
		  << "Displacement length= " << params.param.disp_length
		  <<": Displacement direction= " << params.param.disp_dir << endl;
      break;
    case SRC_TYPE_RAND_Z2_WALL_SOURCE:
      write(xml_out, "RNG", ran_seed) ;
      QDPIO::cout << "Random complex Z(2) wall source" << endl;
      break;
    case SRC_TYPE_RAND_U1_WALL_SOURCE:
      write(xml_out, "RNG", ran_seed) ;
      QDPIO::cout << "Random U(1) wall source" << endl;
      break;
    default:
      QDPIO::cout << "Unknown source_type" << endl;
      QDP_abort(1);
    }

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.write(xml_out, "Input");

    // Write out the config header
    write(xml_out, "Config_info", gauge_xml);

    // Calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);

    // Smear the gauge field if needed
    multi1d<LatticeColorMatrix> u_smr(Nd);
    u_smr = u;

    if (params.param.source_type == SRC_TYPE_SHELL_SOURCE && 
	params.param.link_smear_num > 0)
    {
      QDPIO::cout << "Smear gauge field" << endl;

      int BlkMax = 100;	// Maximum number of blocking/smearing iterations
      Real BlkAccu = 1.0e-5;	// Blocking/smearing accuracy

      for(int i=0; i < params.param.link_smear_num; ++i)
      {
	multi1d<LatticeColorMatrix> u_tmp(Nd);

	for(int mu = 0; mu < Nd; ++mu)
	  if ( mu != params.param.j_decay )
	    APE_Smear(u_smr, u_tmp[mu], mu, 0, 
		      params.param.link_smear_fact, BlkAccu, BlkMax, 
		      params.param.j_decay);
	  else
	    u_tmp[mu] = u_smr[mu];

	u_smr = u_tmp;
      }

      QDPIO::cout << "Gauge field APE-smeared!" << endl;
    }


    //
    // Create a propagator object
    //
    LatticePropagator quark_source;

    //
    // Loop over the source color and spin, creating the source
    // and calling the relevant propagator routines. The QDP
    // terminology is that a propagator is a matrix in color
    // and spin space
    //
    // For this calculation, a smeared source is used. A point
    // source is first constructed and then smeared. If a user
    // only wanted a point source, then remove the smearing stuff
    //
    switch (params.param.source_type)
    {
      //
      // Gauge inv. point or shell sources within some S_WAVE, P_WAVE, state etc.
      //
    case SRC_TYPE_POINT_SOURCE:
    case SRC_TYPE_SHELL_SOURCE:
    case SRC_TYPE_RAND_Z2_WALL_SOURCE:
    case SRC_TYPE_RAND_U1_WALL_SOURCE:
    {
      multi1d<LatticeColorVector> tmp_color_vec(Nc);

      if(params.param.source_type == SRC_TYPE_RAND_Z2_WALL_SOURCE) {
	//multi1d<int> crd(Nd); crd = 0;
	LatticeComplex z;
	LatticeReal r1, r2;
	gaussian(r1);
	//cout << peekSite(r1,crd) << endl;
	r1 /= sqrt(2.)*fabs(r1);
	//cout << peekSite(r1,crd) << endl;
	gaussian(r2);
	//cout << peekSite(r2,crd) << endl;
	r2 /= sqrt(2.)*fabs(r2);
	//cout << peekSite(r2,crd) << endl;
	//z = r1;
	//cout << peekSite(z,crd) << endl;
	//z += timesI(r2); 
	//cout << peekSite(z,crd) << endl;
	z = cmplx(r1,0) + timesI(r2); // cmplx used to avoid bug in older qdp++
	//cout << peekSite(z,crd) << endl;
	//r1 = sqrt(localNorm2(z));
	//cout << peekSite(r1,crd) << endl;
	for(int i=0; i<Nc; i++) {
	  tmp_color_vec[i] = zero;
	  pokeColor(tmp_color_vec[i], z, i);
	}
      } else if(params.param.source_type == SRC_TYPE_RAND_U1_WALL_SOURCE) {
	//multi1d<int> crd(Nd); crd = 0;
	LatticeComplex z;
	LatticeReal r;
	gaussian(z);
	//cout << peekSite(z,crd) << endl;
	r = sqrt(localNorm2(z));
	//cout << peekSite(r,crd) << endl;
	z /= r;
	//cout << peekSite(z,crd) << endl;
	//r = sqrt(localNorm2(z));
	//cout << peekSite(r,crd) << endl;
	for(int i=0; i<Nc; i++) {
	  tmp_color_vec[i] = zero;
	  pokeColor(tmp_color_vec[i], z, i);
	}
      }

      for(int color_source = 0; color_source < Nc; ++color_source)
      {
	QDPIO::cout << "color = " << color_source << endl; 

	LatticeColorVector src_color_vec = zero;

	if( (params.param.source_type == SRC_TYPE_RAND_Z2_WALL_SOURCE) ||
	    (params.param.source_type == SRC_TYPE_RAND_U1_WALL_SOURCE) ) {
	  int mu = params.param.j_decay;
	  int slice = params.param.t_source[params.param.j_decay];
	  src_color_vec = where( Layout::latticeCoordinate(mu) == slice,
				 tmp_color_vec[color_source],
				 LatticeColorVector(zero));
	} else {
	  // Make a point source at coordinates t_source
	  srcfil(src_color_vec, params.param.t_source, color_source);
	}

	// Smear the colour source if specified
	if(params.param.source_type == SRC_TYPE_SHELL_SOURCE)
	{
	  // There should be a call to maksrc2 or some-such for general source smearing

	  // displace the point source first, then smear
	  // displacement has to be taken along negative direction.
	  displacement(u_smr,src_color_vec,
		       (-1)*params.param.disp_length, params.param.disp_dir);

	  if(params.param.wave_state == WAVE_TYPE_P_WAVE)
	    p_src(u_smr, src_color_vec, params.param.direction);

	  if(params.param.wave_state == WAVE_TYPE_D_WAVE)   /* added */
	    d_src(u_smr, src_color_vec, params.param.direction);

	  gausSmear(u_smr, src_color_vec, 
		    params.param.sourceSmearParam.wvf_param, 
		    params.param.sourceSmearParam.wvfIntPar, 
		    params.param.j_decay);

	  laplacian(u_smr, src_color_vec, 
		    params.param.j_decay, params.param.laplace_power);

	  //power = 1 for one laplacian operator
	}

	for(int spin_source = 0; spin_source < Ns; ++spin_source)
	{
	  QDPIO::cout << "spin = " << spin_source << endl; 

	  // Insert a ColorVector into spin index spin_source
	  // This only overwrites sections, so need to initialize first
	  LatticeFermion chi = zero;

	  CvToFerm(src_color_vec, chi, spin_source);
      

	  /*
	   *  Move the source to the appropriate components
	   *  of quark source.
	   */
	  FermToProp(chi, quark_source, color_source, spin_source);
	}
      }
    }
    break;

    case SRC_TYPE_WALL_SOURCE:
    {
      for(int color_source = 0; color_source < Nc; ++color_source)
      {
	for(int spin_source = 0; spin_source < Ns; ++spin_source)
	{
	  // Wall fill a fermion source. Insert it into the propagator source
	  LatticeFermion chi;
	  walfil(chi, 
		 params.param.t_source[params.param.j_decay], 
		 params.param.j_decay, 
		 color_source, spin_source);
	  FermToProp(chi, quark_source, color_source, spin_source);
	}
      }
    }
    break;

    default:
      QDPIO::cout << "Unsupported source type" << endl;
      QDP_abort(1);
    }
  

    // Sanity check - write out the norm2 of the source in the j_decay direction.
    // Use this for any possible verification.
    {
      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, params.param.j_decay);

      multi1d<Double> source_corr = sumMulti(localNorm2(quark_source),
					     phases.getSet());

      push(xml_out, "Source_correlator");
      write(xml_out, "source_corr", source_corr);
      pop(xml_out);
    }
 

    // Now write the source
    try
    {
      QDPIO::cout << "Attempt to update source" << endl;

      XMLBufferWriter file_xml;
      push(file_xml, "make_source");
      int id = 0;    // NEED TO FIX THIS - SOMETHING NON-TRIVIAL NEEDED
      write(file_xml, "id", id);
      pop(file_xml);

      XMLBufferWriter record_xml;
      push(record_xml, "MakeSource");
      write(record_xml, "PropSource", params.param);
      if( (params.param.source_type == SRC_TYPE_RAND_Z2_WALL_SOURCE) ||
          (params.param.source_type == SRC_TYPE_RAND_U1_WALL_SOURCE) ) {
        write(record_xml, "RNG", ran_seed) ;
      }
      write(record_xml, "Config_info", gauge_xml);
      pop(record_xml);
    
      // Store the source
      TheNamedObjMap::Instance().create<LatticePropagator>(params.named_obj.source_id);
      TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.source_id) = quark_source;
      TheNamedObjMap::Instance().get(params.named_obj.source_id).setFileXML(file_xml);
      TheNamedObjMap::Instance().get(params.named_obj.source_id).setRecordXML(record_xml);

      QDPIO::cout << "Source successfully update" << endl;
    }
    catch (std::bad_cast)
    {
      QDPIO::cerr << InlineMakeSourceEnv::name << ": dynamic cast error" 
		  << endl;
      QDP_abort(1);
    }
    catch (const string& e) 
    {
      QDPIO::cerr << InlineMakeSourceEnv::name << ": error message: " << e << endl;
      QDP_abort(1);
    }
    
    pop(xml_out);  // make_source

    snoop.stop();
    QDPIO::cout << InlineMakeSourceEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << InlineMakeSourceEnv::name << ": ran successfully" << endl;

    END_CODE();
  } 

};
