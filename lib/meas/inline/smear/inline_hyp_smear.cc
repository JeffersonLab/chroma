// $Id: inline_hyp_smear.cc,v 1.2 2005-04-10 17:05:33 edwards Exp $
/*! \file
 *  \brief Inline Hyp smearing
 */

#include "meas/inline/smear/inline_hyp_smear.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/smear/hyp_smear.h"
#include "util/info/proginfo.h"
#include "util/gauge/unit_check.h"
#include "io/writeszin.h"

#include <sys/time.h>   // for timings

namespace Chroma 
{ 
  namespace InlineHypSmearEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineHypSmear(InlineHypSmearParams(xml_in, path));
    }

    const std::string name = "HYP_SMEAR";
    const bool registered = TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
  };


  //! Target file
  void read(XMLReader& xml, const string& path, InlineHypSmearParams::Hyp_t& input)
  {
    XMLReader inputtop(xml, path);

    input.cfg_type = CFG_TYPE_SZIN;
    read(inputtop, "hyp_file", input.hyp_file);
  }


  //! Target file
  void write(XMLWriter& xml, const string& path, const InlineHypSmearParams::Hyp_t& input)
  {
    push(xml, path);
    write(xml, "hyp_file", input.hyp_file);
    pop(xml);
  }


  //! Parameters for running code
  void read(XMLReader& xml, const string& path, InlineHypSmearParams::Param_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    switch (version) 
    {
    case 3:
      /**************************************************************************/

      /* this version allows a variable num_smear */
      read(paramtop, "num_smear", param.num_smear);

      if( param.num_smear < 0 )
      {
	QDP_error_exit( "hypsmear.cc: invalid number of hyp smearing iterations, num_smear = %d", param.num_smear );
      }

      break;

    case 2:
      /**************************************************************************/

      /* this version only allows num_smear = 1 */
      param.num_smear = 1;

      break;

    default :
      /**************************************************************************/

      QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
      QDP_abort(1);
    }

    read(paramtop, "alpha1", param.alpha1);
    read(paramtop, "alpha2", param.alpha2);
    read(paramtop, "alpha3", param.alpha3);

    read(paramtop, "nrow", param.nrow);

    /*
     *  Now information about whether to truncate the configuration
     */
    read(paramtop, "trunc", param.trunc);
    switch(param.trunc){
    case 1:
      read(paramtop, "t_start", param.t_start);
      read(paramtop, "t_end", param.t_end);
      read(paramtop, "j_decay", param.j_decay);
      break;
    default:
      break;
    }
  }

  //! Parameters
  void write(XMLWriter& xml, const string& path, const InlineHypSmearParams::Param_t& param)
  {
    push(xml, path);

    int version = 3;
    write(xml, "version", version);

    /* this version allows a variable num_smear */
    write(xml, "num_smear", param.num_smear);

    if( param.num_smear < 0 )
    {
      QDP_error_exit("hypsmear: invalid number of hyp smearing iterations, num_smear = %d", param.num_smear );
    }

    write(xml, "alpha1", param.alpha1);
    write(xml, "alpha2", param.alpha2);
    write(xml, "alpha3", param.alpha3);

    write(xml, "nrow", param.nrow);

    /*
     *  Now information about whether to truncate the configuration
     */
    write(xml, "trunc", param.trunc);
    switch(param.trunc){
    case 1:
      write(xml, "t_start", param.t_start);
      write(xml, "t_end", param.t_end);
      write(xml, "j_decay", param.j_decay);
      break;
    default:
      break;
    }
  }



  // Param stuff
  InlineHypSmearParams::InlineHypSmearParams() { frequency = 0; }

  InlineHypSmearParams::InlineHypSmearParams(XMLReader& xml_in, const std::string& path) 
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

      // Read in the hyp file info
      read(paramtop, "Hyp", hyp);
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  void
  InlineHypSmearParams::write(XMLWriter& xml, const std::string& path) 
  {
    push(xml, path);

    Chroma::write(xml, "Param", param);
    Chroma::write(xml, "Hyp", hyp);

    pop(xml);
  }


  void 
  InlineHypSmear::operator()(const multi1d<LatticeColorMatrix>& u,
			     XMLBufferWriter& gauge_xml,
			     unsigned long update_no,
			     XMLWriter& xml_out) 
  {
    push(xml_out, "hypsmear");
    write(xml_out, "update_no", update_no);
    
    QDPIO::cout << " HYPSMEAR: HYP smearing of gauge config" << endl;

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.write(xml_out, "Input");

    // Write out the config header
    write(xml_out, "Config_info", gauge_xml);

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 1);
    pop(xml_out);


    // Check if the gauge field configuration is unitarized
    clock_t t1 = clock();
    unitarityCheck(u);
    clock_t t2 = clock();
    QDPIO::cout << "Unitarity took " << (double)((int)(t2)-(int)(t1))/(double)(CLOCKS_PER_SEC) << " secs" << endl;
  

    // Calculate some gauge invariant observables just for info.
    t1 = clock();
    MesPlq(xml_out, "Observables", u);
    t2 = clock();
    QDPIO::cout << "Plaquette took " << (double)((int)(t2)-(int)(t1))/(double)(CLOCKS_PER_SEC) << " secs" << endl;


    // Now hyp smear
    multi1d<LatticeColorMatrix> u_hyp(Nd);

    Real BlkAccu = 1.0e-5;
    int BlkMax = 100;

    t1 = clock();
    if( params.param.num_smear > 0 )
    {
      for( int n = 0; n < params.param.num_smear; n ++ )
      {
	Hyp_Smear(u, u_hyp, 
		  params.param.alpha1, params.param.alpha2, params.param.alpha3, 
		  BlkAccu, BlkMax);
      }
    }
    else
    {
      u_hyp = u;
    }
    t2 = clock();
    QDPIO::cout << "Hypsmear took " << (double)((int)(t2)-(int)(t1))/(double)(CLOCKS_PER_SEC) << " secs" << endl;

    // Calculate some gauge invariant observables just for info.

    MesPlq(xml_out, "HYP_observables", u_hyp);

    // Now write the configuration to disk

    QDPIO::cout << "call szin trunc-er" << endl;

    t1 = clock();

    switch (params.hyp.cfg_type) 
    {
    case CFG_TYPE_SZIN :
    {
      SzinGauge_t szin_out;

      switch(params.param.trunc)
      {
      case 1:
	QDPIO::cout << "Call writeSzinTrunc" << endl;
	writeSzinTrunc(szin_out, u_hyp, params.param.j_decay,
		       params.param.t_start, params.param.t_end,
		       params.hyp.hyp_file);
	break;
      default:
	QDPIO::cout << "Call writeSzin" << endl;
	writeSzin(szin_out, u_hyp,
		  params.hyp.hyp_file);
	break;
      }
      break;
    }
    default :
      QDP_error_exit("Configuration type is unsupported.");
    }

    t2 = clock();
    QDPIO::cout << "Gauge write took " << (double)((int)(t2)-(int)(t1))/(double)(CLOCKS_PER_SEC) << " secs" << endl;

    pop(xml_out);

    END_CODE();
  } 

};
