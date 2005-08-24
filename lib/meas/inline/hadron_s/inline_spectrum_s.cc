// $Id: inline_spectrum_s.cc,v 1.2 2005-08-24 16:39:23 mcneile Exp $
/*! \file
 * \brief Inline construction of staggered spectrum
 *
 * Spectrum calculations
 */

#include "meas/inline/hadron_s/inline_spectrum_s.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "meas/smear/ape_smear.h"
#include "meas/smear/sink_smear2.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "io/param_io.h"

// staggered stuff
#include "handle.h"
#include "state.h"
#include "actions/ferm/fermbcs/fermbcs.h"
#include "actions/ferm/fermacts/fermacts_s.h"

#include "meas/inline/make_xml_file.h"

namespace Chroma 
{ 
  namespace InlineSpectrumEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineSpectrum_s(InlineSpectrumParams_s(xml_in, path));
    }

    const std::string name = "SPECTRUM_S";
    const bool registered = TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
  };



  //! Reader for parameters
  void read(XMLReader& xml, const string& path, InlineSpectrumParams_s::Param_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    //    switch (version) 


    read(paramtop, "Meson_local", param.Meson_local);
    read(paramtop, "Baryon_local", param.Baryon_local);
    read(paramtop, "disconnected_local", param.disconnected_local);

    read(paramtop, "boundary", param.boundary);
    read(paramtop, "t_srce", param.t_srce);
    read(paramtop, "nrow", param.nrow);
  }


  //! Writer for parameters
  void write(XMLWriter& xml, const string& path, const InlineSpectrumParams_s::Param_t& param)
  {
    push(xml, path);

    int version = 1;
    write(xml, "version", version);

    write(xml, "Meson_local", param.Meson_local);
    write(xml, "Baryon_local", param.Baryon_local);
    write(xml, "disconnected_local", param.disconnected_local);
    write(xml, "boundary", param.boundary);
    write(xml, "nrow", param.nrow);
    write(xml, "t_srce", param.t_srce);


    pop(xml);
  }


  //! Propagator generation params input
  void read(XMLReader& xml, const string& path, InlineSpectrumParams_s::Quark_Prop_t& input)
  {
    XMLReader inputtop(xml, path);

    //    read(inputtop, "prop_inversion", input.prop_param);
    read(inputtop, "Mass", input.Mass);
    read(inputtop, "u0", input.u0);
    input.invParam.invType = CG_INVERTER;   
    read(inputtop, "RsdCG", input.invParam.RsdCG);
    read(inputtop, "MaxCG", input.invParam.MaxCG);

  }

  //! Propagator output
  void write(XMLWriter& xml, const string& path, 
           const InlineSpectrumParams_s::Quark_Prop_t& input)
  {
    push(xml, path);
    write(xml, "Mass", input.Mass);
    write(xml,"Inverter",input.invParam.invType);
    write(xml,"RsdCG", input.invParam.RsdCG);
    write(xml,"MaxCG", input.invParam.MaxCG);
    pop(xml);
  }


  // Param stuff
  InlineSpectrumParams_s::InlineSpectrumParams_s() { frequency = 0; }

  InlineSpectrumParams_s::InlineSpectrumParams_s(XMLReader& xml_in, const std::string& path) 
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
      read(paramtop, "Inversion", prop_param);

      // Possible alternate XML file pattern
      if (paramtop.count("xml_file") != 0) 
      {
	read(paramtop, "xml_file", xml_file);
      }
    }
    catch(const std::string& e) 
    {
      QDPIO::cerr << "Caught Exception reading XML: " << e << endl;
      QDP_abort(1);
    }
  }


  void
  InlineSpectrumParams_s::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    Chroma::write(xml_out, "Param", param);
    Chroma::write(xml_out, "Inversion", prop_param);
    QDP::write(xml_out, "xml_file", xml_file);

    pop(xml_out);
  }


  // Function call
  void 
  InlineSpectrum_s::operator()(const multi1d<LatticeColorMatrix>& u,
			     XMLBufferWriter& gauge_xml,
			     unsigned long update_no,
			     XMLWriter& xml_out) 
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "spectrum_s");
      write(xml_out, "update_no", update_no);
      write(xml_out, "xml_file", xml_file);
      pop(xml_out);

      XMLFileWriter xml(xml_file);
      func(u, gauge_xml, update_no, xml);
    }
    else
    {
      func(u, gauge_xml, update_no, xml_out);
    }
  }


  // Real work done here
  void 
  InlineSpectrum_s::func(const multi1d<LatticeColorMatrix>& u,
		       XMLBufferWriter& gauge_xml,
		       unsigned long update_no,
		       XMLWriter& xml_out) 
  {
    QDPIO::cout << "SPECTRUM_S: Spectroscopy for Staggered-like fermions" << endl;
    QDPIO::cout << "Gauge group: SU(" << Nc << ")" << endl;

    push(xml_out, "spectrum_s");
    write(xml_out, "update_no", update_no);


    /*
     * Sanity checks
     */


    QDPIO::cout << "Volume: " << params.param.nrow[0];
    for (int i=1; i<Nd; ++i) {
      QDPIO::cout << " x " << params.param.nrow[i];
    }
    QDPIO::cout << endl;

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.write(xml_out, "Input");

    // Write out the config info
    write(xml_out, "Config_info", gauge_xml);

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 12);
    pop(xml_out);


    // First calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);


    // generate the local quark propagator

    // Create a fermion BC
    Handle< FermBC<LatticeStaggeredFermion> >  fbc(new SimpleFermBC<LatticeStaggeredFermion>(params.param.boundary));

    //
    // Initialize fermion action
    //
    AsqtadFermAct S_f(fbc, params.prop_param.Mass, 
		      params.prop_param.u0);
    Handle<const ConnectState > state(S_f.createState(u));
    Handle<const SystemSolver<LatticeStaggeredFermion> > qprop(S_f.qprop(state,
									 params.prop_param.invParam));



    pop(xml_out);  // spectrum_s
    QDPIO::cout << "Staggered spectroscopy ran successfully" << endl;

    END_CODE();
  } 

};
