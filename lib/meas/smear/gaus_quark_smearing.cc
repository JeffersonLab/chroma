// $Id: gaus_quark_smearing.cc,v 2.3 2005-11-07 22:46:46 edwards Exp $
/*! \file
 *  \brief Gaussian smearing of color vector
 */

#include "chromabase.h"

#include "meas/smear/quark_smearing_factory.h"
#include "meas/smear/gaus_quark_smearing.h"
#include "meas/smear/gaus_smear.h"

namespace Chroma 
{

  //! Hooks to register the class
  namespace GausQuarkSmearingEnv
  {
    //! Callback function
    QuarkSmearing<LatticePropagator>* createProp(XMLReader& xml_in,
						 const std::string& path)
    {
      return new GausQuarkSmearing<LatticePropagator>(GausQuarkSmearingParams(xml_in, path));
    }

    //! Callback function
    QuarkSmearing<LatticeFermion>* createFerm(XMLReader& xml_in,
					      const std::string& path)
    {
      return new GausQuarkSmearing<LatticeFermion>(GausQuarkSmearingParams(xml_in, path));
    }
    
    //! Callback function
    QuarkSmearing<LatticeColorVector>* createColorVec(XMLReader& xml_in,
						      const std::string& path)
    {
      return new GausQuarkSmearing<LatticeColorVector>(GausQuarkSmearingParams(xml_in, path));
    }
    
    //! Name to be used
    const std::string name = "GAUGE_INV_GAUSSIAN";

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= Chroma::ThePropSmearingFactory::Instance().registerObject(name, createProp);
      foo &= Chroma::TheFermSmearingFactory::Instance().registerObject(name, createFerm);
      foo &= Chroma::TheColorVecSmearingFactory::Instance().registerObject(name, createColorVec);
      return foo;
    }

    //! Register the source construction
    const bool registered = registerAll();
  }


  //! Parameters for running code
  GausQuarkSmearingParams::GausQuarkSmearingParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "wvf_param", param.wvf_param);
    read(paramtop, "wvfIntPar", param.wvfIntPar);
    read(paramtop, "no_smear_dir", param.no_smear_dir);
  }

  // Read parameters
  void read(XMLReader& xml, const string& path, GausQuarkSmearingParams& param)
  {
    GausQuarkSmearingParams tmp(xml, path);
    param = tmp;
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const string& path, const GausQuarkSmearingParams::Param_t& param)
  {
    push(xml, path);
    
    write(xml, "wvf_param", param.wvf_param);
    write(xml, "wvfIntPar", param.wvfIntPar);
    write(xml, "no_smear_dir", param.no_smear_dir);

    pop(xml);
  }


  //! Smear the quark
  void
  GausQuarkSmearing<LatticePropagator>::operator()(LatticePropagator& quark,
						   const multi1d<LatticeColorMatrix>& u) const
  {
    gausSmear(u, quark, params.param.wvf_param, params.param.wvfIntPar, params.param.no_smear_dir);
  }

  //! Smear the quark
  void
  GausQuarkSmearing<LatticeFermion>::operator()(LatticeFermion& quark,
						const multi1d<LatticeColorMatrix>& u) const
  {
    gausSmear(u, quark, params.param.wvf_param, params.param.wvfIntPar, params.param.no_smear_dir);
  }

  //! Smear the color-vector
  void
  GausQuarkSmearing<LatticeColorVector>::operator()(LatticeColorVector& quark,
						    const multi1d<LatticeColorMatrix>& u) const
  {
    gausSmear(u, quark, params.param.wvf_param, params.param.wvfIntPar, params.param.no_smear_dir);
  }


}  // end namespace Chroma
