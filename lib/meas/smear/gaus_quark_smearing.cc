// $Id: gaus_quark_smearing.cc,v 2.5 2005-11-16 02:34:58 edwards Exp $
/*! \file
 *  \brief Gaussian smearing of color vector
 */

#include "chromabase.h"

#include "meas/smear/quark_smearing_factory.h"
#include "meas/smear/gaus_quark_smearing.h"
#include "meas/smear/gaus_smear.h"

namespace Chroma 
{

  // Read parameters
  void read(XMLReader& xml, const string& path, GausQuarkSmearingEnv::Params& param)
  {
    GausQuarkSmearingEnv::Params tmp(xml, path);
    param = tmp;
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const string& path, const GausQuarkSmearingEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Hooks to register the class
  namespace GausQuarkSmearingEnv
  {
    //! Callback function
    QuarkSmearing<LatticePropagator>* createProp(XMLReader& xml_in,
						 const std::string& path)
    {
      return new QuarkSmear<LatticePropagator>(Params(xml_in, path));
    }

    //! Callback function
    QuarkSmearing<LatticeFermion>* createFerm(XMLReader& xml_in,
					      const std::string& path)
    {
      return new QuarkSmear<LatticeFermion>(Params(xml_in, path));
    }
    
    //! Callback function
    QuarkSmearing<LatticeColorVector>* createColorVec(XMLReader& xml_in,
						      const std::string& path)
    {
      return new QuarkSmear<LatticeColorVector>(Params(xml_in, path));
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


    //! Parameters for running code
    Params::Params(XMLReader& xml, const string& path)
    {
      XMLReader paramtop(xml, path);

      read(paramtop, "wvf_param", wvf_param);
      read(paramtop, "wvfIntPar", wvfIntPar);
      read(paramtop, "no_smear_dir", no_smear_dir);
    }


    //! Parameters for running code
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);
    
      write(xml, "wvf_kind", GausQuarkSmearingEnv::name);
      write(xml, "wvf_param", wvf_param);
      write(xml, "wvfIntPar", wvfIntPar);
      write(xml, "no_smear_dir", no_smear_dir);

      pop(xml);
    }


    //! Smear the quark
    void
    QuarkSmear<LatticePropagator>::operator()(LatticePropagator& quark,
					      const multi1d<LatticeColorMatrix>& u) const
    {
      gausSmear(u, quark, params.wvf_param, params.wvfIntPar, params.no_smear_dir);
    }

    //! Smear the quark
    void
    QuarkSmear<LatticeFermion>::operator()(LatticeFermion& quark,
					   const multi1d<LatticeColorMatrix>& u) const
    {
      gausSmear(u, quark, params.wvf_param, params.wvfIntPar, params.no_smear_dir);
    }

    //! Smear the color-vector
    void
    QuarkSmear<LatticeColorVector>::operator()(LatticeColorVector& quark,
					       const multi1d<LatticeColorMatrix>& u) const
    {
      gausSmear(u, quark, params.wvf_param, params.wvfIntPar, params.no_smear_dir);
    }

  }  // end namespace
}  // end namespace Chroma

