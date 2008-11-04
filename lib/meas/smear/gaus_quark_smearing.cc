// $Id: gaus_quark_smearing.cc,v 3.3 2008-11-04 18:43:57 edwards Exp $
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
    namespace
    {
      //! Callback function
      QuarkSmearing<LatticePropagator>* createProp(XMLReader& xml_in,
						   const std::string& path)
      {
	return new QuarkSmear<LatticePropagator>(Params(xml_in, path));
      }

      //! Callback function
      QuarkSmearing<LatticeStaggeredPropagator>* createStagProp(XMLReader& xml_in,
								const std::string& path)
      {
	return new QuarkSmear<LatticeStaggeredPropagator>(Params(xml_in, path));
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
    
      //! Local registration flag
      bool registered = false;

      //! Name to be used
      const std::string name = "GAUGE_INV_GAUSSIAN";
    }

    //! Return the name
    std::string getName() {return name;}

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::ThePropSmearingFactory::Instance().registerObject(name, createProp);
	success &= Chroma::TheStagPropSmearingFactory::Instance().registerObject(name, createStagProp);
	success &= Chroma::TheFermSmearingFactory::Instance().registerObject(name, createFerm);
	success &= Chroma::TheColorVecSmearingFactory::Instance().registerObject(name, createColorVec);
	registered = true;
      }
      return success;
    }


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
    
      write(xml, "wvf_kind", GausQuarkSmearingEnv::getName());
      write(xml, "wvf_param", wvf_param);
      write(xml, "wvfIntPar", wvfIntPar);
      write(xml, "no_smear_dir", no_smear_dir);

      pop(xml);
    }


    //! Smear the quark
    template<>
    void
    QuarkSmear<LatticePropagator>::operator()(LatticePropagator& quark,
					      const multi1d<LatticeColorMatrix>& u) const
    {
      gausSmear(u, quark, params.wvf_param, params.wvfIntPar, params.no_smear_dir);
    }

    //! Smear the quark
    template<>
    void
    QuarkSmear<LatticeStaggeredPropagator>::operator()(LatticeStaggeredPropagator& quark,
						       const multi1d<LatticeColorMatrix>& u) const
    {
      gausSmear(u, quark, params.wvf_param, params.wvfIntPar, params.no_smear_dir);
    }

    //! Smear the quark
    template<>
    void
    QuarkSmear<LatticeFermion>::operator()(LatticeFermion& quark,
					   const multi1d<LatticeColorMatrix>& u) const
    {
      gausSmear(u, quark, params.wvf_param, params.wvfIntPar, params.no_smear_dir);
    }

    //! Smear the color-vector
    template<>
    void
    QuarkSmear<LatticeColorVector>::operator()(LatticeColorVector& quark,
					       const multi1d<LatticeColorMatrix>& u) const
    {
      gausSmear(u, quark, params.wvf_param, params.wvfIntPar, params.no_smear_dir);
    }

  }  // end namespace
}  // end namespace Chroma

