// $Id: jacobi_quark_smearing.cc,v 3.1 2009-02-20 15:10:24 edwards Exp $
/*! \file
 *  \brief Jacobi smearing of color vector
 */

#include "chromabase.h"

#include "meas/smear/quark_smearing_factory.h"
#include "meas/smear/jacobi_quark_smearing.h"
#include "meas/smear/jacobi_smear.h"

namespace Chroma 
{

  // Read parameters
  void read(XMLReader& xml, const string& path, JacobiQuarkSmearingEnv::Params& param)
  {
    JacobiQuarkSmearingEnv::Params tmp(xml, path);
    param = tmp;
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const string& path, const JacobiQuarkSmearingEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Hooks to register the class
  namespace JacobiQuarkSmearingEnv
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
      const std::string name = "GAUGE_INV_JACOBI";
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

      read(paramtop, "wvf_param", kappa);
      read(paramtop, "wvfIntPar", iter);
      read(paramtop, "no_smear_dir", no_smear_dir);
    }


    //! Parameters for running code
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);
    
      write(xml, "wvf_kind", JacobiQuarkSmearingEnv::getName());
      write(xml, "wvf_param", kappa);
      write(xml, "wvfIntPar", iter);
      write(xml, "no_smear_dir", no_smear_dir);

      pop(xml);
    }


    //! Smear the quark
    template<>
    void
    QuarkSmear<LatticePropagator>::operator()(LatticePropagator& quark,
					      const multi1d<LatticeColorMatrix>& u) const
    {
      jacobiSmear(u, quark, params.kappa, params.iter, params.no_smear_dir);
    }

    //! Smear the quark
    template<>
    void
    QuarkSmear<LatticeStaggeredPropagator>::operator()(LatticeStaggeredPropagator& quark,
						       const multi1d<LatticeColorMatrix>& u) const
    {
      jacobiSmear(u, quark, params.kappa, params.iter, params.no_smear_dir);
    }

    //! Smear the quark
    template<>
    void
    QuarkSmear<LatticeFermion>::operator()(LatticeFermion& quark,
					   const multi1d<LatticeColorMatrix>& u) const
    {
      jacobiSmear(u, quark, params.kappa, params.iter, params.no_smear_dir);
    }

    //! Smear the color-vector
    template<>
    void
    QuarkSmear<LatticeColorVector>::operator()(LatticeColorVector& quark,
					       const multi1d<LatticeColorMatrix>& u) const
    {
      jacobiSmear(u, quark, params.kappa, params.iter, params.no_smear_dir);
    }

  }  // end namespace
}  // end namespace Chroma

