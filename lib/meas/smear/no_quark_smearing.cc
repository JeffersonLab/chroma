// $Id: no_quark_smearing.cc,v 3.4 2008-11-04 18:43:58 edwards Exp $
/*! \file
 *  \brief No quark smearing
 */

#include "chromabase.h"

#include "meas/smear/quark_smearing_factory.h"
#include "meas/smear/no_quark_smearing.h"

namespace Chroma 
{

  // Read parameters
  void read(XMLReader& xml, const string& path, NoQuarkSmearingEnv::Params& param)
  {
    NoQuarkSmearingEnv::Params tmp(xml, path);
    param = tmp;
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const string& path, const NoQuarkSmearingEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Hooks to register the class
  namespace NoQuarkSmearingEnv
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
      const std::string name = "NONE";
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
    }


    //! Parameters for running code
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);
    
      write(xml, "wvf_kind", NoQuarkSmearingEnv::getName());

      pop(xml);
    }


    //! Do not smear the quark
    template<>
    void
    QuarkSmear<LatticePropagator>::operator()(LatticePropagator& quark,
					      const multi1d<LatticeColorMatrix>& u) const {}

    //! Do not smear the quark
    template<>
    void
    QuarkSmear<LatticeStaggeredPropagator>::operator()(LatticeStaggeredPropagator& quark,
						       const multi1d<LatticeColorMatrix>& u) const {}

    //! Do no smear the quark
    template<>
    void
    QuarkSmear<LatticeFermion>::operator()(LatticeFermion& quark,
					   const multi1d<LatticeColorMatrix>& u) const {}

    //! Do not smear the color-vector
    template<>
    void
    QuarkSmear<LatticeColorVector>::operator()(LatticeColorVector& quark,
					       const multi1d<LatticeColorMatrix>& u) const {}

  }  // end namespace
}  // end namespace Chroma

