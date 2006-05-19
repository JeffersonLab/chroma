// $Id: no_quark_smearing.cc,v 3.1 2006-05-19 15:05:01 edwards Exp $
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
    const std::string name = "NONE";

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
    }


    //! Parameters for running code
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);
    
      write(xml, "wvf_kind", NoQuarkSmearingEnv::name);

      pop(xml);
    }


    //! Do not smear the quark
    template<>
    void
    QuarkSmear<LatticePropagator>::operator()(LatticePropagator& quark,
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

