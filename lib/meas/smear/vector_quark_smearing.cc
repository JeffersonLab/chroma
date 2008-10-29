// $Id: vector_quark_smearing.cc,v 3.1 2008-10-29 19:42:37 jbulava Exp $
/*! \file
 *  \brief Gaussian smearing of color vector
 */

#include "chromabase.h"

#include "meas/smear/quark_smearing_factory.h"
#include "meas/smear/vector_quark_smearing.h"
#include "meas/smear/vector_smear.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{

  // Read parameters
  void read(XMLReader& xml, const string& path, VectorQuarkSmearingEnv::Params& param)
  {
    VectorQuarkSmearingEnv::Params tmp(xml, path);
    param = tmp;
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const string& path, const VectorQuarkSmearingEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Hooks to register the class
  namespace VectorQuarkSmearingEnv
  {
    
    //! Callback function
    QuarkSmearing<LatticeFermion>* createFerm(XMLReader& xml_in,
					      const std::string& path)
    {
      return new VectorQuarkSmear<LatticeFermion>(Params(xml_in, path));
    }
    
    //! Callback function
    QuarkSmearing<LatticeColorVector>* createColorVec(XMLReader& xml_in,
						      const std::string& path)
    {
      return new VectorQuarkSmear<LatticeColorVector>(Params(xml_in, path));
    }
    
    //! Name to be used
    const std::string name = "VECTOR_SMEAR";

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	
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

      read(paramtop, "subset_vecs_id", vecs_id);
      read(paramtop, "sigma", sigma);
      read(paramtop, "no_smear_dir", no_smear_dir);
    }


    //! Parameters for running code
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);
    
      write(xml, "wvf_kind", VectorQuarkSmearingEnv::name);
      write(xml, "subset_vecs_id", vecs_id);
      write(xml, "sigma", sigma);
      write(xml, "no_smear_dir", no_smear_dir);

      pop(xml);
    }

    //! Smear the quark
    template<>
    void
    VectorQuarkSmear<LatticeFermion>::operator()(LatticeFermion& quark,
					   const multi1d<LatticeColorMatrix>& u) const
    {
      vectorSmear(quark, vecs, params.sigma, params.no_smear_dir);
    }

    //! Smear the color-vector
    template<>
    void
    VectorQuarkSmear<LatticeColorVector>::operator()(LatticeColorVector& quark,
					       const multi1d<LatticeColorMatrix>& u) const
    {
      vectorSmear(quark, vecs, params.sigma, params.no_smear_dir);
    }

  }  // end namespace
}  // end namespace Chroma

