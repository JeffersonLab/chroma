#include "chromabase.h"
#include "update/molecdyn/monomial/read_rat_approx.h"
#include "update/molecdyn/monomial/rat_approx_factory.h"
#include "update/molecdyn/monomial/rat_approx_aggregate.h"

namespace Chroma { 

  namespace ReadRatApproxEnv
  {
    //! Callback function
    RationalApprox* createApprox(XMLReader& xml_in,
				 const std::string& path)
    {
      return new ReadRatApprox(Params(xml_in, path));
    }

    
    //! Name to be used
    const std::string name = "READ_COEFFS";

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheRationalApproxFactory::Instance().registerObject(name, createApprox);
	registered = true;
      }
      return success;
    }

    //! Params for Remez type rational approximation
    /*! @ingroup monomial */
    Params::Params(XMLReader& in, const std::string& path) { 
      
      try { 
	XMLReader paramtop(in,path);

       
	XMLReader pfe_in(paramtop, "PFECoeffs");
	read(pfe_in, "norm", pfe.norm);
	read(pfe_in, "res", pfe.res);
	read(pfe_in, "pole", pfe.pole);
	XMLReader ipfe_in(paramtop, "IPFECoeffs");
	read(ipfe_in, "norm", ipfe.norm);
	read(ipfe_in, "res", ipfe.res);
	read(ipfe_in, "pole", ipfe.pole);
      }
      catch(const std::string& e) { 
	QDPIO::cout << "Caught Exception reading XML" << e << endl;
	QDP_abort(1);
      }
    }

    void Params::writeXML(XMLWriter& out, const std::string& path) const {
      push(out, path);
      push(out, "PFECoeffs");
      write(out, "norm", pfe.norm);
      write(out, "res", pfe.res);
      write(out, "pole", pfe.pole);
      pop(out);
      push(out, "IPFECoeffs");
      write(out, "norm", ipfe.norm);
      write(out, "res",  ipfe.res);
      write(out, "pole", ipfe.pole);
      pop(out);

    }


    
  }; // Namespace   

  void read(XMLReader &xml, const string& path, ReadRatApproxEnv::Params& param) {
    ReadRatApproxEnv::Params tmp(xml, path);
    param=tmp;
  }

 //! Write Parameters
  void write(XMLWriter& xml, const std::string& path, const ReadRatApproxEnv::Params& param) 
  {
    param.writeXML(xml, path);
  }


  void ReadRatApprox::operator()(RemezCoeff_t& pfe, RemezCoeff_t& ipfe) const 
  {
    pfe=params.pfe;
    ipfe=params.ipfe;
  }

};
