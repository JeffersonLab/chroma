/*! @file
 * @brief Remez-type rational approximation
 */

#include "update/molecdyn/monomial/remez_rat_approx.h"

#include "update/molecdyn/monomial/rat_approx_factory.h"
#include "update/molecdyn/monomial/rat_approx_aggregate.h"

#include "update/molecdyn/monomial/remez.h"

namespace Chroma 
{ 

  //! Remez param
  void read(XMLReader& xml, const std::string& path, RemezRatApproxEnv::Params& param)
  {
    RemezRatApproxEnv::Params tmp(xml, path);
    param = tmp;
  }

  //! Write Parameters
  void write(XMLWriter& xml, const std::string& path, const RemezRatApproxEnv::Params& param) 
  {
    param.writeXML(xml, path);
  }

  //! Hooks to register the class
  namespace RemezRatApproxEnv
  {
    //! Callback function
    RationalApprox* createApprox(XMLReader& xml_in,
				 const std::string& path)
    {
      return new RatApprox(Params(xml_in, path));
    }

    //! Name to be used
    const std::string name = "REMEZ";

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


    //! Parameters for running code
    Params::Params(XMLReader& xml, const std::string& path)
    {
      XMLReader paramtop(xml, path);

      read(paramtop, "numPower", numPower);
      read(paramtop, "denPower", denPower);
      read(paramtop, "lowerMin", lowerMin);
      read(paramtop, "upperMax", upperMax);
      read(paramtop, "degree", degree);

      if (paramtop.count("digitPrecision") != 0)
	read(paramtop, "digitPrecision", digitPrecision);
      else
	digitPrecision = 50;
    }


    //! Parameters for running code
    void Params::writeXML(XMLWriter& xml, const std::string& path) const
    {
      push(xml, path);

      write(xml, "ratApproxType", RemezRatApproxEnv::name);
      write(xml, "numPower", numPower);
      write(xml, "denPower", denPower);
      write(xml, "lowerMin", lowerMin);
      write(xml, "upperMax", upperMax);
      write(xml, "degree", degree);
      write(xml, "digitPrecision", digitPrecision);
      
      pop(xml);
    }


    // Produce the partial-fraction-expansion (PFE) and its inverse (IPFE)
    void RatApprox::operator()(RemezCoeff_t& pfe, RemezCoeff_t& ipfe) const
    {
      START_CODE();

      unsigned long prec = abs(params.digitPrecision);
      unsigned long power_num = abs(params.numPower);
      unsigned long power_den = abs(params.denPower);

      QDPIO::cout << "GenApprox: Numerator  : " << params.numPower << "  Denominator: " << params.denPower << std::endl;
      QDPIO::cout << "           Action Degree " << params.degree << std::endl; 

      if (params.denPower <= 0)
      {
	QDPIO::cerr << name << ": invalid params" << std::endl;
	QDP_abort(1);
      }

      // Find approx to  x^abs(params.numPower/params.denPower)
      QDPIO::cout << "Compute partial fraction expansion" << std::endl;
      QDPIO::cout << "Numerator Power=" << power_num << " Denominator Power=" << power_den << std::endl;
      Remez  remez(params.lowerMin, params.upperMax, prec);
      remez.generateApprox(params.degree, power_num, power_den);

      if (params.numPower > 0)
      {
	// Find approx to  x^(params.numPower/params.denPower)
	QDPIO::cout << "Sign = +1" << std::endl;

	pfe = remez.getPFE();
	ipfe = remez.getIPFE();
      }
      else
      {
	// Find approx to  x^(-params.numPower/params.denPower)
	QDPIO::cout << "Sign = -1" << std::endl;

	pfe = remez.getIPFE();
	ipfe = remez.getPFE();
      }

      END_CODE();
    }


  }  // end namespace

  
} //end namespace Chroma


