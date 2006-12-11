// $Id: baryon_spinmat_funcmap_w.cc,v 3.1 2006-12-11 17:20:34 edwards Exp $
/*! \file
 *  \brief All baryon spin matrix contraction objects
 */

#include "meas/hadron/baryon_spinmat_funcmap_w.h"
#include "meas/hadron/barspinmat_w.h"

namespace Chroma
{

  // Registration aggregator
  namespace BaryonSpinMatrixEnv
  {
    /*!
     * \ingroup hadron
     *
     * @{
     */
    namespace
    {
      //! Local registration flag
      bool registered = false;


      //! NR = (1/2)* ( 1 + g_4 )
      SpinMatrix createNR(XMLReader& xml_in,
			  const std::string& path)
      {
	return BaryonSpinMats::NR();
      }

      //! NRnegPar = (1/2)* ( 1 - g_4 )
      SpinMatrix createNRnegPar(XMLReader& xml_in,
				const std::string& path)
      {
	return BaryonSpinMats::NRnegPar();
      }

      //! C = Gamma(10)
      SpinMatrix createC(XMLReader& xml_in,
			 const std::string& path)
      {
	return BaryonSpinMats::C();
      }

      //! C NR = (1/2)*C * ( 1 + g_4 )
      SpinMatrix createCNR(XMLReader& xml_in,
			   const std::string& path)
      {
	return BaryonSpinMats::CNR();
      }

      //! C gamma_5 gamma_4 = - Gamma(13)
      SpinMatrix createCg5g4(XMLReader& xml_in,
			     const std::string& path)
      {
	return BaryonSpinMats::Cg5g4();
      }

      //! C g_k = C gamma_k
      SpinMatrix createCgk(XMLReader& xml_in,
			   const std::string& path)
      {
	XMLReader paramtop(xml_in, path);

	int k;
	read(paramtop, "k", k);

	return BaryonSpinMats::Cgk(k);
      }

      //! C g4 g_k = C gamma_4 gamma_k
      SpinMatrix createCg4gk(XMLReader& xml_in,
			     const std::string& path)
      {
	XMLReader paramtop(xml_in, path);

	int k;
	read(paramtop, "k", k);

	return BaryonSpinMats::Cg4gk(k);
      }

      //! C g_k NR = C gamma_k (1/2)(1 + gamma_4)
      SpinMatrix createCgkNR(XMLReader& xml_in,
			     const std::string& path)
      {
	XMLReader paramtop(xml_in, path);

	int k;
	read(paramtop, "k", k);
	return BaryonSpinMats::CgkNR(k);
      }

      //! C g_5 = C gamma_5 = Gamma(5)
      SpinMatrix createCg5(XMLReader& xml_in,
			   const std::string& path)
      {
	return BaryonSpinMats::Cg5();
      }
      
      //! C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 )
      SpinMatrix createCg5NR(XMLReader& xml_in,
			     const std::string& path)
      {
	return BaryonSpinMats::Cg5NR();
      }

      //! C g_5 NR = (1/2)*C gamma_5 * ( 1 - g_4 )
      SpinMatrix createCg5NRnegPar(XMLReader& xml_in,
				   const std::string& path)
      {
	return BaryonSpinMats::Cg5NRnegPar();
      }

      //! C gamma_- = Cgm = (C gamma_-)^T
      SpinMatrix createCgm(XMLReader& xml_in,
			   const std::string& path)
      {
	return BaryonSpinMats::Cgm();
      }

      //! C gamma_4 gamma_- = Cg4m
      SpinMatrix createCg4m(XMLReader& xml_in,
			    const std::string& path)
      {
	return BaryonSpinMats::Cg4m();
      }

      //! C gamma_- NR = CgmNR = C gamma_- (1/2)(1 + gamma_4)
      SpinMatrix createCgmNR(XMLReader& xml_in,
			     const std::string& path)
      {
	return BaryonSpinMats::CgmNR();
      }



      //! T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2
      SpinMatrix createTunpol(XMLReader& xml_in,
			      const std::string& path)
      {
	return BaryonSpinMats::Tunpol();
      }

      //! T = (1 + gamma_4) / 2 = (1 - Gamma(8)) / 2
      SpinMatrix createTunpolNegPar(XMLReader& xml_in,
				    const std::string& path)
      {
	return BaryonSpinMats::TunpolNegPar();
      }

      //! T = \Sigma_3 (1 + gamma_4) / 2 = -i (Gamma(3) + Gamma(11)) / 2
      SpinMatrix createTpol(XMLReader& xml_in,
			    const std::string& path)
      {
	return BaryonSpinMats::Tpol();
      }

      //! T = \Sigma_3 (1 - gamma_4) / 2 = -i (-Gamma(3) + Gamma(11)) / 2
      SpinMatrix createTpolNegPar(XMLReader& xml_in,
				  const std::string& path)
      {
	return BaryonSpinMats::TpolNegPar();
      }

      //! T = (1 + \Sigma_3)*(1 + gamma_4) / 2   = (1 + Gamma(8) - i G(3) - i G(11)) / 2
      SpinMatrix createTmixed(XMLReader& xml_in,
			      const std::string& path)
      {
	return BaryonSpinMats::Tmixed();
      }

      //! T = (1 - \Sigma_3)*(1 - gamma_4) / 2   = (1 - Gamma(8) - i G(3) + i G(11)) / 2
      // Need to flip the spin for time reversal
      SpinMatrix createTmixedNegPar(XMLReader& xml_in,
				    const std::string& path)
      {
	return BaryonSpinMats::TmixedNegPar();
      }
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheBarSpinMatFuncMap::Instance().registerFunction(string("NR"), 
								     createNR);

	success &= TheBarSpinMatFuncMap::Instance().registerFunction(string("NRnegPar"), 
								     createNRnegPar);

	success &= TheBarSpinMatFuncMap::Instance().registerFunction(string("C"), 
								     createC);

	success &= TheBarSpinMatFuncMap::Instance().registerFunction(string("CNR"), 
								     createCNR);

	success &= TheBarSpinMatFuncMap::Instance().registerFunction(string("Cg5g4"), 
								     createCg5g4);

	success &= TheBarSpinMatFuncMap::Instance().registerFunction(string("Cgk"), 
								     createCgk);

	success &= TheBarSpinMatFuncMap::Instance().registerFunction(string("Cg4gk"), 
								     createCg4gk);

	success &= TheBarSpinMatFuncMap::Instance().registerFunction(string("CgkNR"), 
								     createCgkNR);

	success &= TheBarSpinMatFuncMap::Instance().registerFunction(string("Cg5"), 
								     createCg5);

	success &= TheBarSpinMatFuncMap::Instance().registerFunction(string("Cg5NR"), 
								     createCg5NR);

	success &= TheBarSpinMatFuncMap::Instance().registerFunction(string("Cg5NRnegPar"), 
								     createCg5NRnegPar);

	success &= TheBarSpinMatFuncMap::Instance().registerFunction(string("Cgm"), 
								     createCgm);

	success &= TheBarSpinMatFuncMap::Instance().registerFunction(string("Cg4m"), 
								     createCg4m);

	success &= TheBarSpinMatFuncMap::Instance().registerFunction(string("CgmNR"), 
								     createCgmNR);

	success &= TheBarSpinMatFuncMap::Instance().registerFunction(string("Tunpol"), 
								     createTunpol);

	success &= TheBarSpinMatFuncMap::Instance().registerFunction(string("TunpolNegPar"), 
								     createTunpolNegPar);

	success &= TheBarSpinMatFuncMap::Instance().registerFunction(string("Tpol"), 
								     createTpol);

	success &= TheBarSpinMatFuncMap::Instance().registerFunction(string("TpolNegPar"), 
								     createTpolNegPar);

	success &= TheBarSpinMatFuncMap::Instance().registerFunction(string("Tmixed"), 
								     createTmixed);

	success &= TheBarSpinMatFuncMap::Instance().registerFunction(string("TmixedNegPar"), 
								     createTmixedNegPar);

	registered = true;
      }
      return success;
    }

    /*! @} */   // end of group qprop
  }

}
