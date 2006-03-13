// $Id: schr_triv_gaugebc.cc,v 2.1 2006-03-13 05:19:01 edwards Exp $
/*! \file
 *  \brief Schroedinger functional trivial gauge BC
 */

#include "actions/gauge/gaugebcs/schr_triv_gaugebc.h"
#include "actions/gauge/gaugebcs/gaugebc_factory.h"

namespace Chroma 
{

  namespace SchrTrivialGaugeBCEnv 
  { 
    //! Callback function to register with the factory
    GaugeBC* createGaugeBC(XMLReader& xml, const string& path)
    {
      return new SchrTrivialGaugeBC(SchrGaugeBCParams(xml, path));
    }

    const std::string name = "SCHROEDINGER_TRIVIAL_GAUGEBC";
    const bool registered = TheGaugeBCFactory::Instance().registerObject(name,
									 createGaugeBC);
  }


  // Only full constructor
  SchrTrivialGaugeBC::SchrTrivialGaugeBC(const SchrGaugeBCParams& p) : 
    param(p)
  {
    // Initialize the phases
    initPhases();

    // Initialize the boundary fields
    initBnd(fld, mask);
  }


  // Initialize the phases
  void SchrTrivialGaugeBC::initPhases()
  {
    phases.lower.resize(Nc);
    phases.upper.resize(Nc);

    switch (Nc)
    {
    case 2:
      /* lower and upper boundary are trivial */
      phases.lower = QDP::zero;
      phases.upper = QDP::zero;
      break;

    case 3:
      /* lower and upper boundary are trivial */
      phases.lower = QDP::zero;
      phases.upper = QDP::zero;
      break;

    default:
      QDPIO::cerr << __func__ << ": unsupport number of colors" << endl;
      QDP_abort(1);
    }
  }


}
