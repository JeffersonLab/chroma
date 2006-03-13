// $Id: schr_coupling_gaugebc.cc,v 2.1 2006-03-13 05:19:01 edwards Exp $
/*! \file
 *  \brief Schroedinger BC - use for coupling determinations
 */

#include "actions/gauge/gaugebcs/schr_coupling_gaugebc.h"
#include "actions/gauge/gaugebcs/gaugebc_factory.h"

namespace Chroma 
{

  namespace SchrCouplingGaugeBCEnv 
  { 
    //! Callback function to register with the factory
    GaugeBC* createGaugeBC(XMLReader& xml, const string& path)
    {
      return new SchrCouplingGaugeBC(SchrGaugeBCParams(xml, path));
    }

    const std::string name = "SCHROEDINGER_COUPLING_GAUGEBC";
    const bool registered = TheGaugeBCFactory::Instance().registerObject(name,
									 createGaugeBC);
  }


  // Only full constructor
  SchrCouplingGaugeBC::SchrCouplingGaugeBC(const SchrGaugeBCParams& p) : 
    param(p)
  {
    // Initialize the phases
    initPhases();

    // Initialize the boundary fields
    initBnd(fld, mask);
  }


  // Initialize the phases
  void SchrCouplingGaugeBC::initPhases()
  {
    phases.lower.resize(Nc);
    phases.upper.resize(Nc);

    Real ftmp = Chroma::twopi * 0.5 * SchrPhiMult();

    switch (Nc)
    {
    case 2:                     /*  SF coupling boundary conditions:  */
      ftmp /= Real(4);

      /* lower boundary */
      phases.lower[0] = -ftmp;
      phases.lower[1] =  ftmp;

      /* upper boundary */
      phases.upper[0] = -Real(3) * ftmp;
      phases.upper[1] =  Real(3) * ftmp;
      break;

    case 3:                     /*  SF coupling boundary conditions:  */
      ftmp /= Real(3);

      /* lower boundary */
      phases.lower[0] = -ftmp;
      phases.lower[1] = 0;
      phases.lower[2] =  ftmp;

      /* upper boundary */
      phases.upper[0] = -Real(3) * ftmp;
      phases.upper[1] =            ftmp;
      phases.upper[2] =  Real(2) * ftmp;
      break;

    default:
      QDPIO::cerr << __func__ << ": unsupport number of colors" << endl;
      QDP_abort(1);
    }
  }


}
