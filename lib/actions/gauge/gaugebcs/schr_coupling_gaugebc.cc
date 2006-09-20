// $Id: schr_coupling_gaugebc.cc,v 3.1 2006-09-20 20:28:01 edwards Exp $
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
    GaugeBC< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >* createGaugeBC(XMLReader& xml, 
										       const string& path)
    {
      return new SchrCouplingGaugeBC(SchrGaugeBCParams(xml, path));
    }

    const std::string name = "SCHROEDINGER_COUPLING_GAUGEBC";

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheGaugeBCFactory::Instance().registerObject(name, createGaugeBC);
	registered = true;
      }
      return success;
    }
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
