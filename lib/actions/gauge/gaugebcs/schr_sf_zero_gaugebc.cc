// $Id: schr_sf_zero_gaugebc.cc,v 3.1 2006-09-21 18:43:26 edwards Exp $
/*! \file
 *  \brief Schroedinger BC - happens to zero out gauge fields in bc_dir
 */

#include "actions/gauge/gaugebcs/schr_sf_zero_gaugebc.h"
#include "actions/gauge/gaugebcs/gaugebc_factory.h"

namespace Chroma 
{

  namespace SchrSFZeroGaugeBCEnv 
  { 
    //! Callback function to register with the factory
    GaugeBC< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >* createGaugeBC(XMLReader& xml, 
										       const string& path)
    {
      return new SchrSFZeroGaugeBC(SchrGaugeBCParams(xml, path));
    }

    const std::string name = "SCHROEDINGER_ZERO_GAUGEBC";

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
  SchrSFZeroGaugeBC::SchrSFZeroGaugeBC(const SchrGaugeBCParams& p) : 
    param(p)
  {
    // Initialize the phases
    initPhases();

    // Initialize the boundary fields
    initBnd(fld, mask);
  }


  // Initialize the phases
  void SchrSFZeroGaugeBC::initPhases()
  {
    phases.lower.resize(Nc);
    phases.upper.resize(Nc);

    phases.lower = QDP::zero;
    phases.upper = QDP::zero;
  }


  // Modify U fields in place
  void SchrSFZeroGaugeBC::modify(multi1d<LatticeColorMatrix>& u) const
  {
    START_CODE();

    LatticeColorMatrix z = QDP::zero;

    for(int mu=0; mu < u.size(); ++mu)
      copymask(u[mu], lSFmask()[mu], z);

    END_CODE();
  }

  // Zero some gauge-like field in place on the masked links
  void SchrSFZeroGaugeBC::zero(multi1d<LatticeColorMatrix>& ds_u) const
  {
    START_CODE();

    LatticeColorMatrix z = QDP::zero;

    for(int mu=0; mu < ds_u.size(); ++mu)
      copymask(ds_u[mu], lSFmask()[mu], z);

    END_CODE();
  }
}
