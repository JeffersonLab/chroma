// $Id: schr_dirich_gaugebc.cc,v 2.1 2006-03-14 04:49:54 edwards Exp $
/*! \file
 *  \brief Schroedinger BC - dirichlet gauge BC
 */

#include "actions/gauge/gaugebcs/schr_dirich_gaugebc.h"
#include "actions/gauge/gaugebcs/gaugebc_factory.h"

namespace Chroma 
{

  namespace SchrDirichletGaugeBCEnv 
  { 
    //! Callback function to register with the factory
    GaugeBC* createGaugeBC(XMLReader& xml, const string& path)
    {
      return new SchrDirichletGaugeBC(SchrGaugeBCParams(xml, path));
    }

    const std::string name = "SCHROEDINGER_DIRICHLET_GAUGEBC";
    const bool registered = TheGaugeBCFactory::Instance().registerObject(name,
									 createGaugeBC);
  }


  // Only full constructor
  SchrDirichletGaugeBC::SchrDirichletGaugeBC(const SchrGaugeBCParams& p) : 
    param(p)
  {
    START_CODE();

    fld.resize(Nd);
    mask.resize(Nd);

    int j_decay = p.decay_dir;
    int igluetmp = p.loop_extent;

    /* Dirichlet boundary condiditons in all Nd directions for the fermions */
    LatticeBoolean lbtest = false;

    /* The fermion mask is set on x,y,z,t = 0 and L-1 */
    for(int mu = 0; mu < Nd; mu++)
    {
      LatticeInteger litmp = Layout::latticeCoordinate(mu);
      lbtest |= (litmp == 0);
      lbtest |= (litmp == (QDP::Layout::lattSize()[mu]-1));
    }
//    lSFmaskF = lbtest;

    /* The gauge mask is set on x,y,z,t = 0 and L-1 for all mu and mask(cb,mu) on L-2 */
    for(int mu = 0; mu < Nd; mu++)
    {
      LatticeBoolean lbtmp = (Layout::latticeCoordinate(mu) == (QDP::Layout::lattSize()[mu]-2));
      mask[mu] = lbtest | lbtmp;
    }

    /* The boundary fields are always zero */
    fld = QDP::zero;

    END_CODE();
  }


}
