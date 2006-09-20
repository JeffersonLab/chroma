// $Id: schr_chromomag_gaugebc.cc,v 3.1 2006-09-20 20:28:01 edwards Exp $
/*! \file
 *  \brief Schroedinger BC - chromo-magnetic gauge BC
 */

#include "actions/gauge/gaugebcs/schr_chromomag_gaugebc.h"
#include "actions/gauge/gaugebcs/gaugebc_factory.h"

namespace Chroma 
{

  namespace SchrChromoMagGaugeBCEnv 
  { 
    //! Callback function to register with the factory
    GaugeBC< multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >* createGaugeBC(XMLReader& xml, 
										       const string& path)
    {
      return new SchrChromoMagGaugeBC(SchrGaugeBCParams(xml, path));
    }

    const std::string name = "SCHROEDINGER_CHROMOMAG_GAUGEBC";

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
  SchrChromoMagGaugeBC::SchrChromoMagGaugeBC(const SchrGaugeBCParams& p) : 
    param(p)
  {
    START_CODE();

    fld.resize(Nd);
    mask.resize(Nd);

    int j_decay = p.decay_dir;
    int igluetmp = p.loop_extent;

    /* First construct the masks indicating the boundaries */
    LatticeInteger litmp = Layout::latticeCoordinate(j_decay);

    LatticeBoolean lbtest;
    switch (igluetmp)
    {
    case 2:
      /*  if (coord(j_decay) == 0 || coord(j_decay) == latt_cb_size(j_decay)-2) then */
      lbtest = (litmp == 0);
      lbtest |= (litmp == (QDP::Layout::lattSize()[j_decay]-2));
      /*  endif */
      break;

    case 1:
      lbtest = false;
      break;

    default:
      QDPIO::cerr << "SchrSFGaugeBC: unsupported igluetmp = " << igluetmp << endl;
      QDP_abort(1);
    }

    /*  if (coord(j_decay) == latt_cb_size(j_decay)-1) then */
    lbtest |= (litmp == (QDP::Layout::lattSize()[j_decay]-1));
    /*  endif */
  
    mask[j_decay] = lbtest;

    /* The remaining spatial directions */
    switch(igluetmp)
    {
    case 2:
      /*  if (coord(j_decay) == 1) then */
      lbtest |= (litmp ==1);
      /*  endif */
      break;

    case 1:
      break;
    }

    /*  if (coord(j_decay) == 0) then */
    lbtest |= (litmp == 0);
    /*  endif */

    for(int mu = 0; mu < Nd; ++mu)
    {
      if (mu == j_decay) continue;

      mask[mu] = lbtest;
    }
//    maskF = lbtest;

    
    /*
     * Construct the boundary fields 
     */

    /* This is a constant chromomagnetic backround field */
    if (Nd < 3)
      QDP_error_exit("Not enough dimensions for chromomagnetic backround field");

    if (Nc < 2)
      QDP_error_exit("Unsupported number of colors");

    int var_dir;
    if (j_decay == (Nd-1))
      var_dir = Nd - 2;
    else
      var_dir = Nd - 1;

    for(int mu = 1; mu < Nd; ++mu)
      fld[mu] = 1;

    Real ftmp = Chroma::twopi * p.SchrPhiMult / Real(QDP::Layout::lattSize()[var_dir]);
    LatticeReal lftmp = ftmp * Layout::latticeCoordinate(var_dir);

    fld[0] = 1.0;
    pokeColor(fld[0], cmplx(cos(lftmp), sin(lftmp)), 0, 0);
    pokeColor(fld[0], cmplx(cos(lftmp),-sin(lftmp)), 0, 0);

    END_CODE();
  }


}
