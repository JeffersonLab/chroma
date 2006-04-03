// $Id: schr_sf_gaugebc.cc,v 3.0 2006-04-03 04:58:54 edwards Exp $
/*! \file
 *  \brief Schroedinger functional base class
 */

#include "gaugebc.h"

#include "actions/gauge/gaugebcs/schr_sf_gaugebc.h"

namespace Chroma 
{

  //! Construct the mask and boundary fields
  void SchrSFGaugeBC::initBnd(multi1d<LatticeColorMatrix>& SFBndFld,
			      multi1d<LatticeBoolean>& lSFmask) const
  {
    START_CODE();

    SFBndFld.resize(Nd);
    lSFmask.resize(Nd);

    int j_decay = getDir();
    int igluetmp = getMaxExtent();

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
  
    lSFmask[j_decay] = lbtest;

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

      lSFmask[mu] = lbtest;
    }
//    lSFmaskF = lbtest;

    
    /*
     * Construct the boundary fields 
     */

    /* In the j_decay direction they are 1 */
    SFBndFld[j_decay] = 1;

    /* They are uniform in the spatial direction */
    /* Here I am assuming a uniform diagonal boundary and general Nc */
                    
    /* Construct the coordinates used */
    int tval;
    switch (igluetmp)
    {
    case 1:
      tval = QDP::Layout::lattSize()[j_decay] - 1;      /* construct  T-x0 */
      break;

    case 2:
      tval = QDP::Layout::lattSize()[j_decay] - 2;      /* construct  T-x0 */
      break;

    default:
      QDP_error_exit("SchrSFGaugeBC: cannot get here");
    }

    LatticeReal lftmp0(litmp);
    LatticeReal lftmp1(tval - litmp);

    /* Construct the spatial phases. Only this part is specific to Nc==2,3 */
    Phases_t phases = getPhases();
    multi1d<Real> phi0 = phases.lower;
    multi1d<Real> phi1 = phases.upper;

#if 0
    //
    // zap from here
    //
    Real ftmp = Chroma::twopi * 0.5 * SchrPhiMult();
    switch (Nc)
    {
    case 2:
      switch (SchrFun)
      {
      case 1:
	/* lower and upper boundary are trivial */
	phi0 = 0;
	phi1 = 0;
	break;

      case 2:                     /*  SF nonperturbative improvement bc's */
	ftmp /= TO_REAL(4);

	/* lower boundary */
	phi0[0] = -ftmp;
	phi0[1] =  ftmp;

	/* upper boundary */
	phi1[0] = -TO_REAL(3) * ftmp;
	phi1[1] =  TO_REAL(3) * ftmp;
	break;

      case 3:                     /*  SF coupling boundary conditions:  */
	ftmp /= TO_REAL(4);

	/* lower boundary */
	phi0[0] = -ftmp;
	phi0[1] =  ftmp;

	/* upper boundary */
	phi1[0] = -TO_REAL(3) * ftmp;
	phi1[1] =  TO_REAL(3) * ftmp;
	break;

      default:
	QDP_error_exit("Unsupported value", SchrFun);
      }
      break;

    case 3:
      switch (SchrFun)
      {
      case 1:
	/* lower and upper boundary are trivial */
	phi0 = 0;
	phi1 = 0;
	break;

      case 2:                     /*  SF nonperturbative improvement bc's */
	ftmp /= TO_REAL(6);

	/* lower boundary */
	phi0[0] = -ftmp;
	phi0[1] = 0;
	phi0[2] = ftmp;

	/* upper boundary */
	phi1[0] = -TO_REAL(5) * ftmp;
	phi1[1] = TO_REAL(2) * ftmp;
	phi1[2] = TO_REAL(3) * ftmp;
	break;

      case 3:                     /*  SF coupling boundary conditions:  */
	ftmp /= TO_REAL(3);

	/* lower boundary */
	phi0[0] = -ftmp;
	phi0[1] = 0;
	phi0[2] =  ftmp;

	/* upper boundary */
	phi1[0] = -TO_REAL(3) * ftmp;
	phi1[1] =               ftmp;
	phi1[2] =  TO_REAL(2) * ftmp;
	break;

      default:
	QDP_error_exit("Unsupported value", SchrFun);
      }
      break;

    default:
      QDP_error_exit("Unsupport number of colors", Nc);
    }
    //
    // zap to here
    //
#endif


    /* Normalize. This is generic */
    multi3d<Real> phi(Nd, Nc, 2);

    for(int mu = 0; mu < Nd; ++mu)
    {
      if (mu == j_decay) continue;

      Real ftmp = Real(1) / Real(tval * QDP::Layout::lattSize()[mu]);

      for(int i = 0; i < Nc; ++i)
      {
	phi[mu][i][0] = phi0[i] * ftmp;
	phi[mu][i][1] = phi1[i] * ftmp;
      }
    }

    /*  push(xml_out,"Setfsbc");
	write(xml_out, "tval", tval);
	write(xml_out, "phi0", phi0);
	write(xml_out, "phi1", phi1);
	write(xml_out, "phi", phi);
	write(xml_out, "lftmp0", lftmp0);
	write(xml_out, "lftmp1", lftmp1);
	pop(xml_out);  */

    /* Exponentiate the diagonal phases into the boundary fields */
              
    for(int mu = 0; mu < Nd; ++mu)
    {
      if (mu == j_decay) continue;

      SFBndFld[mu] = QDP::zero;

      LatticeReal lftmp;
      for(int i = 0; i < Nc; ++i)
      {
	lftmp  = lftmp0 * phi[mu][i][1];
	lftmp += lftmp1 * phi[mu][i][0];
	
	pokeColor(SFBndFld[mu], cmplx(cos(lftmp),sin(lftmp)), i, i);
      }
    }

    END_CODE();
  }

}
