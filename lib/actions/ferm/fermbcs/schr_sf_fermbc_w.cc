// $Id: schr_sf_fermbc_w.cc,v 3.1 2006-04-10 21:21:21 edwards Exp $
/*! \file
 *  \brief Schroedinger functional base class
 */

#include "actions/ferm/fermbcs/schr_sf_fermbc_w.h"

namespace Chroma 
{

  //! Starting slice in decay direction
  int SchrSFFermBC::getDecayMin() const
  {
    int tmin;
    int j_decay = getDir();

    switch (getMaxExtent())
    {
    case 1:
      tmin = 0;
      break;

    case 2:
      tmin = 1;
      break;

    default:
      QDPIO::cerr << __func__ << ": unsupport max extent" << endl;
      QDP_abort(1);
    }

    return tmin;
  }


  //! Ending slice in decay direction
  int SchrSFFermBC::getDecayMax() const
  {
    int tmax;
    int j_decay = getDir();

    switch (getMaxExtent())
    {
    case 1:
      tmax = QDP::Layout::lattSize()[j_decay] - 2;
      break;

    case 2:
      tmax = QDP::Layout::lattSize()[j_decay] - 3;
      break;

    default:
      QDPIO::cerr << __func__ << ": unsupport max extent" << endl;
      QDP_abort(1);
    }

    return tmax;
  }


  //! Construct the mask and boundary fields
  void SchrSFFermBC::initBnd(multi1d<LatticeColorMatrix>& SFBndFld,
			     multi1d<LatticeBoolean>& lSFmask,
			     LatticeBoolean& lSFmaskF,
			     const multi1d<LatticeColorMatrix>& SFBndFldG,
			     const multi1d<LatticeBoolean>& lSFmaskG) const
  {
    START_CODE();

    SFBndFld.resize(Nd);
    lSFmask.resize(Nd);

    int j_decay = getDir();

    // Set the fermion mask to one of the non-decay gauge masks
    for(int mu = 0; mu < Nd; ++mu)
    {
      if (mu == j_decay) continue;

      lSFmaskF = lSFmaskG[mu];
      break;
    }

    /* Set the fermion phases */
    /* Multiply the phase factors. It is all done here */
    /*   SFBndFldF[mu] *= exp(i*2*pi*theta(mu)/L(mu))  */
    for(int mu = 0, i = 0; mu < Nd; ++mu)
    {
      SFBndFld[mu] = SFBndFldG[mu];
      lSFmask[mu] = lSFmaskG[mu];

      if (mu != j_decay)
      {
	Real ftmp = Chroma::twopi * getTheta()[i] / Real(QDP::Layout::lattSize()[mu]);
	SFBndFld[mu] *= cmplx(cos(ftmp),sin(ftmp));
	++i;
      }
    }
    
    END_CODE();
  }

}
