// -*- C++ -*-
// $Id: schr_sf_fermbc_w.h,v 3.1 2006-04-10 21:21:21 edwards Exp $
/*! @file
 * @brief Schroedinger ferm boundary conditions
 */

#ifndef __schr_sf_fermbc_w_h__
#define __schr_sf_fermbc_w_h__

#include "actions/ferm/fermbcs/schroedinger_fermbc_w.h"

namespace Chroma
{

  //! Abstract class for SOME Schroedinger ferm BC
  /*! @ingroup fermbcs
   *
   *  Schroedinger BC implies periodic in dirs orthog to decay dir, and some
   *  kind of fixed BC in the decay dir.
   */
  class SchrSFFermBC : public SchrFermBC
  {
  public:
    //! Virtual destructor
    virtual ~SchrSFFermBC() {}

    //! Decay direction
    virtual int getDir() const = 0;

    //! Starting slice in decay direction
    virtual int getDecayMin() const;

    //! Ending slice in decay direction
    virtual int getDecayMax() const;

  protected:
    //! Construct the mask and boundary fields
    virtual void initBnd(multi1d<LatticeColorMatrix>& SFBndFld,
			 multi1d<LatticeBoolean>& lSFmask,
			 LatticeBoolean& lSFmaskF,
			 const multi1d<LatticeColorMatrix>& SFBndFldG,
			 const multi1d<LatticeBoolean>& lSFmaskG) const;

    //! Get the angles on the boundaries
    virtual const multi1d<Real>& getTheta() const = 0;

    //! Maximum plaquette size. This is what knows about 1x1 plaq or 1x2 rect.
    /*! \return 1 for 1x1 plaq or 2 for 1x2 rect in decay_dir */
    virtual int getMaxExtent() const = 0;
  };

}



#endif
