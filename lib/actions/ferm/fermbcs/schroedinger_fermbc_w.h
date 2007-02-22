// -*- C++ -*-
// $Id: schroedinger_fermbc_w.h,v 3.3 2007-02-22 21:11:45 bjoo Exp $
/*! @file
 * @brief Fermion action boundary conditions
 */

#ifndef __schroedinger_fermbc_w_h__
#define __schroedinger_fermbc_w_h__

#include "fermbc.h"

namespace Chroma
{

  //! Abstract class for all gauge action boundary conditions with Schroedinger BC
  /*! @ingroup fermbcs
   *
   *  Schroedinger BC implies periodic in dirs orthog to decay dir, and some
   *  kind of fixed BC in the decay dir.
   */
  class SchrFermBC : public FermBC<LatticeFermion,
		                   multi1d<LatticeColorMatrix>, 
		                   multi1d<LatticeColorMatrix> >
  {
  public:
    //! Virtual destructor to help with cleanup;
    virtual ~SchrFermBC() {}

    //! Modify U fields according to the fermion BC in place
    virtual void modify(multi1d<LatticeColorMatrix>& u) const
    {
      START_CODE();

      for(int mu=0; mu < u.size(); ++mu)
	copymask(u[mu], lSFmask()[mu], SFBndFld()[mu]);
    
      END_CODE();
    }

    //! Modify fermion fields in place
    virtual void modifyF(LatticeFermion& psi) const
    {
      START_CODE();

      copymask(psi, lSFmaskF(), LatticeFermion(QDP::zero));
      
      END_CODE();
    }

    //! Modify fermion fields in place under a subset
    virtual void modifyF(LatticeFermion& psi, 
			 const Subset& s) const
    {
      START_CODE();

      // Ooops, this is ignoring the subset!! Need to fix
      // but I don't think this is an error though
      copymask(psi, lSFmaskF(), LatticeFermion(QDP::zero));
    
      END_CODE();
    }

    //! Modify fermion fields in place
    virtual void modifyF(multi1d<LatticeFermion>& psi) const
    {
      QDP_error_exit("not implemented");
    }
    
    //! Modify fermion fields in place under a subset
    virtual void modifyF(multi1d<LatticeFermion>& psi, 
			 const Subset& s) const
    {
      QDP_error_exit("not implemented");
    }

    //! Zero some gauge-like field in place on the masked links
    virtual void zero(multi1d<LatticeColorMatrix>& ds_u) const
    {
      START_CODE();

      LatticeColorMatrix z = QDP::zero;

      for(int mu=0; mu < ds_u.size(); ++mu)
	copymask(ds_u[mu], lSFmask()[mu], z);
    
      END_CODE();
    }

    //! Says if there are fermion non-trivial 
    bool nontrivialP() const {return true;}

    //! Decay direction
    virtual int getDir() const = 0;

    //! Starting slice in decay direction
    virtual int getDecayMin() const = 0;

    //! Ending slice in decay direction
    virtual int getDecayMax() const = 0;

  protected:
    //! Mask which lattice sites have fixed fermion fields
    virtual const LatticeBoolean& lSFmaskF() const = 0;

    //! Mask which lattice sites have fixed gauge links
    virtual const multi1d<LatticeBoolean>& lSFmask() const = 0;

    //! Fixed gauge links on only the lSFmask() sites
    virtual const multi1d<LatticeColorMatrix>& SFBndFld() const = 0;

    //! Maximum plaquette size. This is what knows about 1x1 plaq or 1x2 rect.
    /*! \return 1 for 1x1 plaq or 2 for 1x2 rect in decay_dir */
    virtual int getMaxExtent() const = 0;
  };

}


#endif
