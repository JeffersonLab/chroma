// -*- C++ -*-
// $Id: schr_dirich_fermbc_w.h,v 2.1 2006-03-16 03:00:12 edwards Exp $
/*! \file
 *  \brief Schroedinger BC - dirichlet ferm BC
 */

#ifndef __schr_dirich_fermbc_w_h__
#define __schr_dirich_fermbc_w_h__

#include "actions/gauge/gaugebcs/schr_dirich_gaugebc.h"
#include "actions/ferm/fermbcs/schr_sf_fermbc_w.h"
#include "actions/ferm/fermbcs/schr_fermbc_params_w.h"

namespace Chroma 
{ 

  /*! @ingroup fermbcs */
  namespace SchrDirichletFermBCEnv 
  {
    extern const std::string name;
    extern const bool registered;
  }


  //! Concrete class for Schroedinger BC - use for nonpertubative tuning
  /*! @ingroup fermbcs
   *
   *  Schroedinger BC for ferm actions
   */
  template<class T>
  class SchrDirichletFermBC : public SchrSFFermBC<T>
  {
  public:
    //! Only full constructor
    SchrDirichletFermBC(const SchrDirichletGaugeBC& gaugebc, 
			const SchrFermBCParams& p) : param(p)
    {
      START_CODE();

      // Initialize boundary fields
      initBnd(fld, mask, maskF, gaugebc.SFBndFld(), gaugebc.lSFmask());

      END_CODE();
    }

    //! Destructor is automatic
    ~SchrDirichletFermBC() {}

    //! Decay direction
    int getDir() const {return param.decay_dir;}

  protected:
    //! Construct the mask and boundary fields
    void initBnd(multi1d<LatticeColorMatrix>& SFBndFl,
		 multi1d<LatticeBoolean>& lSFmas,
		 LatticeBoolean& lSFmasF,
		 const multi1d<LatticeColorMatrix>& SFBndFldG,
		 const multi1d<LatticeBoolean>& lSFmaskG) const
    {
      START_CODE();

      SFBndFl = SFBndFldG;
      lSFmas = lSFmaskG;

      /* Dirichlet boundary condiditons in all Nd directions for the fermions */
      LatticeBoolean lbtest = false;

      /* The fermion mask is set on x,y,z,t = 0 and L-1 */
      for(int mu = 0; mu < Nd; mu++)
      {
	LatticeInteger litmp = Layout::latticeCoordinate(mu);
	lbtest |= (litmp == 0);
	lbtest |= (litmp == (QDP::Layout::lattSize()[mu]-1));
      }
      lSFmasF = lbtest;

      END_CODE();
    }

  protected:
    //! Mask which lattice sites have fixed ferm sites
    const LatticeBoolean& lSFmaskF() const {return maskF;}

    //! Mask which lattice sites have fixed gauge links
    const multi1d<LatticeBoolean>& lSFmask() const {return mask;}

    //! Fixed gauge links on only the lSFmask() sites
    const multi1d<LatticeColorMatrix>& SFBndFld() const {return fld;}

  protected:
    //! Get the angles on the boundaries
    const multi1d<Real>& getTheta() const {return param.theta;}

  private:
    // Hide default constuctor
    SchrDirichletFermBC() {}
    void operator=(const SchrDirichletFermBC&) {}

  private:
    SchrFermBCParams             param;
    LatticeBoolean               maskF;
    multi1d<LatticeBoolean>      mask;
    multi1d<LatticeColorMatrix>  fld;
  };
  
}

#endif
