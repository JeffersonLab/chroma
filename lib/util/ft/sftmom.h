// -*- C++ -*-
//  $Id: sftmom.h,v 3.9 2009-02-17 16:33:38 edwards Exp $
/*! \file
 *  \brief Fourier transform phase factor support
 */

#ifndef __sftmom_h__
#define __sftmom_h__

#include "chromabase.h"

namespace Chroma 
{

  //! Param struct for SftMom
  /*!
   * \ingroup ft
   */
  struct SftMomParams_t
  {
    SftMomParams_t();                 /*!< Default constructor in sftmom.cc */
    int           mom2_max;           /*!< (mom - mom_origin)^2 <= mom2_max */
    multi1d<int>  mom_offset;         /*!< Origin for the momentum */
    bool          avg_equiv_mom;      /*!< average over equivalent momenta */
    multi1d<int>  origin_offset;      /*<! Coordinate offset of the origin. Used to fix phase factor */
    int           decay_dir;          /*!< Decay direction */
  };


  //! Fourier transform phase factor support
  /*!
   * \ingroup ft
   */
  class SftMom
  {
  public:
    //! Constructor about origin
    SftMom(int mom2_max, bool avg_equiv_mom_=false, int j_decay=-1);
    
		//! Constructor about origin, with a list of momenta 
    SftMom(const multi2d<int> & moms , int j_decay=-1);

    //! Construct around some fixed origin_offset
    SftMom(int mom2_max, multi1d<int> origin_offset_,
	   bool avg_equiv_mom_=false, int j_decay=-1) ;

    //! Construct around some fixed origin_offset and mom_offset
    SftMom(int mom2_max, multi1d<int> origin_offset_, multi1d<int> mom_offset_,
	   bool avg_equiv_mom_=false, int j_decay=-1) 
      { init(mom2_max, origin_offset_, mom_offset_, avg_equiv_mom_, j_decay); }

    //! General constructor
    SftMom(const SftMomParams_t& p)
      { init(p.mom2_max, p.origin_offset, p.mom_offset, p.avg_equiv_mom, p.decay_dir); }

    //! The set to be used in sumMulti
    const Set& getSet() const { return sft_set; }

    //! Number of momenta
    int numMom() const { return num_mom; }

    //! Number of subsets - length in decay direction
    int numSubsets() const { return sft_set.numSubsets(); }

    //! Number of sites in each subset
    int numSites() const;

    //! Decay direction
    int getDir() const { return decay_dir; }

    //! Are momenta averaged?
    bool getAvg() const { return avg_equiv_mom; }

    //! Momentum offset
    multi1d<int> getMomOffset() const { return mom_offset; }

    //! Convert momenta id to actual array of momenta
    multi1d<int> numToMom(int mom_num) const { return mom_list[mom_num]; }

    //! Convert array of momenta to momenta id
    /*! \return id in [0,numMom()-1] or -1 if not in list */
    int momToNum(const multi1d<int>& mom_in) const;

    //! Canonically order an array of momenta
    /*! \return abs(mom[0]) >= abs(mom[1]) >= ... >= abs(mom[mu]) >= ... >= 0 */
    multi1d<int> canonicalOrder(const multi1d<int>& mom) const;

    //! Return the phase for this particular momenta id
    const LatticeComplex& operator[](int mom_num) const
      { return phases[mom_num]; }

    //! Return the the multiplicity for this momenta id.
    /*! Only nonzero if momentum averaging is turned on */
    int multiplicity(int mom_num) const
      { return mom_degen[mom_num]; }

    //! Do a sumMulti(cf*phases,getSet())
    multi2d<DComplex> sft(const LatticeComplex& cf) const;

    //! Do a sum(cf*phases,getSet()[my_subset])
    multi2d<DComplex> sft(const LatticeComplex& cf, int subset_color) const;

    //! Do a sumMulti(cf*phases,getSet())
    multi2d<DComplex> sft(const LatticeReal& cf) const;

    //! Do a sumMulti(cf*phases,getSet()[my_subset])
    multi2d<DComplex> sft(const LatticeReal& cf, int subset_color) const;

#if BASE_PRECISION==32
    multi2d<DComplex> sft(const LatticeComplexD& cf) const;
    //! Do a sum(cf*phases,getSet()[my_subset])
    multi2d<DComplex> sft(const LatticeComplexD& cf, int subset_color) const;
#endif

  private:
    SftMom() {} // hide default constructor

    void init(int mom2_max, multi1d<int> origin_offset, multi1d<int> mom_offset,
	      bool avg_mom_=false, int j_decay=-1);

    multi2d<int> mom_list;
    bool         avg_equiv_mom;
    int          decay_dir;
    int          num_mom;
    multi1d<int> origin_offset;
    multi1d<int> mom_offset;
    multi1d<LatticeComplex> phases;
    multi1d<int> mom_degen;
    Set sft_set;
  };

}  // end namespace Chroma

#endif
