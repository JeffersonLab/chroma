// -*- C++ -*-
//  $Id: sftmom.h,v 2.1 2005-09-29 15:53:42 edwards Exp $
/*! \file
 *  \brief Fourier transform phase factor support
 */

#ifndef SFTMOM_INCLUDE
#define SFTMOM_INCLUDE

namespace Chroma {

//! Fourier transform phase factor support
/*!
 * \ingroup ft
 */
class SftMom
{
public:
  //! Constructor about origin
  SftMom(int mom2_max, bool avg_equiv_mom_=false, int j_decay=-1);

  //! Construct around some fixed mom_offset
  SftMom(int mom2_max, multi1d<int> mom_offset_,
         bool avg_equiv_mom_=false, int j_decay=-1) 
  { init(mom2_max, mom_offset_, avg_equiv_mom_, j_decay); }

  //! The set to be used in sumMulti
  const UnorderedSet& getSet() const { return sft_set; }

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

  //! Do a sumMulti(cf*phases,getSubset())
  multi2d<DComplex> sft(const LatticeComplex& cf) const;

  //! Do a sumMulti(cf*phases,getSubset())
  multi2d<DComplex> sft(const LatticeReal& cf) const;

private:
  SftMom() {} // hide default constructor

  void init(int mom2_max, multi1d<int> mom_offset,
            bool avg_mom_=false, int j_decay=-1);

  multi2d<int> mom_list;
  bool         avg_equiv_mom;
  int          decay_dir;
  int          num_mom;
  multi1d<int> mom_offset;
  multi1d<LatticeComplex> phases;
  UnorderedSet sft_set;
};

}  // end namespace Chroma

#endif
