//  $Id: sftmom.h,v 1.10 2005-01-14 18:42:38 edwards Exp $
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
  SftMom(int mom2_max, bool avg_equiv_mom=false, int j_decay=-1);

  //! Construct around some fixed mom_offset
  SftMom(int mom2_max, multi1d<int> mom_offset,
         bool avg_equiv_mom=false, int j_decay=-1) 
  { init(mom2_max, mom_offset, avg_equiv_mom, j_decay); }

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

  //! Convert momenta id to actual array of momenta
  multi1d<int> numToMom(int mom_num) const { return mom_list[mom_num]; }

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
            bool avg_equiv_mom=false, int j_decay=-1);

  multi2d<int> mom_list;
  int decay_dir;
  int num_mom;
  multi1d<LatticeComplex> phases;
  UnorderedSet sft_set;
};

}  // end namespace Chroma

#endif
