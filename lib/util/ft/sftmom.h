//  $Id: sftmom.h,v 1.2 2003-03-14 17:13:44 flemingg Exp $
//  $Log: sftmom.h,v $
//  Revision 1.2  2003-03-14 17:13:44  flemingg
//  SftMom::sft() now works.
//
//  Revision 1.1  2003/03/14 05:06:06  flemingg
//  Initial version of SftMom class
//

#ifndef SFTMOM_INCLUDE
#define SFTMOM_INCLUDE

class SftMom
{
public:
  SftMom(int mom2_max, bool avg_equiv_mom, int j_decay) ;

  const Set& getSubset() const { return sft_subsets ; }

  int momToNum(const multi1d<int>& mom) ; // GTF: could be private ???

  int numMom() { return num_mom ; }

  int numSubsets() const { return sft_subsets.numSubsets() ; }

  multi1d<int> numToMom(int mom_num) const { return mom_list[mom_num] ; }

  const LatticeComplex& operator[](int mom_num) const
    { return phases[mom_num] ; }

  multi2d<DComplex> sft(const LatticeComplex& cf) ;

private:
  SftMom() {} // hide default constructor

  multi2d<int> mom_list ;

  int num_mom ;

  multi1d<LatticeComplex> phases ;

  Set sft_subsets ;

} ;

#endif
