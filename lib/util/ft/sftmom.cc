//  $Id: sftmom.cc,v 1.4 2003-04-01 02:46:43 edwards Exp $
//  $Log: sftmom.cc,v $
//  Revision 1.4  2003-04-01 02:46:43  edwards
//  Added const qual.
//
//  Revision 1.3  2003/03/20 19:34:25  flemingg
//  Evolved formfac_w.cc to use SftMom class, which included some bug fixes
//  in features in SftMom which had been previously untested and evolution
//  of the corresponding test program.
//
//  Revision 1.2  2003/03/14 17:13:44  flemingg
//  SftMom::sft() now works.
//
//  Revision 1.1  2003/03/14 05:06:06  flemingg
//  Initial version of SftMom class
//

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "proto.h"                 // part of QDP++, for crtesn()

using namespace QDP;


//! Function object used for constructing the time-slice set
class TimeSliceFunc : public SetFunc
{
public:
  TimeSliceFunc(int dir): dir_decay(dir) {}

  int operator() (const multi1d<int>& coordinate) const ;

  int numSubsets() const ;

private:
  TimeSliceFunc() {}  // hide default constructor

  int dir_decay;
};

int
TimeSliceFunc::operator() (const multi1d<int>& coordinate) const
{
  if ((dir_decay<0)||(dir_decay>=Nd)) {
    return 0 ;
  } else {
    return coordinate[dir_decay] ;
  }
}

int
TimeSliceFunc::numSubsets() const
{
  if ((dir_decay<0)||(dir_decay>=Nd)) {
    return 1 ;
  } else {
    return Layout::lattSize()[dir_decay] ;
  }
}


SftMom::SftMom(int mom2_max, bool avg_equiv_mom, int j_decay)
{
  multi1d<int> mom_offset ;

  if ((j_decay<0)||(j_decay>=Nd)) {
    mom_offset.resize(Nd) ;
  } else {
    mom_offset.resize(Nd-1) ;
  }
  mom_offset = 0 ;

  init(mom2_max, mom_offset, avg_equiv_mom, j_decay) ;
}

void
SftMom::init(int mom2_max, multi1d<int> mom_offset,
             bool avg_equiv_mom, int j_decay)
{
  // Averaging over equivalent momenta is only allowed if 
  // mom_offset is zero.
  if (avg_equiv_mom) {
    bool zero_offset = true ;

    for (int mu=0; mu < mom_offset.size(); ++mu) {
      if (mom_offset[mu] != 0) {
        zero_offset = false ;
      }
    }

    if (!zero_offset) {
      QDP_error_exit("SftMom: averaging not allowed with non-zero offset") ;
    }
  }

  sft_subsets.make(TimeSliceFunc(j_decay)) ;

  // determine the number of momenta with mom^2 <= (mom_max)^2
  // If avg_equiv_mom is true then only consider momenta with
  // mom[0] >= mom[1] >= ... >= mom[mu] >= ... >= 0
  multi1d<int> mom_size ;
  if ((j_decay<0)||(j_decay>=Nd)) {
    mom_size.resize(Nd) ;
  } else {
    mom_size.resize(Nd-1) ;
  }

  int L;
  int mom_vol = 1;

  for (L=1; L*L <= mom2_max; ++L) ;

  for(int mu=0; mu < mom_size.size(); ++mu) {
    if (avg_equiv_mom) {  
      mom_vol      *= L;
      mom_size[mu]  = L;
    } else {
      mom_vol      *= (2*L) + 1;
      mom_size[mu]  = (2*L) + 1;
    }
  }

  num_mom = 0;

  for(int n=0; n < mom_vol; ++n) {
    multi1d<int> mom = crtesn(n, mom_size);

    int mom2 = 0 ;

    for(int mu=0; mu < mom_size.size(); ++mu) {
      if (!avg_equiv_mom) mom[mu] -= L ;
      mom2 += mom[mu]*mom[mu];
    }

    if (mom2 > mom2_max) {
      continue;
    } else if (avg_equiv_mom) {
      // Ensure mom[0] >= mom[1] >= ... >= mom[mu] >= ... >= 0
      bool skip = false ;
      for(int mu=0; mu < mom_size.size()-1; ++mu)
        for(int nu=mu+1; nu < mom_size.size(); ++nu)
          if (mom[nu] > mom[mu]) skip=true;

      if (!skip) ++num_mom ;
    } else {
      ++num_mom ;
    }
  }

  // After all that shenanigans just to get num_mom, resize the momentum list
  mom_list.resize(num_mom, mom_size.size()) ;

  // Now we do exactly the same thing we did when counting, except this time
  // we can acutally fill the list
  int mom_num = 0;

  for(int n=0; n < mom_vol; ++n)
  {
    multi1d<int> mom = crtesn(n, mom_size);

    int mom2 = 0 ;

    for(int mu=0; mu < mom_size.size(); ++mu) {
      if (!avg_equiv_mom) mom[mu] -= L ;
      mom2 += mom[mu]*mom[mu];
    }

    if (mom2 > mom2_max) {
      continue;
    } else if (avg_equiv_mom) {
      // Ensure mom[0] >= mom[1] >= ... >= mom[mu] >= ... >= 0
      bool skip = false ;
      for(int mu=0; mu < mom_size.size()-1; ++mu)
        for(int nu=mu+1; nu < mom_size.size(); ++nu)
          if (mom[nu] > mom[mu]) skip = true ;

      if (!skip) mom_list[mom_num++] = mom ;
    } else {
      for (int mu=0; mu < mom_size.size(); ++mu) {
        mom_list[mom_num][mu] = mom_offset[mu] + mom[mu]  ;
      }
      ++mom_num ;
    }
  }

  // Now resize and initialize the Fourier phase table.  Then, loop over
  // allowed momenta, optionally averaging over equivalent momenta.
  phases.resize(num_mom) ;
  phases = 0. ;

  // Coordinates for sink momenta
  multi1d<LatticeInteger> my_coord(Nd);
  for (int mu=0; mu < Nd; ++mu)
    my_coord[mu] = Layout::latticeCoordinate(mu);

  // If averaging over equivalent momenta, we need redo mom_size and mom_vol
  // to allow both positive and negative momentum components
  multi1d<int> mom_degen ;

  if (avg_equiv_mom) {
    mom_vol = 1 ;

    for (int mu=0; mu < mom_size.size(); ++mu) {
      mom_vol      *= (2*L) + 1 ;
      mom_size[mu]  = (2*L) + 1 ;
    }

    // Keep track of |mom| degeneracy for averaging
    mom_degen.resize(num_mom) ;
    mom_degen = 0 ;
  }

  // reset mom_num
  mom_num = 0 ;

  for (int n=0; n < mom_vol; ++n) {
    multi1d<int> mom = crtesn(n, mom_size) ;

    int mom2 = 0 ;

    for(int mu=0; mu < mom_size.size(); ++mu) {
      mom[mu] -= L ;
      mom2 += mom[mu]*mom[mu];
    }

    // skip when (mom)^2 > (mom_max)^2
    if (mom2 > mom2_max) continue;

    // At this point, if (avg_equiv_mom == true) then we need to determine
    // mom_num by a fairly time consuming process.
    // If (avg_equiv_mom == false) then (mom == mom_list[mom_num])
    // (will double check this) and the momentum offset can be applied.

    if (avg_equiv_mom) {

      // determine mom_num for entering table mom_list
      // first step: make all the compontents positive
      multi1d<int> mom_tmp = mom ;
      for (int mu=0; mu < mom_tmp.size(); ++mu)
        if (mom_tmp[mu] < 0) mom_tmp[mu] = -mom_tmp[mu];

      // second step: sort the components
      // (using insertion sort, so we hope mom_tmp.size() <~ 10)

      // Initially, the first item is considered sorted.  mu divides mom
      // into a sorted region (<mu) and an unsorted one (>=mu)
      for (int mu=1; mu < mom_tmp.size(); ++mu) {
        // Select the item at the beginning of the unsorted region
        int v = mom_tmp[mu];
        // Work backwards, finding where v should go
        int nu = mu;
        // If this element is less than v, move it up one
        while (mom_tmp[nu-1] < v) {
          mom_tmp[nu] = mom_tmp[nu-1] ;
          --nu ;
          if (nu < 1) break ;
        }
        // Stopped when mom_tmp[nu-1] >= v, so put v at postion nu
        mom_tmp[nu] = v ;
      }

      // mom_tmp should now contain a momentum that appears in mom_list.
      // scan through list until we find a match.

      mom_num = -1 ;

      for(int k=0; k < num_mom; ++k) {
        bool match = true ;
        for (int mu=0; mu < mom_tmp.size(); ++mu) {
          if (mom_list[k][mu] != mom_tmp[mu]) {
            match = false ;
            break;
          }
        }
        if (match) {
          mom_num = k ;
          break ;
        }
      }

      if (mom_num < 0) {
        QDP_error_exit("SftMom: mom_num < 0. Shouldn't be here.\n") ;
      }

      // increment degeneracy for this mom_num
      ++(mom_degen[mom_num]) ;
    } else /* (avg_equiv_mom == false) */ {

      // apply momentum offset
      for (int mu=0; mu < mom_size.size(); ++mu) {
        mom[mu] += mom_offset[mu] ;
      }

      // double check that (mom == mom_list[n])
      // this check should never fail and could be removed in the future
      for (int mu=0; mu < mom_size.size(); ++mu) {
        if (mom[mu] != mom_list[mom_num][mu]) {
          // Should never get here !!!
          QDP_error_exit("SftMom: (mom != mom_list[mom_num])\n") ;
        }
      }
    } // end if (avg_equiv_mom)

    LatticeReal p_dot_x ;
    p_dot_x = 0. ;

    int j = 0;
    for(int mu = 0; mu < Nd; ++mu) {
      const Real twopi = 6.283185307179586476925286;

      if (mu == j_decay) continue ;

      p_dot_x += LatticeReal(my_coord[mu]) * twopi * Real(mom[j]) /
                   Layout::lattSize()[mu];
      ++j ;
    } // end for(mu)

    phases[mom_num] += cmplx(cos(p_dot_x), sin(p_dot_x)) ;

    // increment mom_num for next valid momenta
    ++mom_num ;

  } // end for (int n=0; n < mom_vol; ++n)

  // Finish averaging
  if (avg_equiv_mom) {
    for (int mom_num=0; mom_num < num_mom; ++mom_num)
      phases[mom_num] /= mom_degen[mom_num] ;
  }
}


// int
// SftMom::momToNum(const multi1d<int>& mom_in)
// {
//   multi1d<int> mom = mom_in ;
// 
//   // make all the compontents positive
//   for (int mu=0; mu < mom.size(); ++mu)
//     if (mom[mu] < 0) mom[mu] = -mom[mu];
// 
//   // sort the components (Insertion sort: hope mom.size() <~ 10)
//   // Initially, the first item is considered sorted.  mu divides mom
//   // into a sorted region (<mu) and an unsorted one (>=mu)
//   for (int mu=1; mu < mom.size(); ++mu) {
//     // Select the item at the beginning of the unsorted region
//     int v = mom[mu];
//     // Work backwards, finding where v should go
//     int nu = mu;
//     // If this element is less than v, move it up one
//     while (mom[nu-1] < v) {
//       mom[nu] = mom[nu-1] ;
//       --nu ;
//       if (nu < 1) break ;
//     }
//     // Stopped when mom[nu-1] >= v, so put v at postion nu
//     mom[nu] = v ;
//   }
// 
//   for(int mom_num=0; mom_num < num_mom; ++mom_num) {
//     bool match = true ;
//     for (int mu=0; mu < mom.size(); ++mu)
//       if (mom_list[mom_num][mu] != mom[mu]) {
//         match = false ;
//         break;
//       }
//     if (match) return mom_num ;
//   }
//   return -1;
// }

multi2d<DComplex>
SftMom::sft(const LatticeComplex& cf) const
{
  multi2d<DComplex> hsum(num_mom, sft_subsets.numSubsets()) ;

  for (int mom_num=0; mom_num < num_mom; ++mom_num)
    hsum[mom_num] = sumMulti(phases[mom_num]*cf, sft_subsets) ;

  return hsum ;
}
