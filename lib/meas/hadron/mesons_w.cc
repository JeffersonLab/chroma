//  $Id: mesons_w.cc,v 1.5 2003-03-06 00:30:14 flemingg Exp $
//  $Log: mesons_w.cc,v $
//  Revision 1.5  2003-03-06 00:30:14  flemingg
//  Complete rewrite of lib/meas/hadron/mesons_w.cc, including a simple test
//  program in mainprogs/tests built with 'make check' and various other
//  changes to autoconf/make files to support this rewrite.
//

#include "chromabase.h"
#include "meas/hadron/mesons_w.h"
#include "proto.h"                 // part of QDP++, for crtesn()

using namespace QDP;

//! Function object used for constructing the time-slice set
class TimeSliceFunc : public SetFunc
{
public:
  TimeSliceFunc(int dir): dir_decay(dir) {}

  int operator() (const multi1d<int>& coordinate) const {return coordinate[dir_decay];}
  int numSubsets() const {return Layout::lattSize()[dir_decay];}

  int dir_decay;

private:
  TimeSliceFunc() {}  // hide default constructor
};


class MomList
{
public:
  MomList(int mom2_max) ;

  int momToNum(const multi1d<int>& mom);

  int numMom() {return num_mom ;}

  multi1d<int>& numToMom(int num);

  ~MomList() ;

private:
  MomList() {} // hide default constructor

  multi2d<int> *mom_list ;

  multi1d<int> *mom ;

  int num_mom;
};

 
//! Meson 2-pt functions
/* This routine is specific to Wilson fermions!
 *
 * Construct meson propagators
 * The two propagators can be identical or different.
 *
 * \param quark_prop_1 -- first quark propagator ( Read )
 * \param quark_prop_2 -- second (anti-) quark propagator ( Read )
 * \param t_source -- cartesian coordinates of the source ( Read )
 * \param sink_mom2_max -- max sink hadron mom squared ( Read )
 * \param j_decay -- direction of the exponential decay ( Read )
 * \param nml -- namelist file object ( Read )
 *
 *        ____
 *        \
 * m(t) =  >  < m(t_source, 0) m(t + t_source, x) >
 *        /
 *        ----
 *          x
 */

void mesons(const LatticePropagator& quark_prop_1,
            const LatticePropagator& quark_prop_2, 
            const multi1d<int>& t_source,
            int sink_mom2_max,
            int j_decay,
            NmlWriter& nml)
{
  // Create the time-slice set
  Set timeslice;
  timeslice.make(TimeSliceFunc(j_decay));

  // Length of lattice in j_decay direction
  int length = timeslice.numSubsets();

  // Create the MomList object
  MomList list(sink_mom2_max) ;

  // Meson correlations function.
  LatticeComplex corr_fn;

  // Create 2-D array for the meson correlation function on each timeslice
  // for each unique momentum.
  multi2d<Double> hsum(list.numMom(), length) ;

  // Coordinates for sink momenta
  multi1d<LatticeInteger> my_coord(Nd);
  for (int mu=0; mu < Nd; ++mu)
    my_coord[mu] = Layout::latticeCoordinate(mu);

  // Construct the anti-quark propagator from quark_prop_2
  int t0 = t_source[j_decay];
  int G5 = Ns*Ns-1;
  LatticePropagator anti_quark_prop =  Gamma(G5) * quark_prop_2 * Gamma(G5);

  // Loop over gamma matrix insertions
  for (int gamma_value=0; gamma_value < (Ns*Ns); ++gamma_value) {

    // Construct the meson correlation function
    corr_fn = trace(adj(anti_quark_prop) * Gamma(gamma_value) * quark_prop_1 *
                    Gamma(gamma_value));

    // Loop over allowed sink momenta: (sink_mom)^2 <= sink_mom2_max.
    // Do this by constructing a L^(Nd-1) grid in momenta centered about the
    // origin. Loop lexicographically over all the "sites" (momenta value)
    // and toss out ones that are too large.
    // 
    // NOTE: spatial anisotropy is no allowed here
    int Ndm1 = Nd-1 ;
    multi1d<int> sink_mom_size(Ndm1);
    int L ;
    int sink_mom_vol = 1 ;
  
    for (L=0; (L+1)*(L+1) <= sink_mom2_max; ++L) ;

    for (int mu=0; mu < Ndm1; ++mu) {
      sink_mom_vol *= (2*L)+1 ;
      sink_mom_size[mu] = (2*L)+1 ;
    }

    // Keep track of |sink_mom| degeneracy for averaging
    multi1d<int> sink_mom_degen(list.numMom()) ;

    sink_mom_degen = 0 ;
    hsum = 0. ;

    for (int n=0; n < sink_mom_vol; ++n) {
      multi1d<int> sink_mom = crtesn(n, sink_mom_size) ;

      int sink_mom2 = 0 ;
      for (int mu=0; mu<Ndm1; ++mu) {
        sink_mom[mu] -= L ;
        sink_mom2 += sink_mom[mu] * sink_mom[mu] ;
      }

      // skip when (sink_mom)^2 > (sink_mom_max)^2
      if (sink_mom2 > sink_mom2_max) continue ;

      int sink_mom_num = list.momToNum(sink_mom) ;
      sink_mom_degen[sink_mom_num]++ ;

      LatticeReal p_dot_x(float(0.0)) ;

      int j = 0;
      for(int mu = 0; mu < Nd; ++mu) {
        const Real twopi = 6.283185307179586476925286;

        if (mu == j_decay) continue ;

        p_dot_x += LatticeReal(my_coord[mu]) * twopi * Real(sink_mom[j]) /
                     Layout::lattSize()[mu];
        j++;
      } // end for(mu)

      // The complex phases exp(i p.x)
      LatticeComplex phasefac = cmplx(cos(p_dot_x), sin(p_dot_x));

      hsum[sink_mom_num] += sumMulti(real(phasefac*corr_fn), timeslice) ;

    } // end for(n)

    for (int sink_mom_num=0; sink_mom_num < list.numMom(); ++sink_mom_num) {
      multi1d<Double> mesprop(length) ;

      for (int t=0; t < length; ++t) {
        int t_eff = (t - t0 + length) % length ;

        // Finish averaging
        mesprop[t_eff] = hsum[sink_mom_num][t] /
                           Double(sink_mom_degen[sink_mom_num]) ;
      }

      // Print out the results
      push(nml, "Wilson_Mesons") ;
      Write(nml, gamma_value) ;
      Write(nml, j_decay) ;
      write(nml, "sink_mom", list.numToMom(sink_mom_num)) ;
      Write(nml, mesprop) ;
      pop(nml) ;
    }

  } // end for(gamma_value)
}


MomList::MomList(int mom2_max)
{
  int Ndm1 = Nd - 1;
  mom = new multi1d<int>(Ndm1) ;

  // determine the number of unique momenta with mom^2 <= (mom_max)^2
  multi1d<int> mom_size(Ndm1);
  int L;
  int mom_vol = 1;

  for (L=1; L*L <= mom2_max; ++L) ;

  for(int mu=0; mu<Ndm1; ++mu)
  {
    mom_vol *= L;
    mom_size[mu] = L;
  }

  num_mom = 0;

  for(int n=0; n < mom_vol; ++n)
  {
    *mom = crtesn(n, mom_size);

    int mom2 = 0 ;

    for(int mu=0; mu<Ndm1; ++mu)
      mom2 += (*mom)[mu]*(*mom)[mu];

    if (mom2 > mom2_max) {
      continue;
    } else {
      // Ensure mom[0] >= mom[1] >= ... >= mom[Ndm1]
      int skip=0;
      for(int mu=0; mu<Ndm1-1; ++mu)
        for(int nu=mu+1; nu<Ndm1; ++nu)
          if ((*mom)[nu] > (*mom)[mu]) skip=1;

      if (!skip) num_mom++ ;
    }
  }

  // After all that shenanigans just to get num_mom, create the momentum list
  mom_list = new multi2d<int>(num_mom,Ndm1) ;

  // Now we do exactly the same thing we did when counting, except this time
  // we can acutally fill the list
  int mom_num = 0;

  for(int n=0; n < mom_vol; ++n)
  {
    *mom = crtesn(n, mom_size);

    int mom2 = 0 ;

    for(int mu=0; mu<Ndm1; ++mu)
      mom2 += (*mom)[mu]*(*mom)[mu];

    if (mom2 > mom2_max) {
      continue;
    } else {
      // Ensure mom[0] >= mom[1] >= ... >= mom[Ndm1]
      int skip=0;
      for(int mu=0; mu<Ndm1-1; ++mu)
        for(int nu=mu+1; nu<Ndm1; ++nu)
          if ((*mom)[nu] > (*mom)[mu]) skip=1;

      if (!skip) (*mom_list)[mom_num++] = *mom;
    }
  }
}


int
MomList::momToNum(const multi1d<int>& mom_in)
{
  int Ndm1 = Nd-1;
  *mom = mom_in ;

  // make all the compontents positive
  for (int mu=0; mu<Ndm1; ++mu)
    if ((*mom)[mu] < 0) (*mom)[mu] = -(*mom)[mu];

  // sort the components (Insertion sort: hope Ndm1 <~ 10)
  // Initially, the first item is considered sorted.  mu divides (*mom)
  // into a sorted region (<mu) and an unsorted one (>=mu)
  for (int mu=1; mu<Ndm1; ++mu) {
    // Select the item at the beginning of the unsorted region
    int v = (*mom)[mu];
    // Work backwards, finding where v should go
    int nu = mu;
    // If this element is less than v, move it up one
    while ((*mom)[nu-1] < v) {
      (*mom)[nu] = (*mom)[nu-1] ;
      --nu ;
      if (nu < 1) break ;
    }
    // Stopped when (*mom)[nu-1] >= v, so put v at postion nu
    (*mom)[nu] = v ;
  }

  for(int mom_num=0; mom_num<num_mom; ++mom_num) {
    int match=1;
    for (int mu=0; mu<Ndm1; ++mu)
      if ((*mom_list)[mom_num][mu] != (*mom)[mu]) {
        match=0;
        break;
      }
    if (match) return mom_num ;
  }
  return -1;
}


multi1d<int>&
MomList::numToMom(int mom_num)
{
  *mom = (*mom_list)[mom_num] ;
  return *mom ;
}


MomList::~MomList()
{
  delete mom_list ;
  delete mom ;
}
