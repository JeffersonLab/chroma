//  $Id: mesons_w.cc,v 1.9 2003-03-14 17:16:13 flemingg Exp $
//  $Log: mesons_w.cc,v $
//  Revision 1.9  2003-03-14 17:16:13  flemingg
//  Variant 1 is now working with SftMom::sft().  In arbitrary units,
//  the relative performance seems to be: V1) 7.5, V2) 10, V3) 100.
//
//  Revision 1.8  2003/03/14 05:14:32  flemingg
//  rewrite of mesons_w.cc to use the new SftMom class.  mesons_w.cc still
//  needs to be cleaned up once the best strategy is resolved.  But for now,
//  the library and test program compiles and runs.
//
//  Revision 1.7  2003/03/06 03:38:35  edwards
//  Added start/end_code.
//
//  Revision 1.6  2003/03/06 02:07:12  flemingg
//  Changed the MomList class to eliminate an unneeded class member.
//
//  Revision 1.5  2003/03/06 00:30:14  flemingg
//  Complete rewrite of lib/meas/hadron/mesons_w.cc, including a simple test
//  program in mainprogs/tests built with 'make check' and various other
//  changes to autoconf/make files to support this rewrite.
//

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/mesons_w.h"
#include "proto.h"                 // part of QDP++, for crtesn()

using namespace QDP;

//! Meson 2-pt functions
/* This routine is specific to Wilson fermions!
 *
 * Construct meson propagators
 * The two propagators can be identical or different.
 *
 * \param quark_prop_1 -- first quark propagator ( Read )
 * \param quark_prop_2 -- second (anti-) quark propagator ( Read )
 * \param t0 -- timeslice coordinate of the source ( Read )
 * \param phases -- object holds list of momenta and Fourier phases ( Read )
 * \param nml -- namelist file object ( Read )
 * \param nml_group -- string used for writing nml data ( Read )
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
            SftMom& phases,
            int t0,
            NmlWriter& nml,
	    char* nml_group)
{
  START_CODE("mesons");

  // Length of lattice in decay direction
  int length = phases.numSubsets() ;

  // Construct the anti-quark propagator from quark_prop_2
  int G5 = Ns*Ns-1;
  LatticePropagator anti_quark_prop =  Gamma(G5) * quark_prop_2 * Gamma(G5);

  // GTF: We're going to consider several working variants here to see
  // which one works best.  Now that Variant 1 is working, it seems to be
  // 25% faster than Variant 2 in the test program.

#if 1    // Variant 1

  // This variant uses the function SftMom::sft() to do all the work
  // computing the Fourier transform of the meson correlation function
  // inside the class SftMom where all the of the Fourier phases and
  // momenta are stored.  It's primary disadvantage is that it
  // requires more memory because it does all of the Fourier transforms
  // at the same time.

  // Loop over gamma matrix insertions
  for (int gamma_value=0; gamma_value < (Ns*Ns); ++gamma_value) {

    // Construct the meson correlation function
    LatticeComplex corr_fn ;
    corr_fn = trace(adj(anti_quark_prop) * Gamma(gamma_value) *
                    quark_prop_1 * Gamma(gamma_value)) ;

    multi2d<DComplex> hsum ;
    hsum = phases.sft(corr_fn) ;

    for (int sink_mom_num=0; sink_mom_num < phases.numMom(); ++sink_mom_num) {
      multi1d<Double> mesprop(phases.numSubsets()) ;

      for (int t=0; t < phases.numSubsets(); ++t) {
        int t_eff = (t - t0 + length) % length ;

        mesprop[t_eff] = real(hsum[sink_mom_num][t]) ;
      }

      // Print out the results
      push(nml, nml_group) ;
      Write(nml, gamma_value) ;
      write(nml, "sink_mom", phases.numToMom(sink_mom_num)) ;
      Write(nml, mesprop) ;
      pop(nml) ;

    } // end for(sink_mom_num)

  } // end for(gamma_value)

#elif 1  // Variant 2

  // This variant does the Fourier transform in this code after getting
  // a copy of the Fourier phases from the SftMom object.  It should require
  // less memory than Variant 1 since the Fourier transform for each momenta
  // is done separately.  The outer loop is gamma matrix insertions, which
  // means that corr_fn is only computed once per gamma_value.  Since the
  // inner loop is sink_mom, this means that the implicit copy done by
  // phases[sink_mom_num] is redone for each gamma_value.  Whether this
  // is cheaper than Variant 3 must be tested.  According to my test code,
  // this variant is about 10 times faster than Variant 3.

  // Loop over gamma matrix insertions
  for (int gamma_value=0; gamma_value < (Ns*Ns); ++gamma_value) {

    // Construct the meson correlation function
    LatticeComplex corr_fn ;
    corr_fn = trace(adj(anti_quark_prop) * Gamma(gamma_value) *
                    quark_prop_1 * Gamma(gamma_value)) ;

    for (int sink_mom_num=0; sink_mom_num < phases.numMom(); ++sink_mom_num) {

      multi1d<Double> hsum ;

      hsum = sumMulti(real(phases[sink_mom_num]*corr_fn), phases.getSubset()) ;

      // This is obviously a copy with simple relabeling which might
      // be eliminated in the future.  It is done to match what was
      // done in the past with NML.
      multi1d<Double> mesprop(hsum.size()) ;

      for (int t=0; t < phases.numSubsets(); ++t) {
        int t_eff = (t - t0 + length) % length ;
        mesprop[t_eff] = hsum[t] ;
      }

      // Print out the results
      push(nml, nml_group) ;
      Write(nml, gamma_value) ;
      write(nml, "sink_mom", phases.numToMom(sink_mom_num)) ;
      Write(nml, mesprop) ;
      pop(nml) ;

    } // end for(sink_mom_num)

  } // end for(gamma_value)

#elif 1  // Variant 3

  // This variant is like Variant 2 with the inner and outer loops swapped.
  // This means that the implicit copy phases[sink_mom_num] is only done
  // as few times as necessary but now the computation of the meson
  // correlation function is repeated SftMom::numMom() times.

  // Loop over sink momenta first
  for (int sink_mom_num=0; sink_mom_num < phases.numMom(); ++sink_mom_num) {

    LatticeComplex phasefac = phases[sink_mom_num] ;

    // Loop over gamma matrix insertions
    for (int gamma_value=0; gamma_value < (Ns*Ns); ++gamma_value) {

      // Construct the meson correlation function
      LatticeComplex corr_fn ;
      corr_fn = trace(adj(anti_quark_prop) * Gamma(gamma_value) *
                      quark_prop_1 * Gamma(gamma_value)) ;

      multi1d<Double> hsum ;

      hsum = sumMulti(real(phasefac*corr_fn), phases.getSubset()) ;

      // This is obviously a copy with simple relabeling which might
      // be eliminated in the future.  It is done to match what was
      // done in the past with NML.
      multi1d<Double> mesprop(hsum.size()) ;

      for (int t=0; t < phases.numSubsets(); ++t) {
        int t_eff = (t - t0 + length) % length ;
        mesprop[t_eff] = hsum[t] ;
      }

      // Print out the results
      push(nml, nml_group) ;
      Write(nml, gamma_value) ;
      write(nml, "sink_mom", phases.numToMom(sink_mom_num)) ;
      Write(nml, mesprop) ;
      pop(nml) ;

    } // end for(gamma_value)

  } // end for(sink_mom_num)

#endif

  END_CODE("mesons");
}
