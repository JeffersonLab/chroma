//  $Id: mescomp_w.cc,v 3.0 2006-04-03 04:59:00 edwards Exp $
/*! \file
 *  \brief Construct all components of a meson propagator
 */

#include "chromabase.h"
#include "meas/hadron/mescomp_w.h"

namespace Chroma 
{
 
  //! Convert generalized correlator object
  void convertMescomp(multi1d<Complex>& mesprop_1d, const multiNd<Complex>& mesprop, 
		      const int j_decay)
  {
    int length = Layout::lattSize()[j_decay]; // Temporal extent of lattice

    multi1d<int> ranks(5);

    mesprop_1d.resize(length*Ns*Ns*Ns*Ns);

    int cnt = 0;
    for(ranks[0]=0; ranks[0] < Ns; ++ranks[0])           // sf_2
      for(ranks[1]=0; ranks[1] < Ns; ++ranks[1])         // sf_1
	for(ranks[2]=0; ranks[2] < Ns; ++ranks[2])       // si_2
	  for(ranks[3]=0; ranks[3] < Ns; ++ranks[3])     // si_1
	    for(ranks[4] = 0; ranks[4] < length; ++ranks[4])
	    {
	      mesprop_1d[cnt++] = mesprop[ranks];
	    }
  }


  //! Construct all components of a meson propagator
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions!
   *
   * In all mesons the colour components are contracted leaving only the
   * spin components.
   *
   * \param mesprop                  meson correlation function (in real space) ( Write )
   * \param quark_propagator_1       quark propagator ( Read )
   * \param quark_propagator_2       quark propagator ( Read )
   * \param phases                   object holds list of momenta ( Read )
   * \param t0                       coordinates of source in decay direction ( Read )
   */

  void mescomp(multiNd<Complex>& mesprop,
	       const LatticePropagator& quark_propagator_1, 
	       const LatticePropagator& quark_propagator_2,
	       const SftMom& phases,
	       int t0)
  {
    START_CODE();

    // Length of lattice in decay direction
    int length  = phases.numSubsets();
    int j_decay = phases.getDir();
    int G5      = Ns*Ns-1;

    multi1d<DComplex> hsum;
  
    multi1d<int> ranks(5);
    ranks = Ns;
    ranks[4] = length;
    mesprop.resize(ranks);

    LatticePropagator antiquark_1 = adj(Gamma(G5) * quark_propagator_1 * Gamma(G5));
    LatticeComplex m_prop;

    for(ranks[0]=0; ranks[0] < Ns; ++ranks[0])           // sf_2
      for(ranks[1]=0; ranks[1] < Ns; ++ranks[1])         // sf_1
	for(ranks[2]=0; ranks[2] < Ns; ++ranks[2])       // si_2
	  for(ranks[3]=0; ranks[3] < Ns; ++ranks[3])     // si_1
	  {
	    // Contract over color indices with antisym tensors
	    m_prop = 
	      traceColor(peekSpin(antiquark_1,
				  ranks[1],ranks[3])   // (sf_1,si_1)
		       * peekSpin(quark_propagator_2,
				  ranks[0],ranks[2])); // (sf_2,si_2)

	    /* Project on zero momentum: Do a slice-wise sum. */
	    hsum = sumMulti(m_prop, phases.getSet());
      
	    for(ranks[4] = 0; ranks[4] < length; ++ranks[4])   // t
	    {
	      int t_eff = (ranks[4] - t0 + length) % length;
	      mesprop[ranks] = hsum[t_eff];
	    }
	  }

    END_CODE();
  }

}  // end namespace Chroma
