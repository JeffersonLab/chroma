//  $Id: barcomp_w.cc,v 1.12 2005-01-14 18:42:35 edwards Exp $
/*! \file
 *  \brief Construct all components of a baryon propagator
 */

#include "chromabase.h"
#include "meas/hadron/barcomp_w.h"

namespace Chroma {
 
//! Construct all components of a baryon propagator
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions!
 *
 * In all baryons the colour components are contracted with the totally
 * antisymmetric 'tensor' eps(a,b,c) = antisym_tensor(a,b,c).
 *
 * \param barprop                  baryon correlation function (in real space) ( Write )
 * \param quark_propagator_1       quark propagator ( Read )
 * \param quark_propagator_2       quark propagator ( Read )
 * \param quark_propagator_3       quark propagator ( Read )
 * \param phases                   object holds list of momenta ( Read )
 * \param t0                       coordinates of source in decay direction ( Read )
 * \param bc_spec                  boundary condition for spectroscopy ( Read )
 */

void barcomp(multiNd<Complex>& barprop,
	     const LatticePropagator& quark_propagator_1, 
	     const LatticePropagator& quark_propagator_2,
	     const LatticePropagator& quark_propagator_3,
	     const SftMom& phases,
	     int t0, int bc_spec)
{
  START_CODE();

  // Length of lattice in decay direction
  int length  = phases.numSubsets();
  int j_decay = phases.getDir();

  multi1d<DComplex> hsum;
  
  multi1d<int> ranks(7);
  ranks = Ns;
  ranks[6] = length;
  barprop.resize(ranks);

  /* Initialize to zero */
  LatticeComplex b_prop;

  for(ranks[0]=0; ranks[0] < Ns; ++ranks[0])           // sf_3
    for(ranks[1]=0; ranks[1] < Ns; ++ranks[1])         // sf_2
      for(ranks[2]=0; ranks[2] < Ns; ++ranks[2])       // sf_1
	for(ranks[3]=0; ranks[3] < Ns; ++ranks[3])     // si_3
	  for(ranks[4]=0; ranks[4] < Ns; ++ranks[4])   // si_2
	    for(ranks[5]=0; ranks[5] < Ns; ++ranks[5]) // si_1
	    {
	      // Contract over color indices with antisym tensors
	      b_prop = 
		colorContract(peekSpin(quark_propagator_1,
				       ranks[5],ranks[2]),  // (si_1,sf_1)
			      peekSpin(quark_propagator_2,
				       ranks[4],ranks[1]),  // (si_2,sf_2)
			      peekSpin(quark_propagator_3,
				       ranks[3],ranks[0])); // (si_3,sf_3)

	      /* Project on zero momentum: Do a slice-wise sum. */
	      hsum = sumMulti(b_prop, phases.getSet());
      
	      for(ranks[6] = 0; ranks[6] < length; ++ranks[6])   // t
	      {
		int t_eff = (ranks[6] - t0 + length) % length;

		barprop[ranks] = 
		  (bc_spec < 0 && (t_eff+t0) >= length) ? -hsum[ranks[6]] : 
		  hsum[ranks[6]];
	      }
	    }

  END_CODE();
}

}  // end namespace Chroma
