//  $Id: barcomp_w.cc,v 1.8 2003-12-16 03:48:42 edwards Exp $
/*! \file
 *  \brief Construct all components of a baryon propagator
 */

#include "chromabase.h"
#include "io/qprop_io.h"
#include "meas/hadron/barcomp_w.h"

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

 
//! Construct all components of a baryon propagator
/*!
 * \ingroup hadron
 *
 * This routine is specific to Wilson fermions!
 *
 * quark_propagator_{1,2,3} -- quark propagators (read)
 * header_{1,2,3} -- second quark propagator ( Read )
 * quark_propagator_3 -- third quark propagator  ( Read )
 * barprop -- baryon correlation function (in real space) ( Write )
 * num_mom -- number of non-zero momenta ( Read )
 * t0      -- cartesian coordinates of the source in the j_decay direction( Read )
 * j_decay -- direction of the exponential decay ( Read )
 * bc_spec  -- boundary condition for spectroscopy ( Read )
 * file - the filename of the binary output file
 * In all baryons the colour components are contracted with the totally
 * antisymmetric 'tensor' eps(a,b,c) = antisym_tensor(a,b,c).
 */

void barcomp(const LatticePropagator& quark_propagator_1, 
	     const PropHead& head_1,
	     const LatticePropagator& quark_propagator_2,
	     const PropHead& head_2,
	     const LatticePropagator& quark_propagator_3,
	     const PropHead& head_3,
	     int t0, int j_decay, int bc_spec,
	     const char file[],
	     NmlWriter& nml)
{
  // Create the time-slice set
  UnorderedSet timeslice;
  timeslice.make(TimeSliceFunc(j_decay));

  // Length of lattice in j_decay direction
  int length = timeslice.numSubsets();

  multi1d<DComplex> hsum(length);
  
  push(nml, "head_1"); writePropHeadinNml(head_1, nml); pop(nml);
  push(nml, "head_2"); writePropHeadinNml(head_2, nml); pop(nml);
  push(nml, "head_3"); writePropHeadinNml(head_3, nml);

  multi1d<int> ranks(7);
  ranks = Ns;
  ranks[6] = length;
  multiNd<Complex> barprop(ranks);

  /*  START_CODE("barcomp");*/

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
	      hsum = sumMulti(b_prop, timeslice);
      
	      for(ranks[6] = 0; ranks[6] < length; ++ranks[6])   // t
	      {
		int t_eff = (ranks[6] - t0 + length) % length;

		barprop[ranks] = 
		  (bc_spec < 0 && (t_eff+t0) >= length) ? -hsum[ranks[6]] : 
		  hsum[ranks[6]];
	      }
	    }

  /*
   *  Now write out the file
   */

  writeBarcomp(file, barprop, head_1, head_2, head_3,j_decay);

  /*  END_CODE("barcomp");*/
}
