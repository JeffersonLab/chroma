//  $Id: barcomp_w.cc,v 3.3 2009-03-19 17:17:20 mcneile Exp $
/*! \file
 *  \brief Construct all components of a baryon propagator
 */


#include "meas/hadron/barcomp_w.h"

namespace Chroma 
{
 
  //! Serialize generalized object
  multi1d<ComplexF> QQQSparse_t::serialize()
  { 
    // Dreadful hack - we will use 3 complexs to hold 6 ints. So, the size
    // of the array is the size of the correlators plus the size of the complexs
    multi1d<ComplexF> barprop_1d(corrs.size()*(3 + length));

    int cnt = 0;
    for(int s=0; s < corrs.size(); ++s)
    {
      const QQQElem_t& elem = corrs[s];
      const QQQSpinIndices_t& sp = elem.source_sink;

      // Dreadful hack - use a complex to hold an int
      barprop_1d[cnt++] = cmplx(Real(sp.source[0]), Real(sp.source[1]));
      barprop_1d[cnt++] = cmplx(Real(sp.source[2]), Real(sp.sink[0]));
      barprop_1d[cnt++] = cmplx(Real(sp.sink[1]), Real(sp.sink[2]));

      for(int t=0; t < length; ++t)
      {
	barprop_1d[cnt++] = elem.corr[t];
      }
    }

    if (cnt != barprop_1d.size())
    {
      QDPIO::cerr << "Size mismatch in QQQSparse_t serialization" << endl;
      QDP_abort(1);
    }

    return barprop_1d;
  }



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
   * \param spin_indices             holds list of source/sink spin indices ( Read )
   * \param phases                   object holds list of momenta ( Read )
   * \param t0                       coordinates of source in decay direction ( Read )
   * \param bc_spec                  boundary condition for spectroscopy ( Read )
   */

  void barcompSparse(QQQSparse_t& barprop,
		     const LatticePropagator& quark_propagator_1, 
		     const LatticePropagator& quark_propagator_2,
		     const LatticePropagator& quark_propagator_3,
		     const multi1d<QQQSpinIndices_t> spin_indices,
		     const SftMom& phases,
		     int t0, int bc_spec)
  {
    START_CODE();
    if ( Nc != 3 ){    /* Code is specific to Ns=4 and Nc=3. */
      QDPIO::cerr<<" code only works for Nc=3 and Ns=4\n";
      QDP_abort(111) ;
    }
#if QDP_NC == 3


    // Length of lattice in decay direction
    int length  = phases.numSubsets();

    // This object is composed of the spin indices and the correlators
    barprop.length = length;
    barprop.corrs.resize(spin_indices.size());

    // Temporaries
    multi1d<DComplex> hsum;
    LatticeComplex b_prop;

    // Loop over all source/sink spin indices
    for(int i=0; i < spin_indices.size(); ++i)
    {
      const QQQSpinIndices_t& d = spin_indices[i];
      barprop.corrs[i].source_sink = d;

      // Contract over color indices with antisym tensors
      b_prop = 
	colorContract(peekSpin(quark_propagator_1,
			       d.source[0],d.sink[0]),  // (si_1,sf_1)
		      peekSpin(quark_propagator_2,
			       d.source[1],d.sink[1]),  // (si_2,sf_2)
		      peekSpin(quark_propagator_3,
			       d.source[2],d.sink[2])); // (si_3,sf_3)

      /* Project on zero momentum: Do a slice-wise sum. */
      hsum = sumMulti(b_prop, phases.getSet());
      
      barprop.corrs[i].corr.resize(length);
      for(int t=0; t < length; ++t)   // t
      {
	int t_eff = (t - t0 + length) % length;
	barprop.corrs[i].corr[t] = (bc_spec < 0 && (t_eff+t0) >= length) ? -hsum[t] : hsum[t];
      }
    }

#endif
    END_CODE();
  }



  //! Serialize generalized object
  multi1d<ComplexF> QQQDense_t::serialize()
  {
    multi1d<int> ranks(7);

    multi1d<ComplexF> barprop_1d(length*Ns*Ns*Ns*Ns*Ns*Ns);

    int cnt = 0;
    for(ranks[0]=0; ranks[0] < Ns; ++ranks[0])           // sf_3
      for(ranks[1]=0; ranks[1] < Ns; ++ranks[1])         // sf_2
	for(ranks[2]=0; ranks[2] < Ns; ++ranks[2])       // sf_1
	  for(ranks[3]=0; ranks[3] < Ns; ++ranks[3])     // si_3
	    for(ranks[4]=0; ranks[4] < Ns; ++ranks[4])   // si_2
	      for(ranks[5]=0; ranks[5] < Ns; ++ranks[5]) // si_1
		for(ranks[6] = 0; ranks[6] < length; ++ranks[6])
		{
		  barprop_1d[cnt++] = corrs[ranks];
		}

    return barprop_1d;
  }



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

  void barcomp(QQQDense_t& barprop,
	       const LatticePropagator& quark_propagator_1, 
	       const LatticePropagator& quark_propagator_2,
	       const LatticePropagator& quark_propagator_3,
	       const SftMom& phases,
	       int t0, int bc_spec)
  {
    START_CODE();
    if ( Nc != 3 ){    /* Code is specific to Ns=4 and Nc=3. */
      QDPIO::cerr<<"barcomp code only works for Nc=3 and Ns=4\n";
      QDP_abort(111) ;
    }
#if QDP_NC == 3


    // Length of lattice in decay direction
    int length  = phases.numSubsets();
  
    multi1d<int> ranks(7);
    ranks = Ns;
    ranks[6] = length;

    barprop.length = length;
    barprop.corrs.resize(ranks);

    // Temporaries
    multi1d<DComplex> hsum;
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

		  barprop.corrs[ranks] = 
		    (bc_spec < 0 && (t_eff+t0) >= length) ? -hsum[ranks[6]] : 
		    hsum[ranks[6]];
		}
	      }

#endif
    END_CODE();
  }


}  // end namespace Chroma

