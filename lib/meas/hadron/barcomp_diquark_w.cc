//  $Id: barcomp_diquark_w.cc,v 1.2 2007-02-28 03:28:40 edwards Exp $
/*! \file
 *  \brief Construct all components of a baryon propagator using a diquark
 */

#include "meas/hadron/diquark_w.h"
#include "meas/hadron/barcomp_diquark_w.h"

namespace Chroma 
{
 
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
   * \param diquark                  diquark ( Read )
   * \param quark_propagator_3       quark propagator ( Read )
   * \param spin_indices             holds list of source/sink spin indices ( Read )
   * \param phases                   object holds list of momenta ( Read )
   * \param t0                       coordinates of source in decay direction ( Read )
   * \param bc_spec                  boundary condition for spectroscopy ( Read )
   */

  void barcompDiquarkSparse(QQQSparse_t& barprop,
			    const QQDiquarkContract_t& diquark,
			    const LatticePropagator& quark_propagator_3,
			    const multi1d<QQQSpinIndices_t> spin_indices,
			    const SftMom& phases,
			    int t0, int bc_spec)
  {
    START_CODE();

    // Length of lattice in decay direction
    int length  = phases.numSubsets();

    // This object is composed of the spin indices and the correlators
    barprop.length = length;
    barprop.corrs.resize(spin_indices.size());

    // We need this fast, so at the expense of a lot of memory we will
    // expose all the color/spin indices of each propagator into a temporary
    multi2d< multi2d<LatticeComplex> > qc_3(unpackQuark(quark_propagator_3));

    // Temporaries
    multi1d<DComplex> hsum;
    LatticeComplex b_prop;
    multi1d<int> rnk(6);

    // Diquark-quark contract
    // Loop over all source/sink spin indices
    for(int i=0; i < spin_indices.size(); ++i)
    {
      const QQQSpinIndices_t& d = spin_indices[i];
      barprop.corrs[i].source_sink = d;

      // Here is where the odd-ball transpose is being done.
      // Note: source is first and sink is last. Don't blame RGE...
      rnk[0] = d.source[0];
      rnk[1] = d.source[1];
      rnk[2] = d.sink[0];
      rnk[3] = d.sink[1];

      // Contract over color indices - the final part of the antisym tensors
      b_prop = zero;
      for(rnk[4]=0; rnk[4] < Nc; ++rnk[4])         // color row
	for(rnk[5]=0; rnk[5] < Nc; ++rnk[5])       // color col
	{
	  b_prop += diquark.comp[rnk] * qc_3(d.source[2],d.sink[2])(rnk[4],rnk[5]);
	}

      /* Project on zero momentum: Do a slice-wise sum. */
      hsum = sumMulti(b_prop, phases.getSet());
      
      barprop.corrs[i].corr.resize(length);
      for(int t=0; t < length; ++t)   // t
      {
	int t_eff = (t - t0 + length) % length;
	barprop.corrs[i].corr[t] = (bc_spec < 0 && (t_eff+t0) >= length) ? -hsum[t] : hsum[t];
      }
    }

    END_CODE();
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
   * \param diquark                  diquark ( Read )
   * \param quark_propagator_3       quark propagator ( Read )
   * \param phases                   object holds list of momenta ( Read )
   * \param t0                       coordinates of source in decay direction ( Read )
   * \param bc_spec                  boundary condition for spectroscopy ( Read )
   */

  void barcompDiquarkDense(QQQDense_t& barprop,
			   const QQDiquarkContract_t& diquark,
			   const LatticePropagator& quark_propagator_3,
			   const SftMom& phases,
			   int t0, int bc_spec)
  {
    START_CODE();
 
    // Length of lattice in decay direction
    int length  = phases.numSubsets();
  
    multi1d<int> ranks(7);
    ranks = Ns;
    ranks[6] = length;

    barprop.length = length;
    barprop.corrs.resize(ranks);

    // We need this fast, so at the expense of a lot of memory we will
    // expose all the color/spin indices of each propagator into a temporary
    multi2d< multi2d<LatticeComplex> > qc_3(unpackQuark(quark_propagator_3));

    // Temporaries
    multi1d<DComplex> hsum;
    LatticeComplex b_prop;
    multi1d<int> rnk(6);

    // Diquark-quark contract
    for(ranks[0]=0; ranks[0] < Ns; ++ranks[0])           // sf_3
      for(ranks[1]=0; ranks[1] < Ns; ++ranks[1])         // sf_2
	for(ranks[2]=0; ranks[2] < Ns; ++ranks[2])       // sf_1
	  for(ranks[3]=0; ranks[3] < Ns; ++ranks[3])     // si_3
	    for(ranks[4]=0; ranks[4] < Ns; ++ranks[4])   // si_2
	      for(ranks[5]=0; ranks[5] < Ns; ++ranks[5]) // si_1
	      {
		// Here we carry on the odd-ball transpose.
		// Note: source is first and sink is last. This makes it
		// compatible with barcomp_w.cc . Don't blame RGE...
		rnk[0] = ranks[1];
		rnk[1] = ranks[2];
		rnk[2] = ranks[4];
		rnk[3] = ranks[5];

		// Contract over color indices - the final part of the antisym tensors
		b_prop = zero;
		for(rnk[4]=0; rnk[4] < Nc; ++rnk[4])         // color row
		  for(rnk[5]=0; rnk[5] < Nc; ++rnk[5])       // color col
		  {
		    b_prop += diquark.comp[rnk] * qc_3(ranks[3],ranks[0])(rnk[4],rnk[5]);
		  }

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

    END_CODE();
  }


}  // end namespace Chroma
