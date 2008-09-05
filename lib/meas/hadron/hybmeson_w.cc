// $Id: hybmeson_w.cc,v 3.2 2008-09-05 13:59:18 edwards Exp $
/*! \file
 *  \brief Hybrid meson 2-pt functions
 */

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/hybmeson_w.h"

namespace Chroma 
{

  //! Print the correlator to xml
  static void print_disp(XMLWriter& xml_hyb, const LatticeComplex& corr_fn,
			 const SftMom& phases, int t0)
  {
    int length  = phases.numSubsets();

    multi2d<DComplex> hsum;
    hsum = phases.sft(corr_fn);

    // Loop over sink momenta
    XMLArrayWriter xml_sink_mom(xml_hyb,phases.numMom());
    push(xml_sink_mom, "momenta");

    for (int sink_mom_num=0; sink_mom_num < phases.numMom(); ++sink_mom_num) 
    {
      push(xml_sink_mom);
      write(xml_sink_mom, "sink_mom_num", sink_mom_num);
      write(xml_sink_mom, "sink_mom", phases.numToMom(sink_mom_num));

      multi1d<Complex> mesprop(length);
      for (int t=0; t < length; ++t) 
      {
	int t_eff = (t - t0 + length) % length;
	mesprop[t_eff] = real(hsum[sink_mom_num][t]);
      }

      write(xml_sink_mom, "mesprop", mesprop);
      pop(xml_sink_mom);
    } // end for(sink_mom_num)
 
    pop(xml_sink_mom);
  }



  //! Hybrid meson 2-pt functions
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions!
   *
   * First we construct a hybrid pion and 3 hybrid rho's, followed by
   * an exotic 0^{+-}, an exotic 0^{--} and finally 2*3 exotic 1^{-+}'s.
   *
   * \param f             field strength tensor ( Read )
   * \param u_smr         the SMEARED gauge field, used in constructing the f's
   * \param quark_prop_1  first quark propagator ( Read )
   * \param quark_prop_2  second (anti-) quark propagator ( Read )
   * \param t_source      cartesian coordinates of the source ( Read )
   * \param phases        object holds list of momenta and Fourier phases ( Read )
   * \param xml           xml file object ( Read )
   * \param xml_group     string used for writing xml data ( Read )
   *
   *        ____
   *        \
   * m(t) =  >  < m(t_source, 0) m(t + t_source, x) >
   *        /                    
   *        ----
   *          x 
   */

  void hybmeson(const multi1d<LatticeColorMatrix>& f, 
		const multi1d<LatticeColorMatrix>& u_smr, 
		const LatticePropagator& quark_prop_1,
		const LatticePropagator& quark_prop_2, 
		const SftMom& phases,
		multi1d<int> t_source,
		XMLWriter& xml,
		const string& xml_group)
  {
    START_CODE();

    // Group
    XMLArrayWriter xml_hyb(xml,15);
    push(xml_hyb, xml_group);

    if ( Nd != 4 )		/* Code is specific to Nd=4. */
    {
      END_CODE();
      return;
    }

    // Length of lattice in decay direction
    int length  = phases.numSubsets();
    int j_decay = phases.getDir();
    int t0 = t_source[j_decay];

    // Fill f_source everywhere with f(t_source)
    multi1d<LatticeColorMatrix> f_source(Nd*(Nd-1)/2);
    for(int m = 0; m < Nd*(Nd-1)/2; ++m)
      f_source[m] = peekSite(f[m], t_source);


    // Construct the anti-quark propagator from quark_prop_2
    int G5 = Ns*Ns-1;
    LatticePropagator anti_quark_prop =  Gamma(G5) * quark_prop_2 * Gamma(G5);

    /* Cyclic direction index */
    multi1d<int> kp1(Nd);
    multi2d<int> n_munu(Nd, Nd);
    for(int k = 0; k < Nd; ++k)
    {
      if( k != j_decay )
      {
	int n = k+1;
	if ( n == j_decay ) n++;
	if ( n == Nd ) n = 0;
	kp1[k] = n;
      }
      else
	kp1[k] = -99999;
    }

    /* Index to F_{mu,nu} */
    int icnt = 0;
    for(int m=0; m < Nd-1; ++m)
      for(int n=m+1; n < Nd; ++n)
      {
	n_munu[n][m] = icnt;
	n_munu[m][n] = icnt;
	icnt++;
      }

    // Hybrid pion: \epsilon_{k,j,n} \bar \psi \gamma_k F_{j,n} \psi
    int kv = 0;
    {
      push(xml_hyb);     // next array element
      write(xml_hyb, "kv", kv);

      LatticePropagator q1_prop = 0;
      LatticePropagator q2_prop = 0;
      for(int k = 0; k < Nd; ++k)
	if( k != j_decay )
	{
	  int jm = 1 << k;
	  int j = kp1[k];
	  int n = kp1[j];
	  int m = n_munu[n][j];

	  if( j < n )
	  {
	    q1_prop += f[m] * (Gamma(jm) * quark_prop_1);
	    q2_prop += (anti_quark_prop * Gamma(jm)) * f_source[m];
	  }
	  else
	  {
	    q1_prop -= f[m] * (Gamma(jm) * quark_prop_1);
	    q2_prop -= (anti_quark_prop * Gamma(jm)) * f_source[m];
	  }
	}

      LatticeComplex corr_fn = localInnerProduct(q2_prop, q1_prop);
      print_disp(xml_hyb, corr_fn, phases, t0);

      pop(xml_hyb);
    }

    // Hybrid rho: \epsilon_{k,j,n} \bar \psi \gamma_5 F_{j,n} \psi
    kv = 0;
    for(int k = 0; k < Nd; ++k)
      if( k != j_decay )
      {
	kv++;
      
	push(xml_hyb);     // next array element
	write(xml_hyb, "kv", kv);

	int j = kp1[k];
	int n = kp1[j];
	int m = n_munu[n][j];

	LatticeComplex corr_fn = localInnerProduct(anti_quark_prop * f_source[m], 
						   f[m] * (Gamma(G5) * quark_prop_1 * Gamma(G5)));
	print_disp(xml_hyb, corr_fn, phases, t0);

	pop(xml_hyb);
      }

    // Exotic 0^{+-}: \epsilon_{k,j,n} \bar \psi \gamma_5 \gamma_k F_{j,n} \psi
    kv = 4;
    {
      push(xml_hyb);     // next array element
      write(xml_hyb, "kv", kv);

      LatticePropagator q1_prop = 0;
      LatticePropagator q2_prop = 0;
      for(int k = 0; k < Nd; ++k)
	if( k != j_decay )
	{
	  int m = 1 << k;
	  int jm = m ^ G5;
	  int j = kp1[k];
	  int n = kp1[j];
	  m = n_munu[n][j];

	  if( j < n )
	  {
	    q1_prop += f[m] * (Gamma(jm) * quark_prop_1);
	    q2_prop += (anti_quark_prop * Gamma(jm)) * f_source[m];
	  }
	  else
	  {
	    q1_prop -= f[m] * (Gamma(jm) * quark_prop_1);
	    q2_prop -= (anti_quark_prop * Gamma(jm)) * f_source[m];
	  }
	}
 
      LatticeComplex corr_fn = localInnerProduct(q2_prop, q1_prop);
      print_disp(xml_hyb, corr_fn, phases, t0);

      pop(xml_hyb);
    }


    // Exotic 0^{--}: \bar \psi \gamma_5 \gamma_k F_{j_decay,k} \psi
    kv = 5;
    {
      push(xml_hyb);     // next array element
      write(xml_hyb, "kv", kv);

      LatticePropagator q1_prop = 0;
      LatticePropagator q2_prop = 0;

      for(int k = 0; k < Nd; ++k)
	if( k != j_decay )
	{
	  int m = 1 << k;
	  int jm = m ^ G5;
	  m = n_munu[k][j_decay];

	  if( j_decay < k )
	  {
	    q1_prop += f[m] * (Gamma(jm) * quark_prop_1);
	    q2_prop += (anti_quark_prop * Gamma(jm)) * f_source[m];
	  }
	  else
	  {
	    q1_prop -= f[m] * (Gamma(jm) * quark_prop_1);
	    q2_prop -= (anti_quark_prop * Gamma(jm)) * f_source[m];
	  }
	}

      LatticeComplex corr_fn = localInnerProduct(q2_prop, q1_prop);
      print_disp(xml_hyb, corr_fn, phases, t0);

      pop(xml_hyb);
    }


    // Exotic 1^{-+}_1: \bar \psi \gamma_j F_{j,k} \psi
    kv = 5;
    for(int k = 0; k < Nd; ++k)
      if( k != j_decay )
      {
	kv++;

	push(xml_hyb);     // next array element
	write(xml_hyb, "kv", kv);

	LatticePropagator q1_prop = 0;
	LatticePropagator q2_prop = 0;

	for(int j = 0; j < Nd; ++j)
	  if( j != j_decay && j != k )
	  {
	    int jm = 1 << j;
	    int m = n_munu[k][j];

	    if( j < k )
	    {
	      q1_prop += f[m] * (Gamma(jm) * quark_prop_1);
	      q2_prop += (anti_quark_prop * Gamma(jm)) * f_source[m];
	    }
	    else
	    {
	      q1_prop -= f[m] * (Gamma(jm) * quark_prop_1);
	      q2_prop -= (anti_quark_prop * Gamma(jm)) * f_source[m];
	    }
	  }

	LatticeComplex corr_fn = localInnerProduct(q2_prop, q1_prop);
	print_disp(xml_hyb, corr_fn, phases, t0);
      
	pop(xml_hyb);
      }


    // Exotic 1^{-+}_2: \bar \psi \gamma_{j_decay} F_{j_decay,k} \psi
    kv = 8;
    for(int k = 0; k < Nd; ++k)
      if( k != j_decay )
      {
	kv++;

	push(xml_hyb);     // next array element
	write(xml_hyb, "kv", kv);

	int jm = 1 << j_decay;
	int m = n_munu[k][j_decay];

	LatticeComplex corr_fn = localInnerProduct(anti_quark_prop * f_source[m],
						   f[m] * (Gamma(jm) * quark_prop_1 * Gamma(jm)));
	print_disp(xml_hyb, corr_fn, phases, t0);

	pop(xml_hyb);
      }

    /*
     *  Finally, some spin-exotic mesons in which we have non-local
     *  sinks
     *
     *  This way we may find we get a better signal
     */

    /*
     *  As a first test, we will compute the 1^{-+} operators in which we have
     *
     *  O_1 at the source
     *  O_3 = \bar{psi} \gamma_j F_{jk} \psi at the sink
     *
     *  where at the sink we now write F{jk} in terms of the
     *  L-shaped paths a la Lacock et al, PRD54, 6997.
     *
     *  In their notation, the F_{jk} operator corresponds to the path
     *  (jk) - (kj) + (k jbar) - (jbar k) + (jbar kbar) - (kbar jbar)
     *  + (kbar j) - (j kbar)
     */

    kv = 11;			/* Channel counter */
    for(int k = 0; k < Nd; ++k)
      if( k != j_decay )
      {
	kv++;

	push(xml_hyb);     // next array element
	write(xml_hyb,"kv", kv);

	LatticePropagator q1_prop = 0;
	LatticePropagator q2_prop = 0;

	LatticePropagator tmp_prop1;
	LatticePropagator tmp_prop2;

	for(int j = 0; j < Nd; ++j)
	{
	  if( j != j_decay && j != k )
	  {
	    /*
	     *  First we perform the gamma-matrix algebra
	     */
	    int jm = 1 << j;
	    int m = n_munu[k][j];

	    LatticePropagator tmp_prop1 = Gamma(jm) * quark_prop_1;
	  
	    /*
	     *  Now the simple operation of multiplying by
	     *  F_{jk} at the SOURCE
	     */
	    if( j < k )
	    {
	      q2_prop += (anti_quark_prop * Gamma(jm)) * f_source[m];
	    }
	    else
	    {
	      q2_prop -= (anti_quark_prop * Gamma(jm)) * f_source[m];
	    }

	    /*
	     *  Note that we can now reuse tmp_prop2
	     */

	    /*
	     *  At the sink, we have to do various communication
	     * 
	     *  We will loop over the forward and backward directions
	     *  for the j and k axes respectively
	     *
	     *  There are a total of 8 terms
	     */

	    /*
	     * (jk) - (j kbar)
	     */
	    tmp_prop2 = u_smr[k] * shift(tmp_prop1, FORWARD, k) /* Forward k */
	      - shift(adj(u_smr[k]) * tmp_prop1, BACKWARD, k); /* k - kbar */
	    q1_prop += u_smr[j] * shift(tmp_prop2, FORWARD, j); /* j x (k - kbar) */

	    /*
	     * - [ (kj) - (k jbar)]
	     */
	    tmp_prop2 = u_smr[j] * shift(tmp_prop1, FORWARD, j) /* Forward j */
	      - shift(adj(u_smr[j]) * tmp_prop1, BACKWARD, j); /* j - jbar */
	    q1_prop -= u_smr[k] * shift(tmp_prop2, FORWARD, k); /* - k x (j - jbar) */

	    /*
	     * - [ (jbar k) - (jbar  kbar)]
	     */
	    tmp_prop2 = u_smr[k] * shift(tmp_prop1, FORWARD, k) /* Forward k */
	      - shift(adj(u_smr[k]) * tmp_prop1, BACKWARD, k); /* k - kbar */
	    q1_prop -= shift(adj(u_smr[j]) * tmp_prop2, BACKWARD, j); /* Now the communication */

	    /*
	     * + [ (kbar j) - (kbar  jbar)]
	     */
	    tmp_prop2 = u_smr[j] * shift(tmp_prop1, FORWARD, j) /* Forward j */
	      - shift(adj(u_smr[j]) * tmp_prop1, BACKWARD, j); /* k - kbar */
	    q1_prop += shift(adj(u_smr[k]) * tmp_prop2, BACKWARD, k); /* Now the communication */
	  } /* End if statement */
	} /* End loop over j */


	LatticeComplex corr_fn = localInnerProduct(q2_prop, q1_prop);
	print_disp(xml_hyb, corr_fn, phases, t0);

	pop(xml_hyb);
      }

    pop(xml_hyb);

    END_CODE();
  }

}  // end namespace Chroma
