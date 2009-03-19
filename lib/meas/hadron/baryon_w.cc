// $Id: baryon_w.cc,v 3.1 2009-03-19 17:17:20 mcneile Exp $ 
/*! \file
 *  \brief Baryon 2-pt functions
 */

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/baryon_w.h"
#include "meas/hadron/barspinmat_w.h"

namespace Chroma 
{

  //! Baryon 2-pt functions
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions! 
   *
   * Construct baryon propagators for the Proton and the Delta^+ with
   * degenerate "u" and "d" quarks, as well as the Lambda for, in
   * addition, a degenerate "s" quark. For these degenerate quarks, the
   * Lambda is degenerate with the Proton, but we keep it for compatibility
   * with the sister routine that treats non-degenerate quarks.
   *
   * The routine optionally computes time-charge reversed baryons and adds them
   * in for increased statistics.
   *
   * \param quark_propagator   quark propagator ( Read )
   * \param t0         cartesian coordinates of the source ( Read )
   * \param bc_spec    boundary condition for spectroscopy ( Read )
   * \param time_rev   add in time reversed contribution if true ( Read )
   * \param phases     object holds list of momenta and Fourier phases ( Read )
   * \param xml        xml file object ( Read )
   * \param xml_group  group name for xml data ( Read )
   *
   */

  void baryon(const LatticePropagator& quark_propagator, 
	      const SftMom& phases,
	      int t0, int bc_spec, bool time_rev,
	      XMLWriter& xml,
	      const string& xml_group)
  {
    START_CODE();

    if ( Ns != 4 || Nc != 3 )		/* Code is specific to Ns=4 and Nc=3. */
      return;

    multi3d<DComplex> bardisp1;
    multi3d<DComplex> bardisp2;

    // Forward
    baryon(quark_propagator, phases, bardisp1);

    // Possibly add in a time-reversed contribution
    bool time_revP = (bc_spec*bc_spec == 1) ? time_rev : false;

    if (time_revP)
    {
      /* Time-charge reverse the quark propagators */
      /* S_{CT} = gamma_5 gamma_4 = gamma_1 gamma_2 gamma_3 = Gamma(7) */
      LatticePropagator q1_tmp = - (Gamma(7) * quark_propagator * Gamma(7));

      baryon(q1_tmp, phases, bardisp2);
    }


    int num_baryons = bardisp1.size3();
    int num_mom = bardisp1.size2();
    int length  = bardisp1.size1();

    // Loop over baryons
    XMLArrayWriter xml_bar(xml,num_baryons);
    push(xml_bar, xml_group);

    for(int baryons = 0; baryons < num_baryons; ++baryons)
    {
      push(xml_bar);     // next array element
      write(xml_bar, "baryon_num", baryons);

      // Loop over sink momenta
      XMLArrayWriter xml_sink_mom(xml_bar,num_mom);
      push(xml_sink_mom, "momenta");

      for(int sink_mom_num = 0; sink_mom_num < num_mom; ++sink_mom_num)
      {
	push(xml_sink_mom);
	write(xml_sink_mom, "sink_mom_num", sink_mom_num) ;
	write(xml_sink_mom, "sink_mom", phases.numToMom(sink_mom_num)) ;

	multi1d<Complex> barprop(length);

	/* forward */
	for(int t = 0; t < length; ++t)
	{
	  int t_eff = (t - t0 + length) % length;
	    
	  if ( bc_spec < 0 && (t_eff+t0) >= length)
	    barprop[t_eff] = -bardisp1[baryons][sink_mom_num][t];
	  else
	    barprop[t_eff] =  bardisp1[baryons][sink_mom_num][t];
	}

	if (time_revP)
	{
	  /* backward */
	  for(int t = 0; t < length; ++t)
	  {
	    int t_eff = (length - t + t0) % length;
	
	    if ( bc_spec < 0 && (t_eff-t0) > 0)
	    {
	      barprop[t_eff] -= bardisp2[baryons][sink_mom_num][t];
	      barprop[t_eff] *= 0.5;
	    }
	    else
	    {
	      barprop[t_eff] += bardisp2[baryons][sink_mom_num][t];
	      barprop[t_eff] *= 0.5;
	    }
	  }
	}

	write(xml_sink_mom, "barprop", barprop);
	pop(xml_sink_mom);
      } // end for(sink_mom_num)
 
      pop(xml_sink_mom);
      pop(xml_bar);
    } // end for(gamma_value)

    pop(xml_bar);

    END_CODE();
  }


  //! Nucleon 2-pt
  /*! \ingroup hadron */
  LatticeComplex nucl2pt(const LatticePropagator& quark_propagator,
			 const SpinMatrix& T, const SpinMatrix& sp) 
  {
#if QDP_NC == 3

    LatticePropagator di_quark = quarkContract13(quark_propagator * sp,
						 sp * quark_propagator);
    return LatticeComplex(trace(T * traceColor(quark_propagator * traceSpin(di_quark)))
                        + trace(T * traceColor(quark_propagator * di_quark)));
#else
    LatticeComplex a ; 
    a = zero ;
    return a ;
#endif
  }
	      

  //! Delta 2-pt
  /*! \ingroup hadron */
  LatticeComplex delta2pt(const LatticePropagator& quark_propagator,
			  const SpinMatrix& T, const SpinMatrix& sp) 
  {
#if QDP_NC == 3
    LatticePropagator di_quark = quarkContract13(quark_propagator * sp,
						 sp * quark_propagator);
    return LatticeComplex(trace(T * traceColor(quark_propagator * traceSpin(di_quark)))
		      + 2*trace(T * traceColor(quark_propagator * di_quark)));

#else
    LatticeComplex a ; 
    a = zero ;
    return a ;
#endif


  }



  //! Baryon 2-pt functions
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions! 
   *
   * Construct baryon propagators for the Proton and the Delta^+ with
   * degenerate "u" and "d" quarks, as well as the Lambda for, in
   * addition, a degenerate "s" quark. For these degenerate quarks, the
   * Lambda is degenerate with the Proton, but we keep it for compatibility
   * with the sister routine that treats non-degenerate quarks.
   *
   * \param quark_propagator  quark propagator ( Read )
   * \param barprop    baryon propagator ( Modify )
   * \param phases     object holds list of momenta and Fourier phases ( Read )
   *
   *        ____
   *        \
   * b(t) =  >  < b(t_source, 0) b(t + t_source, x) >
   *        /                    
   *        ----
   *          x

   * For the Proton we take

   * |P_1, s_z=1/2> = (d C gamma_5 u) "u_up"

   * for the Lambda

   * |L_1, s_z=1/2> = 2*(u C gamma_5 d) "s_up" + (s C gamma_5 d) "u_up"
   *                  + (u C gamma_5 s) "d_up"

   * and for the Delta^+

   * |D_1, s_z=3/2> = 2*(d C gamma_- u) "u_up" + (u C gamma_- u) "d_up".

   * We have put "q_up" in quotes, since this is meant in the Dirac basis,
   * not in the 'DeGrand-Rossi' chiral basis used in the program!

   * For all baryons we compute a 'B_2' that differs from the 'B_1' above
   * by insertion of a gamma_4 between C and the gamma_{5,-}.
   * And finally, we also compute the non-relativistic baryons, 'B_3',
   * which up to a factor 1/2 are just the difference B_1 - B_2, as can
   * be seen by projecting to the "upper" components in the Dirac basis,
   * achieved by (1 + gamma_4)/2 q, for quark q.

   * The Proton_k is baryon 3*(k-1), the Lambda_k is baryon 3*(k-1)+1
   * and the Delta^+_k is baryon 3*(k-1)+2. 
   */

  void baryon(const LatticePropagator& quark_propagator, 
	      const SftMom& phases,
	      multi3d<DComplex>& barprop)
  {
    START_CODE();

    // Length of lattice in decay direction
    int length = phases.numSubsets() ;

    if ( Ns != 4 || Nc != 3 )		/* Code is specific to Ns=4 and Nc=3. */
      return;

    // Setup the return stuff
    const int num_baryons = 22;
    int num_mom = phases.numMom();
    barprop.resize(num_baryons,num_mom,length);

    // T_mixed = (1 + \Sigma_3)*(1 + gamma_4) / 2 
    //         = (1 + Gamma(8) - i G(3) - i G(11)) / 2
    SpinMatrix T_mixed = BaryonSpinMats::Tmixed();

    // T_unpol = (1/2)(1 + gamma_4)
    SpinMatrix T_unpol = BaryonSpinMats::Tunpol();

    // C gamma_5 = Gamma(5)
    SpinMatrix Cg5 = BaryonSpinMats::Cg5();

    // C gamma_5 gamma_4 = - Gamma(13)
    SpinMatrix Cg5g4 = BaryonSpinMats::Cg5g4();

    // C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 )
    SpinMatrix Cg5NR = BaryonSpinMats::Cg5NR();

    // C = Gamma(10)
    SpinMatrix C = BaryonSpinMats::C();

    LatticeComplex b_prop;

    // Loop over baryons
    for(int baryons = 0; baryons < num_baryons; ++baryons)
    {
      switch (baryons)
      {
      case 0:
	// Proton_1; use also for Lambda_1!
	// |P_1, s_z=1/2> = (d C gamma_5 u) "u_up", see comments at top
	// C gamma_5 = Gamma(5)
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	b_prop = nucl2pt(quark_propagator, T_mixed, Cg5);
	break;
		  
      case 1:
	// Lambda_1 = 3*Proton_1 (for compatibility with heavy-light routine)
	// |L_1, s_z=1/2> = 2*(u C gamma_5 d) "s_up" + (s C gamma_5 d) "u_up"
	//                  + (u C gamma_5 s) "d_up" , see comments at top   
	// C gamma_5 = Gamma(5)
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	b_prop *= 3.0;
	break;

      case 2:
	// Delta^+_1
	// |D_1, s_z=3/2> = 2*(d C gamma_- u) "u_up" + (u C gamma_- u) "d_up"
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	// Multiply by 3 for compatibility with heavy-light routine
	b_prop = 3.0 * delta2pt(quark_propagator, T_mixed, BaryonSpinMats::Cgm());
	break;

      case 3:
	// Proton_2; use also for Lambda_2!
	// |P_2, s_z=1/2> = (d C gamma_4 gamma_5 u) "u_up" 
	// C gamma_5 gamma_4 = - Gamma(13)
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	b_prop = nucl2pt(quark_propagator, T_mixed, Cg5g4);
	break;

      case 4:
	// Lambda_2 = 3*Proton_2 (for compatibility with heavy-light routine)
	// |L_2, s_z=1/2> = 2*(u C gamma_4 gamma_5 d) "s_up"
	//                  + (s C gamma_4 gamma_5 d) "u_up"
	//                  + (u C gamma_4 gamma_5 s) "d_up"
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	b_prop *= 3.0;
	break;

      case 5:
	// Sigma^{*+}_2
	// |D_2, s_z=3/2> = 2*(d C gamma_4 gamma_- u) "u_up" 
	//                  + (u C gamma_4 gamma_- u) "d_up" 
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
	//            = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	// Multiply by 3 for compatibility with heavy-light routine
	b_prop = 3.0 * delta2pt(quark_propagator, T_mixed, BaryonSpinMats::Cg4m());
	break;

      case 6:
	// Proton^+_3; use also for Lambda_3!
	// |P_3, s_z=1/2> = (d C (1/2)(1 + gamma_4) gamma_5 u) "u_up" 
	// C gamma_5 - C gamma_5 gamma_4 = Gamma(5) + Gamma(13)
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	b_prop = nucl2pt(quark_propagator, T_mixed, Cg5NR);
	break;

      case 7:
	// Lambda_3 = 3*Proton_3 (for compatibility with heavy-light routine)
	// |L_3, s_z=1/2> = 2*(u C (1/2)(1 + gamma_4) gamma_5 d) "s_up"
	//                  + (s C (1/2)(1 + gamma_4) gamma_5 d) "u_up"
	//                  + (u C (1/2)(1 + gamma_4) gamma_5 s) "d_up"
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	b_prop *= 3.0;
	break;

      case 8:
	// Sigma^{*+}_3
	// |D_3, s_z=3/2> = 2*(d C (1/2)(1 + gamma_4) gamma_- d) u) "u_up"
	//                  + (u C (1/2)(1 + gamma_4) gamma_- d) u) "d_up"
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	// Multiply by 3 for compatibility with heavy-light routine
	b_prop = 3.0 * delta2pt(quark_propagator, T_mixed, BaryonSpinMats::CgmNR());

	// Agghh, we have a goofy factor of 4 normalization factor here. The
	// ancient szin way didn't care about norms, so it happily made it
	// 4 times too big. There is a missing 0.5 in the NR normalization
	// in the old szin code.
	// So, we compensate to keep the same normalization
	b_prop *= 4.0;
	break;

      case 9:
	// Proton_4 -- but unpolarised ; use also for Lambda_4!
	// |P_4, s_z=1/2> = (d C gamma_5 u) "u_up", see comments at top
	// C gamma_5 = Gamma(5)
	// Unpolarized:
	// T_unpol = T = (1/2)(1 + gamma_4)
	b_prop = nucl2pt(quark_propagator, T_unpol, Cg5);
	break;

      case 10:
	// Proton_5; use also for Lambda_5!
	// |P_5, s_z=1/2> = (d C gamma_4 gamma_5 u) "u_up", see comments at top
	// C gamma_5 gamma_4 = - Gamma(13)
	// Unpolarized:
	// T_unpol = T = (1/2)(1 + gamma_4)
	b_prop = nucl2pt(quark_propagator, T_mixed, Cg5g4);
	break;
    
      case 11:
	// Proton^+_6; use also for Lambda_6!
	// |P_6, s_z=1/2> = (d C (1/2)(1 + gamma_4) gamma_5 u) "u_up", see comments at top
	// C gamma_5 = Gamma(5)
	// Unpolarized:
	// T_unpol = T = (1/2)(1 + gamma_4)
	b_prop = nucl2pt(quark_propagator, T_unpol, Cg5NR);
	break;

      case 12:
	// Delta_x^+_4 -- unpolarised with explicit gamma_k interpolation
	// |D_4, s_z=3/2> = 2*(d C gamma_1 u) "u_up" + (u C gamma_1 u) "d_up"
	// C gamma_1 = Gamma(10) * Gamma(1) = Gamma(11)
	// Unpolarized:
	// T_unpol = T = (1/2)(1 + gamma_4)
	// Multiply by 3 for compatibility with heavy-light routine
	b_prop = 3.0 * delta2pt(quark_propagator, T_unpol, BaryonSpinMats::Cgk(1));
	break;

      case 13:
	// Delta_y^+_4 -- unpolarised with explicit gamma_k interpolation
	// |D_4, s_z=3/2> = 2*(d C gamma_2 u) "u_up" + (u C gamma_2 u) "d_up"
	// C gamma_2 = Gamma(10) * Gamma(2) = Gamma(8)
	// Unpolarized:
	// T_unpol = T = (1/2)(1 + gamma_4)
	// Multiply by 3 for compatibility with heavy-light routine
	b_prop = 3.0 * delta2pt(quark_propagator, T_unpol, BaryonSpinMats::Cgk(2));
	break;

      case 14:
	// Delta_z^+_4 -- unpolarised with explicit gamma_k interpolation
	// |D_4, s_z=3/2> = 2*(d C gamma_3 u) "u_up" + (u C gamma_3 u) "d_up"
	// C gamma_3 = Gamma(10) * Gamma(4) = Gamma(14)
	// Unpolarized:
	// T_unpol = T = (1/2)(1 + gamma_4)
	// Multiply by 3 for compatibility with heavy-light routine
	b_prop = 3.0 * delta2pt(quark_propagator, T_unpol, BaryonSpinMats::Cgk(3));
	break;

      case 15:
	// Delta_x^+_5 -- unpolarised with explicit gamma_k interpolation
	// |D_5, s_z=3/2> = 2*(d C gamma_4 gamma_1 u) "u_up" + (u C gamma_4 gamma_1 u) "d_up"
	// C gamma_4 gamma_1 = Gamma(10) * Gamma(8) * Gamma(1) = Gamma(3)
	// Unpolarized:
	// T_unpol = T = (1/2)(1 + gamma_4)
	// Multiply by 3 for compatibility with heavy-light routine
	b_prop = 3.0 * delta2pt(quark_propagator, T_unpol, BaryonSpinMats::Cg4gk(1));
	break;

      case 16:
	// Delta_y^+_5 -- unpolarised with explicit gamma_k interpolation
	// |D_4, s_z=3/2> = 2*(d C gamma_4 gamma_2 u) "u_up" + (u C gamma_4 gamma_2 u) "d_up"
	// C gamma_4 gamma_2 = Gamma(10) * Gamma(8) * Gamma(2) = Gamma(0)
	// Unpolarized:
	// T_unpol = T = (1/2)(1 + gamma_4)
	// Multiply by 3 for compatibility with heavy-light routine
	b_prop = 3.0 * delta2pt(quark_propagator, T_unpol, BaryonSpinMats::Cg4gk(2));
	break;

      case 17:
	// Delta_z^+_5 -- unpolarised with explicit gamma_k interpolation
	// |D_4, s_z=3/2> = 2*(d C gamma_4 gamma_3 u) "u_up" + (u C gamma_4 gamma_3 u) "d_up"
	// C gamma_4 gamma_3 = Gamma(10) * Gamma(8) * Gamma(4) = Gamma(6)
	// Unpolarized:
	// T_unpol = T = (1/2)(1 + gamma_4)
	// Multiply by 3 for compatibility with heavy-light routine
	b_prop = 3.0 * delta2pt(quark_propagator, T_unpol, BaryonSpinMats::Cg4gk(3));
	break;

      case 18:
	// Delta_x^+_6 -- unpolarised NR with explicit gamma_k interpolation
	// |D_6, s_z=3/2> = 2*(d C gamma_1 (1/2)(1 + gamma_4) d) u) "u_up"
	//                  + (u C gamma_1 (1/2)(1 + gamma_4) d) u) "d_up"
	// Unpolarized:
	// T_unpol = T = (1/2)(1 + gamma_4)
	// Multiply by 3 for compatibility with heavy-light routine
	b_prop = 3.0 * delta2pt(quark_propagator, T_unpol, BaryonSpinMats::CgkNR(1));
	break;

      case 19:
	// Delta_y^+_6 -- unpolarised NR with explicit gamma_k interpolation
	// |D_6, s_z=3/2> = 2*(d C gamma_2 (1/2)(1 + gamma_4) d) u) "u_up"
	//                  + (u C gamma_2 (1/2)(1 + gamma_4) d) u) "d_up"
	// Unpolarized:
	// T_unpol = T = (1/2)(1 + gamma_4)
	// Multiply by 3 for compatibility with heavy-light routine
	b_prop = 3.0 * delta2pt(quark_propagator, T_unpol, BaryonSpinMats::CgkNR(2));
	break;

      case 20:
	// Delta_z^+_6 -- unpolarised NR with explicit gamma_k interpolation
	// |D_6, s_z=3/2> = 2*(d C gamma_3 (1/2)(1 + gamma_4) d) u) "u_up"
	//                  + (u C gamma_3 (1/2)(1 + gamma_4) d) u) "d_up"
	// Unpolarized:
	// T_unpol = T = (1/2)(1 + gamma_4)
	// Multiply by 3 for compatibility with heavy-light routine
	b_prop = 3.0 * delta2pt(quark_propagator, T_unpol, BaryonSpinMats::CgkNR(3));
	break;

      case 21:
	// Proton_negpar_3; use also for Lambda_negpar_3!
	// |P_7, s_z=1/2> = (d C gamma_5 (1/2)(1 - g_4) u) "u_up", see comments at top
	// C g_5 NR negpar = (1/2)*C gamma_5 * ( 1 - g_4 )
	// T = (1 + \Sigma_3)*(1 - gamma_4) / 2 
	//   = (1 - Gamma(8) + i G(3) - i G(11)) / 2
	b_prop = nucl2pt(quark_propagator, 
			 BaryonSpinMats::TmixedNegPar(), BaryonSpinMats::Cg5NRnegPar());
	break;
		  
      default:
	QDP_error_exit("Unknown baryon: baryons=%d",baryons);
      }

      // Project onto zero and if desired non-zero momentum
      multi2d<DComplex> hsum;
      hsum = phases.sft(b_prop);

      for(int sink_mom_num=0; sink_mom_num < num_mom; ++sink_mom_num) 
	for(int t = 0; t < length; ++t)
	{
	  // NOTE: there is NO  1/2  multiplying hsum
	  barprop[baryons][sink_mom_num][t] = hsum[sink_mom_num][t];
	}

    } // end loop over baryons

    END_CODE();
  }

}  // end namespace Chroma

