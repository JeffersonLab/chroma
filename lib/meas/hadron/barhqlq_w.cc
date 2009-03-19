// $Id: barhqlq_w.cc,v 3.5 2009-03-19 17:17:20 mcneile Exp $
/*! \file
 *  \brief Heavy-light baryon 2-pt functions
 */


#include "meas/hadron/barhqlq_w.h"
#include "meas/hadron/barspinmat_w.h"

namespace Chroma 
{

  //! Baryon 2pt contractions
  /*! \ingroup hadron */
  namespace  Baryon2PtContractions
  {
    //! Sigma 2-pt
    /*! \ingroup hadron */
    LatticeComplex sigma2pt(const LatticePropagator& quark_propagator_1,
			    const LatticePropagator& quark_propagator_2,
			    const SpinMatrix& T, const SpinMatrix& sp) 
    {
#if QDP_NC == 3

      LatticePropagator di_quark = quarkContract13(quark_propagator_1 * sp,
						   sp * quark_propagator_2);
      return LatticeComplex(trace(T * traceColor(quark_propagator_2 * traceSpin(di_quark)))
			    + trace(T * traceColor(quark_propagator_2 * di_quark)));
#endif
    }
	      

    //! Cascade 2-pt
    /*! \ingroup hadron */
    LatticeComplex xi2pt(const LatticePropagator& quark_propagator_1,
			 const LatticePropagator& quark_propagator_2,
			 const SpinMatrix& T, const SpinMatrix& sp) 
    {
#if QDP_NC == 3

      LatticePropagator di_quark = quarkContract13(quark_propagator_1 * sp,
						   sp * quark_propagator_2);
      return LatticeComplex(trace(T * traceColor(quark_propagator_1 * traceSpin(di_quark)))
			    + trace(T * traceColor(quark_propagator_1 * di_quark)));
#endif
    }
	      

    //! Lambda 2-pt
    /*! \ingroup hadron */
    LatticeComplex lambda2pt(const LatticePropagator& quark_propagator_1,
			     const LatticePropagator& quark_propagator_2,
			     const SpinMatrix& T, const SpinMatrix& sp) 
    {
#if QDP_NC == 3

      // WARNING: I'm not convinced the original SZIN version (or this version) is correct!
      LatticePropagator di_quark = quarkContract13(quark_propagator_2 * sp,
						   sp * quark_propagator_2);

      LatticeComplex b_prop  = trace(T * traceColor(quark_propagator_1 * traceSpin(di_quark)))
	+ trace(T * traceColor(quark_propagator_1 * di_quark));

      di_quark = quarkContract13(quark_propagator_2 * sp,
				 sp * quark_propagator_1);
      b_prop += trace(T * traceColor(quark_propagator_2 * di_quark));

      return b_prop;
#endif
    }


    //! Lambda 2-pt
    /*! \ingroup hadron */
    LatticeComplex lambdaNaive2pt(const LatticePropagator& quark_propagator_1,
				  const LatticePropagator& quark_propagator_2,
				  const SpinMatrix& T, const SpinMatrix& sp) 
    {
#if QDP_NC == 3

      LatticePropagator di_quark = quarkContract13(quark_propagator_2 * sp,
						   sp * quark_propagator_2);
      return LatticeComplex(trace(T * traceColor(quark_propagator_1 * traceSpin(di_quark))));
#endif
    }


    //! Delta 2-pt
    /*! \ingroup hadron */
    LatticeComplex sigmast2pt(const LatticePropagator& quark_propagator_1,
			      const LatticePropagator& quark_propagator_2,
			      const SpinMatrix& T, const SpinMatrix& sp) 
    {
#if QDP_NC == 3

      LatticePropagator di_quark = quarkContract13(quark_propagator_1 * sp, 
						   sp * quark_propagator_2);
      LatticeComplex b_prop = trace(T * traceColor(quark_propagator_2 * traceSpin(di_quark)))
	+ trace(T * traceColor(quark_propagator_2 * di_quark));
      
      di_quark = quarkContract13(quark_propagator_2 * sp, 
				 sp * quark_propagator_1);
      b_prop += trace(T * traceColor(quark_propagator_2 * di_quark));
      
      di_quark = quarkContract13(quark_propagator_2 * sp, 
				 sp * quark_propagator_2);
      b_prop += trace(T * traceColor(quark_propagator_1 * di_quark));
      b_prop *= 2;
      b_prop += trace(T * traceColor(quark_propagator_1 * traceSpin(di_quark)));

      return b_prop;
#endif
    }

    //! Delta 2-pt
    /*! \ingroup hadron */
    LatticeComplex sigmast2pt(const LatticePropagator& quark_propagator_1,
			      const LatticePropagator& quark_propagator_2,
			      const SpinMatrix& T, const SpinMatrix& spSRC,
			      const SpinMatrix& spSNK) 
    {
#if QDP_NC == 3

      LatticePropagator di_quark = quarkContract13(quark_propagator_1 * spSRC, 
						   spSNK * quark_propagator_2);
      LatticeComplex b_prop = trace(T * traceColor(quark_propagator_2 * traceSpin(di_quark)))
	+ trace(T * traceColor(quark_propagator_2 * di_quark));
      
      di_quark = quarkContract13(quark_propagator_2 * spSRC, 
				 spSNK * quark_propagator_1);
      b_prop += trace(T * traceColor(quark_propagator_2 * di_quark));
      
      di_quark = quarkContract13(quark_propagator_2 * spSRC, 
				 spSNK * quark_propagator_2);
      b_prop += trace(T * traceColor(quark_propagator_1 * di_quark));
      b_prop *= 2;
      b_prop += trace(T * traceColor(quark_propagator_1 * traceSpin(di_quark)));
      return b_prop;
#endif
    }

  }  // namespace  Baryon2PtContractions


  //! Heavy-light baryon 2-pt functions
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

   * The routine optionally computes time-charge reversed baryons and adds them
   * in for increased statistics.

   * \param propagator_1   "s" quark propagator ( Read )
   * \param propagator_2   "u" quark propagator ( Read )
   * \param t0             cartesian coordinates of the source ( Read )
   * \param bc_spec        boundary condition for spectroscopy ( Read )
   * \param time_rev       add in time reversed contribution if true ( Read )
   * \param phases         object holds list of momenta and Fourier phases ( Read )
   * \param xml            xml file object ( Read )
   * \param xml_group      group name for xml data ( Read )
   *
   */

  void barhqlq(const LatticePropagator& propagator_1, 
	       const LatticePropagator& propagator_2, 
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
    barhqlq(propagator_1, propagator_2, phases, bardisp1);

    // Possibly add in a time-reversed contribution
    bool time_revP = (bc_spec*bc_spec == 1) ? time_rev : false;

    if (time_revP)
    {
      /* Time-charge reverse the quark propagators */
      /* S_{CT} = gamma_5 gamma_4 = gamma_1 gamma_2 gamma_3 = Gamma(7) */
      LatticePropagator q1_tmp = - (Gamma(7) * propagator_1 * Gamma(7));
      LatticePropagator q2_tmp = - (Gamma(7) * propagator_2 * Gamma(7));

      barhqlq(q1_tmp, q2_tmp, phases, bardisp2);
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



  //! Heavy-light baryon 2-pt functions
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions! 
   *
   *###########################################################################
   * WARNING: No symmetrization over the spatial part of the wave functions   #
   *          is performed. Therefore, if this routine is called with         #
   *          "shell-sink" quark propagators of different widths the          #
   *          resulting octet baryons may have admixters of excited           #
   *          decouplet baryons with mixed symmetric spatial wave functions,  #
   *          and vice-versa!!!                                               #
   *###########################################################################

   * Construct heavy-light baryon propagators with two "u" quarks and
   * one separate "s" quark for the Sigma^+, the Lambda and the Sigma^{*+}.
   * In the Lambda we take the "u" and "d" quark as degenerate!

   * The routine also computes time-charge reversed baryons and adds them
   * in for increased statistics.

   * \param quark_propagator_1   "s" quark propagator ( Read )
   * \param quark_propagator_2   "u" quark propagator ( Read )
   * \param barprop              baryon propagator ( Modify )
   * \param phases               object holds list of momenta and Fourier phases ( Read )
   *
   *        ____
   *        \
   * b(t) =  >  < b(t_source, 0) b(t + t_source, x) >
   *        /                    
   *        ----
   *          x

   * For the Sigma^+ we take

   * |S_1, s_z=1/2> = (s C gamma_5 u) "u_up"

   * for the Lambda

   * |L_1, s_z=1/2> = 2*(u C gamma_5 d) "s_up" + (s C gamma_5 d) "u_up"
   *                  + (u C gamma_5 s) "d_up"

   * and for the Sigma^{*+}

   * |S*_1, s_z=3/2> = 2*(s C gamma_- u) "u_up" + (u C gamma_- u) "s_up".

   * We have put "q_up" in quotes, since this is meant in the Dirac basis,
   * not in the 'DeGrand-Rossi' chiral basis used in the program!
   * In gamma_- we ignore a factor sqrt(2).

   * For all baryons we compute a 'B_2' that differs from the 'B_1' above
   * by insertion of a gamma_4 between C and the gamma_{5,-}.
   * And finally, we also compute the non-relativistic baryons, 'B_3',
   * which up to a factor 1/2 are just the difference B_1 - B_2, as can
   * be seen by projecting to the "upper" components in the Dirac basis,
   * achieved by (1 + gamma_4)/2 q, for quark q.

   * The Sigma^+_k is baryon 3*(k-1), the Lambda_k is baryon 3*(k-1)+1
   * and the Sigma^{*+}_k is baryon 3*(k-1)+2.

   * We are using a chiral basis for the Dirac matrices (gamma_5 diagonal).
   * Therefore a spin-up quark in the Dirac basis corresponds to
   * 1/sqrt(2) * ( - q_1 - q_3 ) in this chiral basis. We shall neglect
   * the sign and the 1/sqrt(2) here.
   * The projection on "spin_up" is done with the projector "T". 
   */

  void barhqlq(const LatticePropagator& quark_propagator_1,
	       const LatticePropagator& quark_propagator_2,
	       const SftMom& phases,
	       multi3d<DComplex>& barprop)
  {
    START_CODE();

    // Length of lattice in decay direction
    int length = phases.numSubsets() ;

    if ( Ns != 4 || Nc != 3 )		/* Code is specific to Ns=4 and Nc=3. */
      return;

    // Setup the return stuff
    const int num_baryons = 17;
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

    LatticeComplex b_prop;

    // Loop over baryons
    for(int baryons = 0; baryons < num_baryons; ++baryons)
    {
      LatticePropagator di_quark;

      switch (baryons)
      {
      case 0:
	// Sigma^+_1 (or proton); use also for Lambda_1!
	// |S_1, s_z=1/2> = (s C gamma_5 u) "u_up"
	// C gamma_5 = Gamma(5)
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	b_prop = Baryon2PtContractions::sigma2pt(quark_propagator_1, quark_propagator_2, 
						 T_mixed, Cg5);
	break;

      case 1:
	// Lambda_1
	// |L_1, s_z=1/2> = 2*(u C gamma_5 d) "s_up" + (s C gamma_5 d) "u_up"
	//                  + (u C gamma_5 s) "d_up" , see comments at top   
	// C gamma_5 = Gamma(5)
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	b_prop = Baryon2PtContractions::lambda2pt(quark_propagator_1, quark_propagator_2, 
						  T_mixed, Cg5);
	break;

      case 2:
	// Sigma^{*+}_1
	// |S*_1, s_z=3/2> = 2*(s C gamma_- u) "u_up" + (u C gamma_- u) "s_up"
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	b_prop = Baryon2PtContractions::sigmast2pt(quark_propagator_1, quark_propagator_2, 
						   T_mixed, BaryonSpinMats::Cgm());
	break;

      case 3:
	// Sigma^+_2; use also for Lambda_2!
	// |S_2, s_z=1/2> = (s C gamma_4 gamma_5 u) "u_up"
	// C gamma_5 gamma_4 = - Gamma(13)
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	b_prop = Baryon2PtContractions::sigma2pt(quark_propagator_1, quark_propagator_2, 
						 T_mixed, Cg5g4);
	break;

      case 4:
	// Lambda_2
	// |L_2, s_z=1/2> = 2*(u C gamma_4 gamma_5 d) "s_up"
	//                  + (s C gamma_4 gamma_5 d) "u_up"
	//                  + (u C gamma_4 gamma_5 s) "d_up"
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	b_prop = Baryon2PtContractions::lambda2pt(quark_propagator_1, quark_propagator_2, 
						  T_mixed, Cg5g4);
	break;

      case 5:
	// Sigma^{*+}_2
	// |S*_2, s_z=3/2> = 2*(s C gamma_4 gamma_- u) "u_up" + (u C gamma_4 gamma_- u) "s_up"
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	b_prop = Baryon2PtContractions::sigmast2pt(quark_propagator_1, quark_propagator_2, 
						   T_mixed, BaryonSpinMats::Cg4m());
	break;

      case 6:
	// Sigma^+_3; use also for Lambda_3!
	// |S_3, s_z=1/2> = (s C (1/2)(1 + gamma_4) gamma_5 u) "u_up"
	// C gamma_5 gamma_4 = - Gamma(13)
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	b_prop = Baryon2PtContractions::sigma2pt(quark_propagator_1, quark_propagator_2, 
						 T_mixed, Cg5NR);
	break;

      case 7:
	// Lambda_3
	// |L_3, s_z=1/2> = 2*(u C (1/2)(1 + gamma_4) gamma_5 d) "s_up"
	//                  + (s C (1/2)(1 + gamma_4) gamma_5 d) "u_up"
	//                  + (u C (1/2)(1 + gamma_4) gamma_5 s) "d_up"
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	b_prop = Baryon2PtContractions::lambda2pt(quark_propagator_1, quark_propagator_2, 
						  T_mixed, Cg5NR);
	break;

      case 8:
	// Sigma^{*+}_3
	// |S*_3, s_z=3/2> = 2*(s C (1/2)(1+gamma_4) gamma_- u) "u_up" 
	//                   + (u C (1/2)(1+gamma_4) gamma_- u) "s_up"
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	// Arrgh, goofy CgmNR normalization again from szin code. 
	b_prop = Baryon2PtContractions::sigmast2pt(quark_propagator_1, quark_propagator_2, 
						   T_mixed, BaryonSpinMats::CgmNR());

	// Agghh, we have a goofy factor of 4 normalization factor here. The
	// ancient szin way didn't care about norms, so it happily made it
	// 4 times too big. There is a missing 0.5 in the NR normalization
	// in the old szin code.
	// So, we compensate to keep the same normalization
	b_prop *= 4.0;
	break;


      case 9:
	// Sigma^+_4 -- but unpolarised
	// |S_4, s_z=1/2> = (s C gamma_5 u) "u_up", see comments at top
	// C gamma_5 = Gamma(5)
	// Unpolarized:
	// T_unpol = T = (1/2)(1 + gamma_4)
	b_prop = Baryon2PtContractions::sigma2pt(quark_propagator_1, quark_propagator_2, 
						 T_unpol, Cg5);
	break;

      case 10:
	// Sigma^+_5
	// |S_5, s_z=1/2> = (s C gamma_4 gamma_5 u) "u_up", see comments at top
	// C gamma_5 gamma_4 = - Gamma(13)
	// Unpolarized:
	// T_unpol = T = (1/2)(1 + gamma_4)
	b_prop = Baryon2PtContractions::sigma2pt(quark_propagator_1, quark_propagator_2, 
						 T_unpol, Cg5g4);
	break;
    
      case 11:
	// Sigma^+_6
	// |S_6, s_z=1/2> = (s C (1/2)(1 + gamma_4) gamma_5 u) "u_up", see comments at top
	// C gamma_5 = Gamma(5)
	// Unpolarized:
	// T_unpol = T = (1/2)(1 + gamma_4)
	b_prop = Baryon2PtContractions::sigma2pt(quark_propagator_1, quark_propagator_2, 
						 T_unpol, Cg5NR);
	break;

      case 12:
	// Lambda_4 : naive Lambda interpolating field
	// |L_4 > = (d C gamma_5 u) s
	// C gamma_5 = Gamma(5)
	// UnPolarized:
	// T_unpol = T = (1/2)(1 + gamma_4)
	b_prop = Baryon2PtContractions::lambdaNaive2pt(quark_propagator_1, quark_propagator_2, 
						       T_unpol, Cg5);
	break;
      
      case 13:
	// Xi_1
	// |X_1 > = (s C gamma_5 u) s
	// C gamma_5 = Gamma(5)
	// UnPolarized:
	// T_unpol = T = (1/2)(1 + gamma_4)
	b_prop = Baryon2PtContractions::xi2pt(quark_propagator_1, quark_propagator_2, 
					      T_unpol, Cg5);
	break;

      case 14:
	// Lambda_5 : naive Lambda interpolating field
	// |L_5 > = (d C gamma_5 u) "s_up"
	// C gamma_5 = Gamma(5)
	// UnPolarized: 
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	b_prop = Baryon2PtContractions::lambdaNaive2pt(quark_propagator_1, quark_propagator_2, 
						       T_unpol, Cg5);
	break;
      
      case 15:
	// Xi_2
	// |X_2 > = (s C gamma_5 u) "s_up"
	// C gamma_5 = Gamma(5)
	// UnPolarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	b_prop = Baryon2PtContractions::xi2pt(quark_propagator_1, quark_propagator_2, 
					      T_mixed, Cg5);
	break;

      case 16:
	// Proton_negpar_3; use also for Lambda_negpar_3!
	// |P_7, s_z=1/2> = (d C gamma_5 (1/2)(1 - g_4) u) "u_up", see comments at top
	// C g_5 NR negpar = (1/2)*C gamma_5 * ( 1 - g_4 )
	// T = (1 + \Sigma_3)*(1 - gamma_4) / 2 
	//   = (1 - Gamma(8) + i G(3) - i G(11)) / 2
	b_prop = Baryon2PtContractions::sigma2pt(quark_propagator_1, quark_propagator_2, 
						 BaryonSpinMats::TmixedNegPar(), 
						 BaryonSpinMats::Cg5NRnegPar());
	break;
		  
      default:
	QDP_error_exit("Unknown baryon", baryons);
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


