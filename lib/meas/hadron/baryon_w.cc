// $Id: baryon_w.cc,v 1.4 2003-02-16 04:14:37 edwards Exp $ 
/*! \file
 *  \brief Baryon 2-pt functions
 */

#include "chromabase.h"
#include "meas/hadron/baryon_w.h"

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

 
//! Baryon 2-pt functions
/*!
 * This routine is specific to Wilson fermions! 
 *
 * Construct baryon propagators for the Proton and the Delta^+ with
 * degenerate "u" and "d" quarks, as well as the Lambda for, in
 * addition, a degenerate "s" quark. For these degenerate quarks, the
 * Lambda is degenerate with the Proton, but we keep it for compatibility
 * with the sister routine that treats non-degenerate quarks.

 * The routine also computes time-charge reversed baryons and adds them
 * in for increased statistics.

 * \param quark_propagator -- quark propagator ( Read )
 * \param barprop -- baryon propagator ( Modify )
 * \param bardisp -- baryon props. at non-zero momenta ( Modify )
 * \param num_mom -- number of non-zero momenta ( Read )
 * \param t_source -- cartesian coordinates of the source ( Read )
 * \param j_decay -- direction of the exponential decay ( Read ) 
 * \param bc_spec  -- boundary condition for spectroscopy ( Read )
 *
 * FftP    -- flag for use of fft or sft ( Read )
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

void baryon(LatticePropagator& quark_propagator, 
	    multi2d<Complex>& barprop, 
	    const multi1d<int>& t_source, int j_decay, int bc_spec)
{
  // Create the time-slice set
  Set timeslice;
  timeslice.make(TimeSliceFunc(j_decay));

  // Length of lattice in j_decay direction
  int length = timeslice.numSubsets();

  if ( Ns != 4 || Nc != 3 )		/* Code is specific to Ns=4 and Nc=3. */
    return;

  int t0 = t_source[j_decay];
  
  SpinMatrix Cgm;
  SpinMatrix Cg4m;
  SpinMatrix CgmNR;

  SpinMatrix g_one = 1.0;
  SpinMatrix g_tmp1;

  /* C = Gamma(10) */
  g_tmp1 = 0.5 * (Gamma(2) * g_one  +  timesI(Gamma(1) * g_one));
  Cgm = Gamma(10) * g_tmp1;

  Cg4m = Gamma(10) * (Gamma(8) * g_tmp1);
  CgmNR = Cgm - Cg4m;

  SpinMatrix S_proj = 
    0.5*((g_one + Gamma(8) * g_one) - timesI(Gamma(3) * g_one  +  Gamma(11) * g_one));

  /*Loop over time-charge reversals */
  for(int time_rev = 0; time_rev < 2; ++time_rev)
  {
    LatticeComplex b_prop;

    /* Loop over baryons */
    for(int baryons = 0; baryons < 9; ++baryons)
    {
      LatticePropagator di_quark;

      switch (baryons)
      {
      case 0:
        /* Proton_1; use also for Lambda_1! */
	/* C gamma_5 = Gamma(5) */
	di_quark = quarkContract13(quark_propagator * Gamma(5),
				   Gamma(5) * quark_propagator);
	b_prop = trace(S_proj * traceColor(quark_propagator * traceSpin(di_quark)))
	  + trace(S_proj * traceColor(quark_propagator * di_quark));
	break;
		  
      case 1:
        /* Lambda_1 = 3*Proton_1 (for compatibility with heavy-light routine) */
	b_prop *= 3.0;
	break;

      case 2:
	/* Delta^+_1 */
	di_quark = quarkContract13(quark_propagator * Cgm, 
				   Cgm * quark_propagator);
	b_prop = trace(S_proj * traceColor(quark_propagator * traceSpin(di_quark)))
	  + 2*trace(S_proj * traceColor(quark_propagator * di_quark));

	/* Multiply by 3 for compatibility with heavy-light routine */
	b_prop *= 3.0;
	break;

      case 3:
        /* Proton_2; use also for Lambda_2! */
	/* C gamma_5 gamma_4 = - Gamma(13) */
	di_quark = quarkContract13(quark_propagator * Gamma(13),
				   Gamma(13) * quark_propagator);
	b_prop = trace(S_proj * traceColor(quark_propagator * traceSpin(di_quark)))
	  + trace(S_proj * traceColor(quark_propagator * di_quark));
	break;

      case 4:
        /* Lambda_2 = 3*Proton_2 (for compatibility with heavy-light routine) */
	b_prop *= 3.0;
	break;

      case 5:
        /* Sigma^{*+}_2 */
	di_quark = quarkContract13(quark_propagator * Cg4m,
				   Cg4m * quark_propagator);
	b_prop = trace(S_proj * traceColor(quark_propagator * traceSpin(di_quark)))
	  + 2*trace(S_proj * traceColor(quark_propagator * di_quark));

	/* Multiply by 3 for compatibility with heavy-light routine */
	b_prop *= 3.0;
	break;

      case 6:
        /* Proton^+_3; use also for Lambda_3! */
	/* C gamma_5 - C gamma_5 gamma_4 = Gamma(5) + Gamma(13) */
	di_quark = quarkContract13(quark_propagator * Gamma(5) + quark_propagator * Gamma(13),  
				   Gamma(5) * quark_propagator + Gamma(13) * quark_propagator);
	b_prop = trace(S_proj * traceColor(quark_propagator * traceSpin(di_quark)))
	  + trace(S_proj * traceColor(quark_propagator * di_quark));
	break;

      case 7:
        /* Lambda_3 = 3*Proton_3 (for compatibility with heavy-light routine) */
	b_prop *= 3.0;
	break;

      case 8:
        /* Sigma^{*+}_3 */
	di_quark = quarkContract13(quark_propagator * CgmNR,
				   CgmNR * quark_propagator);
	b_prop = trace(S_proj * traceColor(quark_propagator * traceSpin(di_quark)))
	  + 2*trace(S_proj * traceColor(quark_propagator * di_quark));

	/* Multiply by 3 for compatibility with heavy-light routine */
	b_prop *= 3.0;
	break;

      default:
	QDP_error_exit("Unknown baryon: baryons=%d",baryons);
      }

      /* Project on zero momentum: Do a slice-wise sum. */
      multi1d<DComplex> hsum(length);
      hsum = sumMulti(b_prop, timeslice);

      switch (time_rev)
      {
      case 0:
        /* forward */
        for(int t = 0; t < length; ++t)
        {
          int t_eff = (t - t0 + length) % length;

          if ( bc_spec < 0 && (t_eff+t0) >= length)
          {
            barprop[baryons][t_eff] = -0.5 * Complex(hsum[t]);
          }
          else
            barprop[baryons][t_eff] =  0.5 * Complex(hsum[t]);
        }
	break;

      case 1:
        /* backward */
        for(int t = 0; t < length; ++t)
        {
          int t_eff = (length - t + t0) % length;

          if ( bc_spec < 0 && (t_eff-t0) > 0)
          {
            barprop[baryons][t_eff] -=  0.5 * Complex(hsum[t]);
          }
          else
            barprop[baryons][t_eff] +=  0.5 * Complex(hsum[t]);
        }
      }

#if 0
      /* Project onto non-zero momentum if desired */
      if ( num_mom != 0 )
      {
	Complex disp(num_mom, length);
  
	sftmom(b_prop, disp, FftP, num_mom, j_decay);

	for(int m = 0; m < num_mom; ++m)
	{
	  switch (time_rev)
	  {
          case 0:
            /* forward */
            for(int t = 0; t < length; ++t)
            {
              int t_eff = (t - t0 + length) % length;
              if ( bc_spec < 0 && (t_eff+t0) >= length)
              {
                bardisp[baryons][m][t_eff] = -0.5 * disp(m,t);
              }
              else
                bardisp[baryons][m][t_eff] =  0.5 * disp(m,t);
            }
	    break;

          case 1:
            /* backward */
            for(int t = 0; t < length; ++t)
            {
              int t_eff = (length - t + t0) % length;
              if ( bc_spec < 0 && (t_eff-t0) > 0)
              {
                bardisp[baryons][m][t_eff] -=  0.5 * disp[m][t];
              }
              else
                bardisp[baryons][m][t_eff] +=  0.5 * disp[m][t];
            }
	  }
	}
      }
#endif

    } /* end loop over baryons */

    /* Time-charge reverse the quark propagators */
    /* S_{CT} = gamma_5 gamma_4 = gamma_1 gamma_2 gamma_3 = Gamma(7) */
    LatticePropagator q1_tmp = - (Gamma(7) * quark_propagator);
    quark_propagator = q1_tmp * Gamma(7);
  }
}
