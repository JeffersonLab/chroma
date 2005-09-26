// $Id: barseqsrc_w.cc,v 2.1 2005-09-26 04:48:35 edwards Exp $
/*! \file
 *  \brief Construct baryon sequential sources.
 */

#include "chromabase.h"
#include "meas/hadron/barseqsrc_w.h"
#include "meas/hadron/seqsrc_funcmap_w.h"

namespace Chroma 
{

  // Anonymous namespace
  /*! \ingroup hadron */
  namespace
  {
    //! Check only 1 prop passed
    void check1Args(const char* name, const multi1d<LatticePropagator>& quark_propagators)
    {
      if (quark_propagators.size() != 1)
      {
//      QDPIO::cerr << __func__ << ": expect only 1 prop" << endl;
//      QDP_abort(1);

	ostringstream s;
	s << name << ": expecting 1 prop, instead passed = " << quark_propagators.size() << endl;
	throw s.str();
      }
    }

    //! Check only 2 props passed
    void check2Args(const char* name, const multi1d<LatticePropagator>& quark_propagators)
    {
      if (quark_propagators.size() != 2)
      {
//      QDPIO::cerr << __func__ << ": expect only 2 prop" << endl;
//      QDP_abort(1);

	ostringstream s;
	s << name << ": expecting 2 props, instead passed = " << quark_propagators.size() << endl;
	throw s.str();
      }
    }
  }

  //! Construct baryon sequential sources
  /*!
   * \ingroup hadron
   *
   * \param quark_propagators    array of quark propagators ( Read )
   *
   * \return the sequential source before projection onto the sink
   */

  LatticePropagator barNuclUUnpol(const multi1d<LatticePropagator>& quark_propagators)
  {
    START_CODE();

    check2Args(__func__, quark_propagators);

    LatticePropagator src_prop_tmp;
    LatticePropagator q1_tmp;
    LatticePropagator q2_tmp;
    LatticePropagator di_quark;
  
    SpinMatrix g_one = 1.0;

    /* "\bar u O u" insertion in proton, ie. "(u C gamma_5 d) u" */
    /* T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 */
    /* C gamma_5 = Gamma(5) = - (C gamma_5)^T */
    di_quark = quarkContract24(quark_propagators[0] * Gamma(5), 
			       Gamma(5) * quark_propagators[1]);
    src_prop_tmp = (di_quark + Gamma(8)*di_quark) 
      + traceSpin(di_quark)*(g_one + Gamma(8)*g_one);

    q1_tmp = Gamma(5) * quark_propagators[1] * Gamma(5);
    q2_tmp = quark_propagators[0] + quark_propagators[0]*Gamma(8);
    src_prop_tmp -= quarkContract13(q1_tmp, q2_tmp) + transposeSpin(quarkContract12(q2_tmp, q1_tmp));
    src_prop_tmp *= 0.5;
    return src_prop_tmp;

  }

  LatticePropagator barNuclDUnpol(const multi1d<LatticePropagator>& quark_propagators) 
  {
    START_CODE();

    check1Args(__func__, quark_propagators);

    LatticePropagator src_prop_tmp;
    LatticePropagator q1_tmp;
    LatticePropagator q2_tmp;
  
    SpinMatrix g_one = 1.0;

    /* "\bar d O d" insertion in proton, ie. "(u C gamma_5 d) u" */
    /* T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 */
    /* C gamma_5 = Gamma(5) = - (C gamma_5)^T */
    q2_tmp = quark_propagators[0] * Gamma(5);
    q1_tmp = q2_tmp + Gamma(8) * q2_tmp;
    q2_tmp = Gamma(5) * quark_propagators[0];
    src_prop_tmp = -quarkContract14(q1_tmp, q2_tmp);

    q1_tmp = q2_tmp * Gamma(5);
    q2_tmp = quark_propagators[0] + quark_propagators[0] * Gamma(8);
    src_prop_tmp -= transposeSpin(quarkContract12(q2_tmp, q1_tmp));

    src_prop_tmp *= 0.5;
    return src_prop_tmp;
  }


  LatticePropagator barNuclUPol(const multi1d<LatticePropagator>& quark_propagators) 
  {
    START_CODE();

    check2Args(__func__, quark_propagators);

    LatticePropagator src_prop_tmp;
    LatticePropagator q1_tmp;
    LatticePropagator q2_tmp;
    LatticePropagator di_quark;
    LatticeColorMatrix col_mat;
  
    SpinMatrix g_one = 1.0;

    /* "\bar u O u" insertion in proton, ie. "(u C gamma_5 d) u" */
    /* T = \Sigma_3 (1 + gamma_4) / 2 = -i (Gamma(3) + Gamma(11)) / 2 */
    /* C gamma_5 = Gamma(5) = - (C gamma_5)^T */
    q1_tmp = quark_propagators[0] * Gamma(5);
    q2_tmp = Gamma(5) * quark_propagators[1];
    di_quark = quarkContract24(q1_tmp, q2_tmp);
    src_prop_tmp = (Gamma(3) * di_quark  +  Gamma(11) * di_quark);

    col_mat = traceSpin(di_quark);
    src_prop_tmp += col_mat * (Gamma(3) * g_one + Gamma(11) * g_one);

    q1_tmp = q2_tmp * Gamma(5);
    q2_tmp = quark_propagators[0] * Gamma(3);
    q2_tmp += quark_propagators[0] * Gamma(11);
    src_prop_tmp -= quarkContract13(q1_tmp, q2_tmp) + transposeSpin(quarkContract12(q2_tmp, q1_tmp));

    q1_tmp = 0.5 * timesMinusI(src_prop_tmp);
    src_prop_tmp = q1_tmp;
    return src_prop_tmp;
  }


  LatticePropagator barNuclDPol(const multi1d<LatticePropagator>& quark_propagators) 
  {
    START_CODE();

    check1Args(__func__, quark_propagators);

    LatticePropagator src_prop_tmp;
    LatticePropagator q1_tmp;
    LatticePropagator q2_tmp;
  
    SpinMatrix g_one = 1.0;

    /* "\bar d O d" insertion in proton, ie. "(u C gamma_5 d) u" */
    /* T = \Sigma_3 (1 + gamma_4) / 2 = -i (Gamma(3) + Gamma(11)) / 2 */

    /* C gamma_5 = Gamma(5) = - (C gamma_5)^T */
    q2_tmp = quark_propagators[0] * Gamma(5);
    q1_tmp = Gamma(3) * q2_tmp + Gamma(11) * q2_tmp;
    q2_tmp = Gamma(5) * quark_propagators[0];
    src_prop_tmp = -quarkContract14(q1_tmp, q2_tmp);

    q1_tmp = q2_tmp * Gamma(5);
    q2_tmp = quark_propagators[0] * Gamma(3) + quark_propagators[0] * Gamma(11);
    src_prop_tmp -= transposeSpin(quarkContract12(q2_tmp, q1_tmp));

    q1_tmp = 0.5 * timesMinusI(src_prop_tmp);
    src_prop_tmp = q1_tmp;
    return src_prop_tmp;
  }


  LatticePropagator barDeltaUUnpol(const multi1d<LatticePropagator>& quark_propagators) 
  {
    START_CODE();

    check2Args(__func__, quark_propagators);

    LatticePropagator src_prop_tmp;
    LatticePropagator q1_tmp;
    LatticePropagator q2_tmp;
    LatticePropagator di_quark;
    LatticeColorMatrix col_mat;
  
    SpinMatrix g_one = 1.0;

    /* C = Gamma(10) */
    SpinMatrix Cgm = 0.5 * (Gamma(10) * (Gamma(2) * g_one  +  timesI(Gamma(1) * g_one)));

    /* "\bar u O u" insertion in Delta^+,
       ie. "2*(u C gamma_- d) u + (u C gamma_- u) d" */
    /* T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 */

    /* C gamma_- = Cgm = (C gamma_-)^T */
    q1_tmp = quark_propagators[0] * Cgm;
    q2_tmp = Cgm * quark_propagators[1];
    di_quark = quarkContract24(q1_tmp, q2_tmp);
    src_prop_tmp = di_quark + Gamma(8) * di_quark;

    col_mat = traceSpin(di_quark);
    q1_tmp = 1;
    di_quark = q1_tmp + Gamma(8) * q1_tmp;
    src_prop_tmp += col_mat * di_quark;

    q1_tmp = q2_tmp * Cgm;
    q2_tmp = quark_propagators[0] + quark_propagators[0] * Gamma(8);
    src_prop_tmp += quarkContract13(q1_tmp, q2_tmp) + transposeSpin(quarkContract12(q2_tmp, q1_tmp));

    q1_tmp = Cgm * quark_propagators[0];
    q2_tmp = q1_tmp * Cgm;
    q1_tmp = quark_propagators[1] + Gamma(8) * quark_propagators[1];
    src_prop_tmp += transposeSpin(quarkContract12(q1_tmp, q2_tmp)) + quarkContract24(q1_tmp, q2_tmp);

    q2_tmp = q1_tmp * Cgm;
    q1_tmp = Cgm * quark_propagators[0];
    src_prop_tmp += quarkContract14(q2_tmp, q1_tmp);

    q1_tmp = Cgm * quark_propagators[1];
    q2_tmp = q1_tmp + q1_tmp * Gamma(8);
    q1_tmp = quark_propagators[0] * Cgm;
    src_prop_tmp += quarkContract14(q1_tmp, q2_tmp) + quarkContract13(q1_tmp, q2_tmp);
    src_prop_tmp *= 2;

    return src_prop_tmp;
  }


  LatticePropagator barDeltaDUnpol(const multi1d<LatticePropagator>& quark_propagators) 
  {
    START_CODE();

    check1Args(__func__, quark_propagators);

    LatticePropagator src_prop_tmp;
    LatticePropagator q1_tmp;
    LatticePropagator q2_tmp;
    LatticePropagator di_quark;
    LatticeColorMatrix col_mat;
  
    SpinMatrix g_one = 1.0;

    /* C = Gamma(10) */
    SpinMatrix Cgm = 0.5 * (Gamma(10) * (Gamma(2) * g_one  +  timesI(Gamma(1) * g_one)));

    /* "\bar d O d" insertion in Delta^+,
       ie. "2*(u C gamma_- d) u + (u C gamma_- u) d" */
    /* T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 */

    /* C gamma_- = Cgm = (C gamma_-)^T */
    q2_tmp = quark_propagators[0] * Cgm;
    q1_tmp = q2_tmp + Gamma(8) * q2_tmp;
    q2_tmp = Cgm * quark_propagators[0];
    src_prop_tmp = quarkContract14(q1_tmp, q2_tmp);

    q1_tmp = q2_tmp * Cgm;
    q2_tmp = quark_propagators[0] + quark_propagators[0] * Gamma(8);
    src_prop_tmp += transposeSpin(quarkContract12(q2_tmp, q1_tmp));

    q1_tmp = quark_propagators[0] * Cgm;
    q2_tmp = Cgm * quark_propagators[0];
    di_quark = quarkContract24(q1_tmp, q2_tmp);
    src_prop_tmp += di_quark + Gamma(8) * di_quark;

    di_quark = quarkContract13(q1_tmp, q2_tmp);
    src_prop_tmp += di_quark + di_quark * Gamma(8);
    src_prop_tmp *= 2;

    col_mat = traceSpin(di_quark);
    q1_tmp = 1;
    q2_tmp = q1_tmp + Gamma(8) * q1_tmp;
    src_prop_tmp += col_mat * q2_tmp;

    return src_prop_tmp;
  }


  LatticePropagator barNuclUUnpolNR(const multi1d<LatticePropagator>& quark_propagators) 
  {
    START_CODE();

    check2Args(__func__, quark_propagators);

    LatticePropagator src_prop_tmp;
    LatticePropagator q1_tmp;
    LatticePropagator q2_tmp;
    LatticePropagator di_quark;
    LatticeColorMatrix col_mat;
  
    SpinMatrix g_one = 1.0;

    /* C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 ) */ 
    SpinMatrix Cg5NR = 0.5 * (Gamma(5) * (g_one + (g_one * Gamma(8))));


    /* "\bar u O u" insertion in NR proton, ie. 
     * "(u C gamma_5 (1/2)(1 + gamma_4)  d) u" */
    /* T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 */

    /* Use precomputed Cg5NR = C gamma_5 (1/2) ( 1 + g_4 ) */
    q1_tmp = quark_propagators[0] * Cg5NR;
    q2_tmp = Cg5NR * quark_propagators[1];
    di_quark = quarkContract24(q1_tmp, q2_tmp);
    src_prop_tmp = di_quark + Gamma(8) * di_quark;

    col_mat = traceSpin(di_quark);
    q1_tmp = 1;
    di_quark = q1_tmp + Gamma(8) * q1_tmp;
    src_prop_tmp += col_mat * di_quark;

    q1_tmp = q2_tmp * Cg5NR;
    q2_tmp = quark_propagators[0] + quark_propagators[0] * Gamma(8);
    src_prop_tmp -= quarkContract13(q1_tmp, q2_tmp) + transposeSpin(quarkContract12(q2_tmp, q1_tmp));
    src_prop_tmp *= 0.5;

    return src_prop_tmp;
  }


  LatticePropagator barNuclDUnpolNR(const multi1d<LatticePropagator>& quark_propagators) 
  {
    START_CODE();

    check1Args(__func__, quark_propagators);

    LatticePropagator src_prop_tmp;
    LatticePropagator q1_tmp;
    LatticePropagator q2_tmp;
  
    SpinMatrix g_one = 1.0;

    /* C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 ) */ 
    SpinMatrix Cg5NR = 0.5 * (Gamma(5) * (g_one + (g_one * Gamma(8))));


    /* "\bar d O d" insertion in NR proton, ie. 
     * "(u C gamma_5 (1/2)(1 + gamma_4) d) u" */
    /* T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 */

    /* C gamma_5 (1/2)( 1 + g_4)  = Cg5NR  */
    q2_tmp = quark_propagators[0] * Cg5NR;
    q1_tmp = q2_tmp + Gamma(8) * q2_tmp;
    q2_tmp = Cg5NR * quark_propagators[0];
    src_prop_tmp = -quarkContract14(q1_tmp, q2_tmp);

    q1_tmp = q2_tmp * Cg5NR;
    q2_tmp = quark_propagators[0] + quark_propagators[0] * Gamma(8);
    src_prop_tmp -= transposeSpin(quarkContract12(q2_tmp, q1_tmp));
    src_prop_tmp *= 0.5;

    return src_prop_tmp;
  }


  LatticePropagator barNuclUPolNR(const multi1d<LatticePropagator>& quark_propagators) 
  {
    START_CODE();

    check2Args(__func__, quark_propagators);

    LatticePropagator src_prop_tmp;
    LatticePropagator q1_tmp;
    LatticePropagator q2_tmp;
    LatticePropagator di_quark;
    LatticeColorMatrix col_mat;
  
    SpinMatrix g_one = 1.0;

    /* C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 ) */ 
    SpinMatrix Cg5NR = 0.5 * (Gamma(5) * (g_one + (g_one * Gamma(8))));

    /* "\bar u O u" insertion in NR proton, ie. 
     * "(u C gamma_5 (1/2)(1+gamma_4) d) u" */
    /* T = \Sigma_3 (1 + gamma_4) / 2 = -i (Gamma(3) + Gamma(11)) / 2 */

    /* C gamma_5 (1/2) (1 + gamma_4) = Cg5NR */
    q1_tmp = quark_propagators[0] * Cg5NR;
    q2_tmp = Cg5NR * quark_propagators[1];
    di_quark = quarkContract24(q1_tmp, q2_tmp);
    src_prop_tmp = Gamma(3) * di_quark + Gamma(11) * di_quark;

    col_mat = traceSpin(di_quark);
    q1_tmp = 1;
    di_quark = Gamma(3) * q1_tmp  +  Gamma(11) * q1_tmp;
    src_prop_tmp += col_mat * di_quark;

    q1_tmp = q2_tmp * Cg5NR;
    q2_tmp = quark_propagators[0] * Gamma(3)  +  quark_propagators[0] * Gamma(11);
    src_prop_tmp -= quarkContract13(q1_tmp, q2_tmp) + transposeSpin(quarkContract12(q2_tmp, q1_tmp));

    q1_tmp = 0.5 * timesMinusI(src_prop_tmp);
    src_prop_tmp = q1_tmp;

    return src_prop_tmp;
  }


  LatticePropagator barNuclDPolNR(const multi1d<LatticePropagator>& quark_propagators) 
  {
    START_CODE();

    check1Args(__func__, quark_propagators);

    LatticePropagator src_prop_tmp;
    LatticePropagator q1_tmp;
    LatticePropagator q2_tmp;
  
    SpinMatrix g_one = 1.0;

    /* C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 ) */ 
    SpinMatrix Cg5NR = 0.5 * (Gamma(5) * (g_one + (g_one * Gamma(8))));


    /* "\bar d O d" insertion in proton, ie. "(u C gamma_5 (1/2)(1+gamma_4)d) u" */
    /* T = \Sigma_3 (1 + gamma_4) / 2 = -i (Gamma(3) + Gamma(11)) / 2 */

    /* C gamma_5 (1/2)(1+gamma_4)= Cg5NR */
    q2_tmp = quark_propagators[0] * Cg5NR;
    q1_tmp = Gamma(3) * q2_tmp  +  Gamma(11) * q2_tmp;
    q2_tmp = Cg5NR * quark_propagators[0];
    src_prop_tmp = -quarkContract14(q1_tmp, q2_tmp);

    q1_tmp = q2_tmp * Cg5NR;
    q2_tmp = quark_propagators[0] * Gamma(3)  +  quark_propagators[0] * Gamma(11);
    src_prop_tmp -= transposeSpin(quarkContract12(q2_tmp, q1_tmp));

    q1_tmp = 0.5 * timesMinusI(src_prop_tmp);
    src_prop_tmp = q1_tmp;
      
    return src_prop_tmp;
  }


  LatticePropagator barNuclUMixedNR(const multi1d<LatticePropagator>& quark_propagators) 
  {
    START_CODE();

    check2Args(__func__, quark_propagators);

    LatticePropagator src_prop_tmp;
    LatticePropagator q1_tmp;
    LatticePropagator q2_tmp;
    LatticePropagator di_quark;
    LatticeColorMatrix col_mat;
  
    SpinMatrix g_one = 1.0;

    /* C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 ) */ 
    SpinMatrix Cg5NR = 0.5 * (Gamma(5) * (g_one + (g_one * Gamma(8))));


    /* "\bar u O u" insertion in NR proton, ie. 
     * "(u C gamma_5 (1/2)(1 + gamma_4)  d) u" */
    /* T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
       = (1 + Gamma(8) - i G(3) - i G(11)) / 2 */

    /*
     *  Note that this is precisely src NUCL_U_UNPOL_NR + src NUCL_U_POL_NR, and corresponds
     *  to the combination of propagators used by the MIT group
     */


    /* Use precomputed Cg5NR = C gamma_5 (1/2) ( 1 + g_4 ) */
    q1_tmp = quark_propagators[0] * Cg5NR;
    q2_tmp = Cg5NR * quark_propagators[1];
    di_quark = quarkContract24(q1_tmp, q2_tmp);

    /*
     *  We begin with the polarized piece
     */

    q1_tmp = Gamma(3) * di_quark  +  Gamma(11) * di_quark;
    /*
     *  Now multiply by -i
     */

    src_prop_tmp = timesMinusI(q1_tmp);

    /*
     * Now add the unpolarized piece
     */

    src_prop_tmp += di_quark  +  Gamma(8) * di_quark;

    /*
     *  Now the second term
     */

    col_mat = traceSpin(di_quark);
    q1_tmp = 1;

    /*
     *  First the polarized part
     */

    di_quark = Gamma(3) * q1_tmp  +  Gamma(11) * q1_tmp;
    q1_tmp = timesMinusI(di_quark); /*  Multiply by -i */
    src_prop_tmp += col_mat * q1_tmp;

    /*
     *  Now the unpolarized part
     */

    q1_tmp = 1;
    di_quark = q1_tmp  +  Gamma(8) * q1_tmp;
    src_prop_tmp += col_mat * di_quark;

    /*
     *  The third term...
     */

    q1_tmp = q2_tmp * Cg5NR;

    /*  First the polarized piece */

    di_quark = quark_propagators[0] * Gamma(3)  +  quark_propagators[0] * Gamma(11);
    q2_tmp = timesMinusI(di_quark);	/* Muliply by -i */

    /* Add the unpolarized piece */

    q2_tmp += quark_propagators[0]  +  quark_propagators[0] * Gamma(8);

    src_prop_tmp -= quarkContract13(q1_tmp, q2_tmp) + transposeSpin(quarkContract12(q2_tmp, q1_tmp));

    /*  Finally, multiply everything by 1/2 */

    src_prop_tmp *= 0.5;

    return src_prop_tmp;
  }


  LatticePropagator barNuclDMixedNR(const multi1d<LatticePropagator>& quark_propagators) 
  {
    START_CODE();

    check1Args(__func__, quark_propagators);

    LatticePropagator src_prop_tmp;
    LatticePropagator q1_tmp;
    LatticePropagator q2_tmp;
    LatticePropagator di_quark;
  
    SpinMatrix g_one = 1.0;

    /* C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 ) */ 
    SpinMatrix Cg5NR = 0.5 * (Gamma(5) * (g_one + (g_one * Gamma(8))));


    /* "\bar d O d" insertion in proton, ie. "(u C gamma_5 (1/2)(1+gamma_4)d) u" */
    /* T = \Sigma_3 (1 + gamma_4) / 2 = -i (Gamma(3) + Gamma(11)) / 2 */

    /* C gamma_5 (1/2)(1+gamma_4)= Cg5NR */
    q2_tmp = quark_propagators[0] * Cg5NR;

    /*
     *  First the polarized piece
     */

    di_quark = Gamma(3) * q2_tmp   +  Gamma(11) * q2_tmp;
    q1_tmp = timesMinusI(di_quark);

    /*
     *  Now add the unpolarized piece
     */

    q1_tmp += q2_tmp  +  Gamma(8) * q2_tmp;


    q2_tmp = Cg5NR * quark_propagators[0];
    src_prop_tmp = -quarkContract14(q1_tmp, q2_tmp);

    /* Now the second term */

    q1_tmp = q2_tmp * Cg5NR;

    /*
     *  First the polarized piece
     */

    di_quark = quark_propagators[0] * Gamma(3) + quark_propagators[0] * Gamma(11);
    q2_tmp = timesMinusI(di_quark);

    /*
     *  Now add the unpolarized piece
     */

    q2_tmp += quark_propagators[0] + quark_propagators[0] * Gamma(8);

    src_prop_tmp -= transposeSpin(quarkContract12(q2_tmp, q1_tmp));
    src_prop_tmp *= 0.5;

    return src_prop_tmp;
  }



  // Patch for the quarkContract12 piece in NuclUMixedNR and NuclDMixedNR
  LatticePropagator barNuclPatchMixedNR(const multi1d<LatticePropagator>& quark_propagators) 
  {
    START_CODE();

    check1Args(__func__, quark_propagators);

    LatticePropagator src_prop_tmp;
    LatticePropagator q1_tmp;
    LatticePropagator q2_tmp;
    LatticePropagator di_quark;
  
    SpinMatrix g_one = 1.0;

    /* C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 ) */ 
    SpinMatrix Cg5NR = 0.5 * (Gamma(5) * (g_one + (g_one * Gamma(8))));

    /* "\bar d O d" insertion in proton, ie. "(u C gamma_5 (1/2)(1+gamma_4)d) u" */
    /* T = \Sigma_3 (1 + gamma_4) / 2 = -i (Gamma(3) + Gamma(11)) / 2 */

    /* C gamma_5 (1/2)(1+gamma_4)= Cg5NR */

    q2_tmp = Cg5NR * quark_propagators[0];
    q1_tmp = q2_tmp * Cg5NR;

    di_quark = quark_propagators[0] * Gamma(3) + quark_propagators[0] * Gamma(11);
    q2_tmp = timesMinusI(di_quark);

    q2_tmp += quark_propagators[0] + quark_propagators[0] * Gamma(8);

    di_quark = quarkContract12(q2_tmp, q1_tmp);
    src_prop_tmp = di_quark - transposeSpin(di_quark);   // bad guy - good guy
    src_prop_tmp *= Real(0.5);

    return src_prop_tmp;
  }



  //! Baryon sequential sources
  /*! \ingroup hadron */
  namespace BaryonSeqSourceCallMapEnv
  { 
    bool registerAll(void) 
    {
      bool success = true;
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("NUCL_U_UNPOL"), 
								  barNuclUUnpol);

      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("NUCL_D_UNPOL"), 
								  barNuclDUnpol);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("NUCL_U_POL"),
								  barNuclUPol);

      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("NUCL_D_POL"),
								  barNuclDPol);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("DELTA_U_UNPOL"),
								  barDeltaUUnpol);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("DELTA_D_UNPOL"),
								  barDeltaDUnpol);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("NUCL_U_UNPOL_NONREL"),
								  barNuclUUnpolNR);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("NUCL_D_UNPOL_NONREL"),
								  barNuclDUnpolNR);

      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("NUCL_U_POL_NONREL"),
								  barNuclUPolNR);

      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("NUCL_D_POL_NONREL"),
								  barNuclDPolNR);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("NUCL_U_MIXED_NONREL"),
								  barNuclUMixedNR);

      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("NUCL_D_MIXED_NONREL"),   
								  barNuclDMixedNR);

      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("NUCL_PATCH_MIXED_NONREL"),   
								  barNuclPatchMixedNR);
      
      return success;
    }

    bool registered = registerAll();
  } // namespace BaryonSeqSourceCallMapEnv


}  // end namespace Chroma
