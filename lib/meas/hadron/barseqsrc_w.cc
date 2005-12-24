// $Id: barseqsrc_w.cc,v 2.6 2005-12-24 21:20:14 edwards Exp $
/*! \file
 *  \brief Construct baryon sequential sources.
 */

#include "chromabase.h"
#include "meas/hadron/barseqsrc_w.h"
#include "meas/hadron/seqsrc_funcmap_w.h"
#include "meas/hadron/barspinmat_w.h"

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


  //! Nucleon-Nucleon U piece with general projector and Cg5
  LatticePropagator barNuclUTCg5(const multi1d<LatticePropagator>& quark_propagators,
				 const SpinMatrix& T, const SpinMatrix& Cg5) 
  {
    START_CODE();

    check2Args(__func__, quark_propagators);

    LatticePropagator src_prop_tmp;
    LatticePropagator q1_tmp;
    LatticePropagator q2_tmp;
    LatticePropagator di_quark;
    LatticeColorMatrix col_mat;
  
    /* "\bar u O u" insertion in NR proton, ie. 
     * "(u Cg5 d) u" */
    /* Some generic T */

    // Use precomputed Cg5
    q1_tmp = quark_propagators[0] * Cg5;
    q2_tmp = Cg5 * quark_propagators[1];
    di_quark = quarkContract24(q1_tmp, q2_tmp);

    // First term
    src_prop_tmp = T * di_quark;

    // Now the second term
    src_prop_tmp += traceSpin(di_quark) * T;

    // The third term...
    q1_tmp = q2_tmp * Cg5;
    q2_tmp = quark_propagators[0] * T;

    src_prop_tmp -= quarkContract13(q1_tmp, q2_tmp) + transposeSpin(quarkContract12(q2_tmp, q1_tmp));

    return src_prop_tmp;
  }


  //! Nucleon-Nucleon D piece with general projector and Cg5
  LatticePropagator barNuclDTCg5(const multi1d<LatticePropagator>& quark_propagators,
				 const SpinMatrix& T, const SpinMatrix& Cg5) 
  {
    START_CODE();

    check1Args(__func__, quark_propagators);

    LatticePropagator src_prop_tmp;
    LatticePropagator q1_tmp;
    LatticePropagator q2_tmp;

    /* "\bar d O d" insertion in NR proton, ie. 
     * "(u Cg5 d) u" */
    /* Some generic T */

    // First term
    q2_tmp = quark_propagators[0] * Cg5;
    q1_tmp = T * q2_tmp;

    q2_tmp = Cg5 * quark_propagators[0];
    src_prop_tmp = -quarkContract14(q1_tmp, q2_tmp);

    // Second term
    q1_tmp = q2_tmp * Cg5;
    q2_tmp = quark_propagators[0] * T;

    src_prop_tmp -= transposeSpin(quarkContract12(q2_tmp, q1_tmp));

    return src_prop_tmp;
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

    /* "\bar u O u" insertion in proton, ie. 
     * "(u C gamma_5 d) u" */

    /* C gamma_5 = Gamma(5) = - (C gamma_5)^T */
    /* T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 */

    // Compute the  \bar{u} O u  insertion
    return barNuclUTCg5(quark_propagators, BaryonSpinMats::Tunpol(), BaryonSpinMats::Cg5());
  }


  LatticePropagator barNuclDUnpol(const multi1d<LatticePropagator>& quark_propagators) 
  {
    START_CODE();

    check1Args(__func__, quark_propagators);

    /* "\bar d O d" insertion in proton, ie. 
     * "(u C gamma_5 d) u" */

    /* C gamma_5 = Gamma(5) = - (C gamma_5)^T */
    /* T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 */

    // Compute the  \bar{d} O d  insertion
    return barNuclDTCg5(quark_propagators, BaryonSpinMats::Tunpol(), BaryonSpinMats::Cg5());
  }


  LatticePropagator barNuclUPol(const multi1d<LatticePropagator>& quark_propagators) 
  {
    START_CODE();

    check2Args(__func__, quark_propagators);

    /* "\bar u O u" insertion in proton, ie. 
     * "(u C gamma_5 (1/2)(1 + gamma_4)  d) u" */

    // C g_5 = C gamma_5
    /* T = \Sigma_3 (1 + gamma_4) / 2 = -i (Gamma(3) + Gamma(11)) / 2 */

    // Compute the  \bar{u} O u  insertion
    return barNuclUTCg5(quark_propagators, BaryonSpinMats::Tunpol(), BaryonSpinMats::Cg5());
  }


  LatticePropagator barNuclDPol(const multi1d<LatticePropagator>& quark_propagators) 
  {
    START_CODE();

    check1Args(__func__, quark_propagators);

    /* "\bar u O u" insertion in proton, ie. 
     * "(u C gamma_5 (1/2)(1 + gamma_4)  d) u" */

    // C g_5 = C gamma_5
    /* T = \Sigma_3 (1 + gamma_4) / 2 = -i (Gamma(3) + Gamma(11)) / 2 */

    // Compute the  \bar{u} O u  insertion
    return barNuclDTCg5(quark_propagators, BaryonSpinMats::Tpol(), BaryonSpinMats::Cg5());
  }

  LatticePropagator barNuclUUnpolNR(const multi1d<LatticePropagator>& quark_propagators) 
  {
    START_CODE();

    check2Args(__func__, quark_propagators);

    /* "\bar u O u" insertion in NR proton, ie. 
     * "(u C gamma_5 (1/2)(1 + gamma_4)  d) u" */

    // C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 )
    /* T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 */

    // Compute the  \bar{u} O u  insertion
    return barNuclUTCg5(quark_propagators, BaryonSpinMats::Tunpol(), BaryonSpinMats::Cg5NR());
  }


  LatticePropagator barNuclDUnpolNR(const multi1d<LatticePropagator>& quark_propagators) 
  {
    START_CODE();

    check1Args(__func__, quark_propagators);

    /* "\bar d O d" insertion in NR proton, ie. 
     * "(u C gamma_5 (1/2)(1 + gamma_4)  d) u" */

    // C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 )
    /* T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 */

    // Compute the  \bar{d} O d  insertion
    return barNuclDTCg5(quark_propagators, BaryonSpinMats::Tunpol(), BaryonSpinMats::Cg5NR());
  }


  LatticePropagator barNuclUPolNR(const multi1d<LatticePropagator>& quark_propagators) 
  {
    START_CODE();

    check2Args(__func__, quark_propagators);

    /* "\bar u O u" insertion in NR proton, ie. 
     * "(u C gamma_5 (1/2)(1 + gamma_4)  d) u" */

    // C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 )
    /* T = \Sigma_3 (1 + gamma_4) / 2 = -i (Gamma(3) + Gamma(11)) / 2 */

    // Compute the  \bar{u} O u  insertion
    return barNuclUTCg5(quark_propagators, BaryonSpinMats::Tpol(), BaryonSpinMats::Cg5NR());
  }


  LatticePropagator barNuclDPolNR(const multi1d<LatticePropagator>& quark_propagators) 
  {
    START_CODE();

    check1Args(__func__, quark_propagators);

    /* "\bar d O d" insertion in NR proton, ie. 
     * "(u C gamma_5 (1/2)(1 + gamma_4)  d) u" */

    // C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 )
    /* T = \Sigma_3 (1 + gamma_4) / 2 = -i (Gamma(3) + Gamma(11)) / 2 */

    // Compute the  \bar{d} O d  insertion
    return barNuclDTCg5(quark_propagators, BaryonSpinMats::Tpol(), BaryonSpinMats::Cg5NR());
  }


  //! \bar u O u" insertion in NR proton
  /*!
   * "\bar u O u" insertion in NR proton, ie. 
   * "(u C gamma_5 (1/2)(1 + gamma_4)  d) u" 
   * 
   * $C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 )$
   * 
   * $T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
   *   = (1 + Gamma(8) - i G(3) - i G(11)) / 2$
   */
  LatticePropagator barNuclUMixedNR(const multi1d<LatticePropagator>& quark_propagators) 
  {
    START_CODE();

    check2Args(__func__, quark_propagators);

    /* "\bar u O u" insertion in NR proton, ie. 
     * "(u C gamma_5 (1/2)(1 + gamma_4)  d) u" */

    // C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 )

    // T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
    //   = (1 + Gamma(8) - i G(3) - i G(11)) / 2

    // Compute the  \bar{u} O u  insertion
    return barNuclUTCg5(quark_propagators, BaryonSpinMats::Tmixed(), BaryonSpinMats::Cg5NR());
  }


  LatticePropagator barNuclDMixedNR(const multi1d<LatticePropagator>& quark_propagators) 
  {
    START_CODE();

    check1Args(__func__, quark_propagators);

    /* "\bar d O d" insertion in NR proton, ie. 
     * "(u C gamma_5 (1/2)(1 + gamma_4)  d) u" */

    // C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 )

    // T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
    //   = (1 + Gamma(8) - i G(3) - i G(11)) / 2

    // Compute the  \bar{d} O d  insertion
    return barNuclDTCg5(quark_propagators, BaryonSpinMats::Tmixed(), BaryonSpinMats::Cg5NR());
  }

  //! \bar u O u" insertion in NR proton
  /*!
   * "\bar u O u" insertion in NR proton, ie. 
   * "(u C gamma_5 (1/2)(1 - gamma_4)  d) u" 
   * 
   * $C g_5 NR = (1/2)*C gamma_5 * ( 1 - g_4 )$
   * 
   * $T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
   *   = (1 + Gamma(8) - i G(3) - i G(11)) / 2$
   */
  LatticePropagator barNuclUMixedNRnegPar(const multi1d<LatticePropagator>& quark_propagators) 
  {
    START_CODE();

    check2Args(__func__, quark_propagators);

    /* "\bar u O u" insertion in NR proton, ie. 
     * "(u C gamma_5 (1/2)(1 - gamma_4)  d) u" */

    // C g_5 NR = (1/2)*C gamma_5 * ( 1 - g_4 )

    // T = (1 + \Sigma_3)*(1 - gamma_4) / 2 
    //   = (1 - Gamma(8) + i G(3) - i G(11)) / 2

    // Compute the  \bar{u} O u  insertion
    return barNuclUTCg5(quark_propagators, 
			BaryonSpinMats::TmixedNegPar(), BaryonSpinMats::Cg5NRnegPar());
  }

 LatticePropagator barNuclDMixedNRnegPar(const multi1d<LatticePropagator>& quark_propagators) 
  {
    START_CODE();

    check1Args(__func__, quark_propagators);

    /* "\bar d O d" insertion in NR proton, ie. 
     * "(u C gamma_5 (1/2)(1 - gamma_4)  d) u" */

    // C g_5 NR = (1/2)*C gamma_5 * ( 1 - g_4 )

    // T = (1 + \Sigma_3)*(1 - gamma_4) / 2 
    //   = (1 - Gamma(8) + i G(3) - i G(11)) / 2

    // Compute the  \bar{d} O d  insertion
    return barNuclDTCg5(quark_propagators, 
			BaryonSpinMats::TmixedNegPar(), BaryonSpinMats::Cg5NRnegPar());
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
  
    /* C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 ) */ 
    SpinMatrix Cg5NR = BaryonSpinMats::Cg5NR();

    /* "\bar d O d" insertion in proton, ie. "(u C gamma_5 (1/2)(1+gamma_4)d) u" */
    // T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
    //   = (1 + Gamma(8) - i G(3) - i G(11)) / 2

    /* C gamma_5 (1/2)(1+gamma_4)= Cg5NR */

    q2_tmp = Cg5NR * quark_propagators[0];
    q1_tmp = q2_tmp * Cg5NR;

    q2_tmp = quark_propagators[0] * BaryonSpinMats::Tmixed();

    di_quark = quarkContract12(q2_tmp, q1_tmp);
    src_prop_tmp = di_quark - transposeSpin(di_quark);   // bad guy - good guy

    return src_prop_tmp;
  }


  /** The octet baryon sequantial sources **/
  //! \bar d O d" insertion in NR proton
  /*!
   * "\bar d O d" insertion in NR proton, ie. 
   * "(s C gamma_5 (1/2)(1 + gamma_4)  d) s" 
   * 
   * $C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 )$
   * 
   * $T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
   *   = (1 + Gamma(8) - i G(3) - i G(11)) / 2$
   
   * The d-quark insertion for a Xi baryon: primarily for the
   * Xi- to Xi0 transision 

  */
  LatticePropagator barXiDMixedNR(const multi1d<LatticePropagator>& quark_propagators) 
  {
    START_CODE();
    
    check1Args(__func__, quark_propagators);

    /* "\bar u O u" insertion in NR proton, ie. 
     * "(u C gamma_5 (1/2)(1 + gamma_4)  d) u" */

    // C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 )

    // T = (1 + \Sigma_3)*(1 + gamma_4) / 2 
    //   = (1 + Gamma(8) - i G(3) - i G(11)) / 2

    END_CODE();

    // Compute the  \bar{d} O d  insertion
    // The Xi is just like the proton with up quark replaced with the strange
    // the single quark propagator passed in is just the strange quark propagator
    return barNuclDTCg5(quark_propagators, BaryonSpinMats::Tmixed(), BaryonSpinMats::Cg5NR());
  }
  


  //! "\bar u O u" insertion in Delta^+
  LatticePropagator barDeltaUTsp(const multi1d<LatticePropagator>& quark_propagators,
				 const SpinMatrix& T, const SpinMatrix& sp) 
  {
    START_CODE();

    check2Args(__func__, quark_propagators);

    LatticePropagator src_prop_tmp;
    LatticePropagator q1_tmp;
    LatticePropagator q2_tmp;
    LatticePropagator di_quark;
    LatticeColorMatrix col_mat;
  
    /* "\bar u O u" insertion in Delta^+,
       ie. "2*(u sp d) u + (u sp u) d" */
    // generic T

    q1_tmp = quark_propagators[0] * sp;
    q2_tmp = sp * quark_propagators[1];
    di_quark = quarkContract24(q1_tmp, q2_tmp);
    src_prop_tmp = T * di_quark;

    src_prop_tmp += traceSpin(di_quark) * T;

    q1_tmp = q2_tmp * sp;
    q2_tmp = quark_propagators[0] * T;
    src_prop_tmp += quarkContract13(q1_tmp, q2_tmp) + transposeSpin(quarkContract12(q2_tmp, q1_tmp));

    q1_tmp = sp * quark_propagators[0];
    q2_tmp = q1_tmp * sp;
    q1_tmp = T * quark_propagators[1];
    src_prop_tmp += transposeSpin(quarkContract12(q1_tmp, q2_tmp)) + quarkContract24(q1_tmp, q2_tmp);

    q2_tmp = q1_tmp * sp;
    q1_tmp = sp * quark_propagators[0];
    src_prop_tmp += quarkContract14(q2_tmp, q1_tmp);

    q1_tmp = sp * quark_propagators[1];
    q2_tmp = q1_tmp * T;
    q1_tmp = quark_propagators[0] * sp;
    src_prop_tmp += quarkContract14(q1_tmp, q2_tmp) + quarkContract13(q1_tmp, q2_tmp);
    src_prop_tmp *= 2;

    END_CODE();

    return src_prop_tmp;
  }


  //! "\bar d O d" insertion in Delta^+
  LatticePropagator barDeltaDTsp(const multi1d<LatticePropagator>& quark_propagators,
				 const SpinMatrix& T, const SpinMatrix& sp) 
  {
    START_CODE();

    check1Args(__func__, quark_propagators);

    LatticePropagator src_prop_tmp;
    LatticePropagator q1_tmp;
    LatticePropagator q2_tmp;
    LatticePropagator di_quark;

    /* "\bar d O d" insertion in Delta^+,
       ie. "2*(u sp d) u + (u sp u) d" */
    // generic T

    q2_tmp = quark_propagators[0] * sp;
    q1_tmp = T * q2_tmp;
    q2_tmp = sp * quark_propagators[0];
    src_prop_tmp = quarkContract14(q1_tmp, q2_tmp);

    q1_tmp = q2_tmp * sp;
    q2_tmp = quark_propagators[0] * T;
    src_prop_tmp += transposeSpin(quarkContract12(q2_tmp, q1_tmp));

    q1_tmp = quark_propagators[0] * sp;
    q2_tmp = sp * quark_propagators[0];
    src_prop_tmp += T * quarkContract24(q1_tmp, q2_tmp);

    di_quark = quarkContract13(q1_tmp, q2_tmp);
    src_prop_tmp += di_quark * T;
    src_prop_tmp *= 2;

    src_prop_tmp += traceSpin(di_quark) * T;

    END_CODE();

    return src_prop_tmp;
  }


  //! "\bar u O u" insertion in Delta^+
  LatticePropagator barDeltaUUnpol(const multi1d<LatticePropagator>& quark_propagators)
  {
    check2Args(__func__, quark_propagators);

    /* "\bar u O u" insertion in Delta^+,
       ie. "2*(u sp d) u + (u sp u) d" */
    /* C gamma_- = sp = (C gamma_-)^T */
    /* T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 */

    // Agghh, we have a goofy factor of 2 normalization factor here. The
    // ancient szin way didn't care about norms, so it happily made it
    // 2 times too big. There is a missing 0.5 in T_unpol guy.
    // Since nobody has used this code before, we are switching to a more
    // sane convention and breaking agreement with the old szin code.
    return barDeltaUTsp(quark_propagators, BaryonSpinMats::Tunpol(), BaryonSpinMats::Cgm());
  }


  //! "\bar d O d" insertion in Delta^+
  LatticePropagator barDeltaDUnpol(const multi1d<LatticePropagator>& quark_propagators) 
  {
    START_CODE();

    check1Args(__func__, quark_propagators);

    /* "\bar d O d" insertion in Delta^+,
       ie. "2*(u sp d) u + (u sp u) d" */
    /* C gamma_- = sp = (C gamma_-)^T */
    /* T = (1 + gamma_4) / 2 = (1 + Gamma(8)) / 2 */

    // Agghh, we have a goofy factor of 2 normalization factor here. The
    // ancient szin way didn't care about norms, so it happily made it
    // 2 times too big. There is a missing 0.5 in T_unpol guy.
    // Since nobody has used this code before, we are switching to a more
    // sane convention and breaking agreement with the old szin code.
    return barDeltaDTsp(quark_propagators, BaryonSpinMats::Tunpol(), BaryonSpinMats::Cgm());
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

      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("NUCL_U_MIXED_NONREL_NEGPAR"),
								  barNuclUMixedNRnegPar);

      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("NUCL_D_MIXED_NONREL_NEGPAR"),   
								  barNuclDMixedNRnegPar);

      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("XI_D_MIXED_NONREL"),   
								  barXiDMixedNR);

      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("NUCL_PATCH_MIXED_NONREL"),   
								  barNuclPatchMixedNR);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("DELTA_U_UNPOL"),
								  barDeltaUUnpol);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("DELTA_D_UNPOL"),
								  barDeltaDUnpol);
      
      return success;
    }

    bool registered = registerAll();
  } // namespace BaryonSeqSourceCallMapEnv


}  // end namespace Chroma
