// $Id: mesonseqsrc_w.cc,v 1.8 2005-03-18 05:12:37 edwards Exp $
/*! \file
 *  \brief Construct meson sequential sources.
 */

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/mesonseqsrc_w.h"

namespace Chroma 
{

  //! Construct pion sequential source
  /*!
   * \ingroup hadron
   *
   *  delta(tz-tx) exp(i p.z) \gamma_5 G \gamma_5
   *
   * \param quark_propagator_1   first (u) quark propagator ( Read )
   * \param quark_propagator_2   second (d) quark propagator ( Read )
   * \param quark_propagator_3   third (s) quark propagator ( Read )
   *
   * \return the sequential source before projection onto the sink
   */

  LatticePropagator mesPionSeqSrc(const LatticePropagator& quark_propagator_1, 
				  const LatticePropagator& quark_propagator_2,
				  const LatticePropagator& quark_propagator_3)
  {
    START_CODE();

    LatticePropagator src_prop_tmp;

    /*
     * Note, the seq. source is simply   src_prop_tmp = gamma_5 * U * gamma_5
     * where in this encoding gamma_5 = Gamma(15) .
     * However, to maintain compatibility with the calling  hadseqsrc_w.cc
     * code, the seqsource defined here is just   adj(U) .
     * That's because the hadseqsrc routine (for compatibility with baryons)
     * will construct finally   gamma_5 * adj(src_prop_tmp) * gamma_5 .
     * So, only the adj() is applied here to compensate for the adj()
     * that is subsequently done. The gamma_5 come along for free.
     */
    
//    int G5 = Ns*Ns-1;
//    src_prop_tmp = Gamma(G5) * quark_propagator_1 * Gamma(G5);

    src_prop_tmp = adj(quark_propagator_1);

    END_CODE();

    return src_prop_tmp;
  }



  //! Construct pion_1-source something-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagator_1   first (u) quark propagator ( Read )
   * \param quark_propagator_2   second (d) quark propagator ( Read )
   * \param quark_propagator_3   third (s) quark propagator ( Read )
   * \param insertion            gamma insertion ( Read )
   *
   * \return \gamma_5 * GAMMA^dag * \gamma_5 * G * \gamma_5
   */

  LatticePropagator mesPionXSeqSrc(const LatticePropagator& quark_propagator_1, 
				   const LatticePropagator& quark_propagator_2,
				   const LatticePropagator& quark_propagator_3,
				   int insertion)
  {
    START_CODE();

    LatticePropagator fin;
    int G5 = Ns*Ns-1;

    // src_prop_tmp = Gamma(G5) * Gamma(insertion)^dag * Gamma(G5) * quark_propagator_1 * Gamma(G5)
    // final = adj(Gamma(G5) * src_prop_tmp * Gamma(G5))
    //       = adj(Gamma(insertion)^dag * Gamma(G5) * quark_propagator_1)
    //       = adj(quark_propagator_1) * Gamma(G5) * Gamma(insertion)

    fin = adj(quark_propagator_1) * Gamma(G5) * Gamma(insertion);

    END_CODE();

    return fin;
  }


  //! Construct pion-source a0_1-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagator_1   first (u) quark propagator ( Read )
   * \param quark_propagator_2   second (d) quark propagator ( Read )
   * \param quark_propagator_3   third (s) quark propagator ( Read )
   *
   * \return \gamma_5 * 1 * gamma_5 * G * \gamma_5
   */

  LatticePropagator mesPionA01SeqSrc(const LatticePropagator& quark_propagator_1, 
				     const LatticePropagator& quark_propagator_2,
				     const LatticePropagator& quark_propagator_3)
  {
    return mesPionXSeqSrc(quark_propagator_1,
			  quark_propagator_2,
			  quark_propagator_3,
			  0);   // a0_1 = 1
  }


  //! Construct pion-source rho_x_1-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagator_1   first (u) quark propagator ( Read )
   * \param quark_propagator_2   second (d) quark propagator ( Read )
   * \param quark_propagator_3   third (s) quark propagator ( Read )
   *
   * \return \gamma_5 * gamma_1^dag * gamma_5 * G * \gamma_5
   */

  LatticePropagator mesPionRhoX1SeqSrc(const LatticePropagator& quark_propagator_1, 
				       const LatticePropagator& quark_propagator_2,
				       const LatticePropagator& quark_propagator_3)
  {
    return mesPionXSeqSrc(quark_propagator_1,
			  quark_propagator_2,
			  quark_propagator_3,
			  1);   // rho_x_1 = gamma_1
  }


  //! Construct pion-source rho_y_1-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagator_1   first (u) quark propagator ( Read )
   * \param quark_propagator_2   second (d) quark propagator ( Read )
   * \param quark_propagator_3   third (s) quark propagator ( Read )
   *
   * \return \gamma_5 * gamma_2^dag * gamma_5 * G * \gamma_5
   */

  LatticePropagator mesPionRhoY1SeqSrc(const LatticePropagator& quark_propagator_1, 
				       const LatticePropagator& quark_propagator_2,
				       const LatticePropagator& quark_propagator_3)
  {
    return mesPionXSeqSrc(quark_propagator_1,
			  quark_propagator_2,
			  quark_propagator_3,
			  2);   // rho_y_1 = gamma_2
  }


  //! Construct pion-source b1_z_1-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagator_1   first (u) quark propagator ( Read )
   * \param quark_propagator_2   second (d) quark propagator ( Read )
   * \param quark_propagator_3   third (s) quark propagator ( Read )
   *
   * \return \gamma_5 * (gamma_1*gamma_2)^dag * gamma_2  * gamma_5 * G * \gamma_5
   */

  LatticePropagator mesPionB1Z1SeqSrc(const LatticePropagator& quark_propagator_1, 
				      const LatticePropagator& quark_propagator_2,
				      const LatticePropagator& quark_propagator_3)
  {
    return mesPionXSeqSrc(quark_propagator_1,
			  quark_propagator_2,
			  quark_propagator_3,
			  3);   // b1_z_1 = gamma_1*gamma_2
  }


  //! Construct pion-source rho_z_1-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagator_1   first (u) quark propagator ( Read )
   * \param quark_propagator_2   second (d) quark propagator ( Read )
   * \param quark_propagator_3   third (s) quark propagator ( Read )
   *
   * \return \gamma_5 * gamma_3^dag * gamma_5 * G * \gamma_5
   */

  LatticePropagator mesPionRhoZ1SeqSrc(const LatticePropagator& quark_propagator_1, 
				       const LatticePropagator& quark_propagator_2,
				       const LatticePropagator& quark_propagator_3)
  {
    return mesPionXSeqSrc(quark_propagator_1,
			  quark_propagator_2,
			  quark_propagator_3,
			  4);   // rho_z_1 = gamma_3
  }


  //! Construct pion_1-source b1_y_1-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagator_1   first (u) quark propagator ( Read )
   * \param quark_propagator_2   second (d) quark propagator ( Read )
   * \param quark_propagator_3   third (s) quark propagator ( Read )
   *
   * \return \gamma_5 * (gamma_1*gamma_3)^dag * gamma_5 * G * \gamma_5
   */

  LatticePropagator mesPion1B1Y1SeqSrc(const LatticePropagator& quark_propagator_1, 
				       const LatticePropagator& quark_propagator_2,
				       const LatticePropagator& quark_propagator_3)
  {
    return mesPionXSeqSrc(quark_propagator_1,
			  quark_propagator_2,
			  quark_propagator_3,
			  5);   // b1_y_1 =  gamma_1*gamma_3
  }


  //! Construct pion_1-source b1_x_1-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagator_1   first (u) quark propagator ( Read )
   * \param quark_propagator_2   second (d) quark propagator ( Read )
   * \param quark_propagator_3   third (s) quark propagator ( Read )
   *
   * \return \gamma_5 * (gamma_2*gamma_3)^dag * gamma_5 * G * \gamma_5
   */

  LatticePropagator mesPion1B1X1SeqSrc(const LatticePropagator& quark_propagator_1, 
				       const LatticePropagator& quark_propagator_2,
				       const LatticePropagator& quark_propagator_3)
  {
    return mesPionXSeqSrc(quark_propagator_1,
			  quark_propagator_2,
			  quark_propagator_3,
			  6);   // b1_x_1 =  gamma_2*gamma_3
  }


  //! Construct pion_1-source pion_2-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagator_1   first (u) quark propagator ( Read )
   * \param quark_propagator_2   second (d) quark propagator ( Read )
   * \param quark_propagator_3   third (s) quark propagator ( Read )
   *
   * \return \gamma_5 * (gamma_4*gamma_5)^dag * gamma_5 * G * \gamma_5
   */

  LatticePropagator mesPion1Pion2SeqSrc(const LatticePropagator& quark_propagator_1, 
					const LatticePropagator& quark_propagator_2,
					const LatticePropagator& quark_propagator_3)
  {
    return mesPionXSeqSrc(quark_propagator_1,
			  quark_propagator_2,
			  quark_propagator_3,
			  7);   // pion_2 =  gamma_4*gamma_5
  }


  //! Construct pion-source a0_2-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagator_1   first (u) quark propagator ( Read )
   * \param quark_propagator_2   second (d) quark propagator ( Read )
   * \param quark_propagator_3   third (s) quark propagator ( Read )
   *
   * \return \gamma_5 * (gamma_4)^dag * gamma_5 * G * \gamma_5
   */

  LatticePropagator mesPionA02SeqSrc(const LatticePropagator& quark_propagator_1, 
				     const LatticePropagator& quark_propagator_2,
				     const LatticePropagator& quark_propagator_3)
  {
    return mesPionXSeqSrc(quark_propagator_1,
			  quark_propagator_2,
			  quark_propagator_3,
			  8);   // a0_2 = gamma_4
  }


  //! Construct pion-source rho_x_2-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagator_1   first (u) quark propagator ( Read )
   * \param quark_propagator_2   second (d) quark propagator ( Read )
   * \param quark_propagator_3   third (s) quark propagator ( Read )
   *
   * \return \gamma_5 * (gamma_4*gamma_1)^dag * gamma_5 * G * \gamma_5
   */

  LatticePropagator mesPionRhoX2SeqSrc(const LatticePropagator& quark_propagator_1, 
				       const LatticePropagator& quark_propagator_2,
				       const LatticePropagator& quark_propagator_3)
  {
    return mesPionXSeqSrc(quark_propagator_1,
			  quark_propagator_2,
			  quark_propagator_3,
			  9);   // rho_x_2 = gamma_4 * gamma_1
  }


  //! Construct pion-source rho_y_2-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagator_1   first (u) quark propagator ( Read )
   * \param quark_propagator_2   second (d) quark propagator ( Read )
   * \param quark_propagator_3   third (s) quark propagator ( Read )
   *
   * \return \gamma_5 * (gamma_4*gamma_2)^dag * gamma_5 * G * \gamma_5
   */

  LatticePropagator mesPionRhoY2SeqSrc(const LatticePropagator& quark_propagator_1, 
				       const LatticePropagator& quark_propagator_2,
				       const LatticePropagator& quark_propagator_3)
  {
    return mesPionXSeqSrc(quark_propagator_1,
			  quark_propagator_2,
			  quark_propagator_3,
			  10);   // rho_y_2 = gamma_4 * gamma_2
  }


  //! Construct pion-source a1_z_1-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagator_1   first (u) quark propagator ( Read )
   * \param quark_propagator_2   second (d) quark propagator ( Read )
   * \param quark_propagator_3   third (s) quark propagator ( Read )
   *
   * \return \gamma_5 * (gamma_1*gamma_2*gamma_4)^dag * gamma_5 * G * \gamma_5
   */

  LatticePropagator mesPionA1Z1SeqSrc(const LatticePropagator& quark_propagator_1, 
				      const LatticePropagator& quark_propagator_2,
				      const LatticePropagator& quark_propagator_3)
  {
    return mesPionXSeqSrc(quark_propagator_1,
			  quark_propagator_2,
			  quark_propagator_3,
			  11);   // a1_z_1 = gamma_1 * gamma_2 * gamma_4
  }


  //! Construct pion-source rho_z_2-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagator_1   first (u) quark propagator ( Read )
   * \param quark_propagator_2   second (d) quark propagator ( Read )
   * \param quark_propagator_3   third (s) quark propagator ( Read )
   *
   * \return \gamma_5 * (gamma_4*gamma_3)^dag * gamma_5 * G * \gamma_5
   */

  LatticePropagator mesPionRhoZ2SeqSrc(const LatticePropagator& quark_propagator_1, 
				       const LatticePropagator& quark_propagator_2,
				       const LatticePropagator& quark_propagator_3)
  {
    return mesPionXSeqSrc(quark_propagator_1,
			  quark_propagator_2,
			  quark_propagator_3,
			  12);   // rho_z_2 = gamma_4 * gamma_3
  }


  //! Construct pion-source a1_y_1-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagator_1   first (u) quark propagator ( Read )
   * \param quark_propagator_2   second (d) quark propagator ( Read )
   * \param quark_propagator_3   third (s) quark propagator ( Read )
   *
   * \return \gamma_5 * (gamma_1*gamma_3*gamma_4)^dag * gamma_5 * G * \gamma_5
   */

  LatticePropagator mesPionA1Y1SeqSrc(const LatticePropagator& quark_propagator_1, 
				       const LatticePropagator& quark_propagator_2,
				       const LatticePropagator& quark_propagator_3)
  {
    return mesPionXSeqSrc(quark_propagator_1,
			  quark_propagator_2,
			  quark_propagator_3,
			  13);   // a1_y_1 = gamma_1 * gamma_3 * gamma_4
  }


  //! Construct pion-source a1_x_1-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagator_1   first (u) quark propagator ( Read )
   * \param quark_propagator_2   second (d) quark propagator ( Read )
   * \param quark_propagator_3   third (s) quark propagator ( Read )
   *
   * \return \gamma_5 * (gamma_2*gamma_3*gamma_4)^dag * gamma_5 * G * \gamma_5
   */

  LatticePropagator mesPionA1X1SeqSrc(const LatticePropagator& quark_propagator_1, 
				      const LatticePropagator& quark_propagator_2,
				      const LatticePropagator& quark_propagator_3)
  {
    return mesPionXSeqSrc(quark_propagator_1,
			  quark_propagator_2,
			  quark_propagator_3,
			  14);   // a1_x_1 = gamma_2*gamma_3*gamma_4
  }


  //! Construct pion sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagator_1   first (u) quark propagator ( Read )
   * \param quark_propagator_2   second (d) quark propagator ( Read )
   * \param quark_propagator_3   third (s) quark propagator ( Read )
   *
   * \return \gamma_5 * (gamma_5)^dag * gamma_5 * G * \gamma_5
   */

  LatticePropagator mesPion1Pion1SeqSrc(const LatticePropagator& quark_propagator_1, 
					const LatticePropagator& quark_propagator_2,
					const LatticePropagator& quark_propagator_3)
  {
    return mesPionXSeqSrc(quark_propagator_1,
			  quark_propagator_2,
			  quark_propagator_3,
			  15);   // pion_1 = gamma_5
  }

}  // end namespace Chroma
