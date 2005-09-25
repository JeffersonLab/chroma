// -*- C++ -*-
// $Id: mesonseqsrc_w.h,v 2.0 2005-09-25 21:04:35 edwards Exp $
/*! \file
 *  \brief Construct meson sequential sources.
 */

#ifndef __mesonseqsrc_w_h__
#define __mesonseqsrc_w_h__

namespace Chroma 
{

  LatticePropagator mesPionSeqSrc(const LatticePropagator& quark_propagator_1, 
				  const LatticePropagator& quark_propagator_2,
				  const LatticePropagator& quark_propagator_3);


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
				     const LatticePropagator& quark_propagator_3);


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
				       const LatticePropagator& quark_propagator_3);


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
				       const LatticePropagator& quark_propagator_3);


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
				      const LatticePropagator& quark_propagator_3);


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
				       const LatticePropagator& quark_propagator_3);


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
				       const LatticePropagator& quark_propagator_3);


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
				       const LatticePropagator& quark_propagator_3);


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
					const LatticePropagator& quark_propagator_3);


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
				     const LatticePropagator& quark_propagator_3);


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
				       const LatticePropagator& quark_propagator_3);


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
				       const LatticePropagator& quark_propagator_3);


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
				      const LatticePropagator& quark_propagator_3);


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
				       const LatticePropagator& quark_propagator_3);


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
				      const LatticePropagator& quark_propagator_3);


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
				      const LatticePropagator& quark_propagator_3);


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
					const LatticePropagator& quark_propagator_3);

}  // end namespace Chroma

#endif
