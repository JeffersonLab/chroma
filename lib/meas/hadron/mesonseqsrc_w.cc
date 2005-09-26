// $Id: mesonseqsrc_w.cc,v 2.1 2005-09-26 04:48:35 edwards Exp $
/*! \file
 *  \brief Construct meson sequential sources.
 */

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/mesonseqsrc_w.h"
#include "meas/hadron/seqsrc_funcmap_w.h"

namespace Chroma 
{

  //! Construct pion sequential source
  /*!
   * \ingroup hadron
   *
   *  delta(tz-tx) exp(i p.z) gamma_5 G gamma_5
   *
   * \param quark_propagators    array of quark propagators ( Read )
   *
   * \return the sequential source before projection onto the sink
   */

  LatticePropagator mesPionSeqSrc(const multi1d<LatticePropagator>& quark_propagators) 
  {
    START_CODE();

    LatticePropagator src_prop_tmp;

    if (quark_propagators.size() != 1)
    {
      QDPIO::cerr << __func__ << ": expect only 1 prop" << endl;
//      throw std::string("mesPionSeqSrc: expect only 1 prop");
      QDP_abort(1);
    }

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
//    src_prop_tmp = Gamma(G5) * quark_propagator_0 * Gamma(G5);

    src_prop_tmp = adj(quark_propagators[0]);

    END_CODE();

    return src_prop_tmp;
  }



  //! Construct pion_1-source something-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagators    array of quark propagators ( Read )
   * \param insertion            gamma insertion ( Read )
   *
   * \return gamma_5 * GAMMA^dag * gamma_5 * G * gamma_5
   */

  LatticePropagator mesPionXSeqSrc(const multi1d<LatticePropagator>& quark_propagators,
				   int insertion)
  {
    START_CODE();

    LatticePropagator fin;
    int G5 = Ns*Ns-1;

    if (quark_propagators.size() != 1)
    {
      QDPIO::cerr << __func__ << ": expect only 1 prop" << endl;
//      throw std::string("mesPionSeqSrc: expect only 1 prop");
      QDP_abort(1);
    }

    // src_prop_tmp = Gamma(G5) * Gamma(insertion)^dag * Gamma(G5) * quark_propagator_0 * Gamma(G5)
    // final = adj(Gamma(G5) * src_prop_tmp * Gamma(G5))
    //       = adj(Gamma(insertion)^dag * Gamma(G5) * quark_propagator_0)
    //       = adj(quark_propagator_0) * Gamma(G5) * Gamma(insertion)

    fin = adj(quark_propagators[0]) * Gamma(G5) * Gamma(insertion);

    END_CODE();

    return fin;
  }


  //! Construct pion-source a0_1-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagators    array of quark propagators ( Read )
   *
   * \return gamma_5 * 1 * gamma_5 * G * gamma_5
   */

  LatticePropagator mesPionA01SeqSrc(const multi1d<LatticePropagator>& quark_propagators) 
  {
    return mesPionXSeqSrc(quark_propagators,
			  0);   // a0_1 = 1
  }


  //! Construct pion-source rho_x_1-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagators    array of quark propagators ( Read )
   *
   * \return gamma_5 * gamma_1^dag * gamma_5 * G * gamma_5
   */

  LatticePropagator mesPionRhoX1SeqSrc(const multi1d<LatticePropagator>& quark_propagators) 
  {
    return mesPionXSeqSrc(quark_propagators,
			  1);   // rho_x_1 = gamma_1
  }


  //! Construct pion-source rho_y_1-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagators    array of quark propagators ( Read )
   *
   * \return gamma_5 * gamma_2^dag * gamma_5 * G * gamma_5
   */

  LatticePropagator mesPionRhoY1SeqSrc(const multi1d<LatticePropagator>& quark_propagators) 
  {
    return mesPionXSeqSrc(quark_propagators,
			  2);   // rho_y_1 = gamma_2
  }


  //! Construct pion-source b1_z_1-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagators    array of quark propagators ( Read )
   *
   * \return gamma_5 * (gamma_1*gamma_2)^dag * gamma_2  * gamma_5 * G * gamma_5
   */

  LatticePropagator mesPionB1Z1SeqSrc(const multi1d<LatticePropagator>& quark_propagators) 
  {
    return mesPionXSeqSrc(quark_propagators,
			  3);   // b1_z_1 = gamma_1*gamma_2
  }


  //! Construct pion-source rho_z_1-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagators    array of quark propagators ( Read )
   *
   * \return gamma_5 * gamma_3^dag * gamma_5 * G * gamma_5
   */

  LatticePropagator mesPionRhoZ1SeqSrc(const multi1d<LatticePropagator>& quark_propagators) 
  {
    return mesPionXSeqSrc(quark_propagators,
			  4);   // rho_z_1 = gamma_3
  }


  //! Construct pion_1-source b1_y_1-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagators    array of quark propagators ( Read )
   *
   * \return gamma_5 * (gamma_1*gamma_3)^dag * gamma_5 * G * gamma_5
   */

  LatticePropagator mesPion1B1Y1SeqSrc(const multi1d<LatticePropagator>& quark_propagators) 
  {
    return mesPionXSeqSrc(quark_propagators,
			  5);   // b1_y_1 =  gamma_1*gamma_3
  }


  //! Construct pion_1-source b1_x_1-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagators    array of quark propagators ( Read )
   *
   * \return gamma_5 * (gamma_2*gamma_3)^dag * gamma_5 * G * gamma_5
   */

  LatticePropagator mesPion1B1X1SeqSrc(const multi1d<LatticePropagator>& quark_propagators) 
  {
    return mesPionXSeqSrc(quark_propagators,
			  6);   // b1_x_1 =  gamma_2*gamma_3
  }


  //! Construct pion_1-source pion_2-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagators    array of quark propagators ( Read )
   *
   * \return gamma_5 * (gamma_4*gamma_5)^dag * gamma_5 * G * gamma_5
   */

  LatticePropagator mesPion1Pion2SeqSrc(const multi1d<LatticePropagator>& quark_propagators) 
  {
    return mesPionXSeqSrc(quark_propagators,
			  7);   // pion_2 =  gamma_4*gamma_5
  }


  //! Construct pion-source a0_2-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagators    array of quark propagators ( Read )
   *
   * \return gamma_5 * (gamma_4)^dag * gamma_5 * G * gamma_5
   */

  LatticePropagator mesPionA02SeqSrc(const multi1d<LatticePropagator>& quark_propagators) 
  {
    return mesPionXSeqSrc(quark_propagators,
			  8);   // a0_2 = gamma_4
  }


  //! Construct pion-source rho_x_2-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagators    array of quark propagators ( Read )
   *
   * \return gamma_5 * (gamma_4*gamma_1)^dag * gamma_5 * G * gamma_5
   */

  LatticePropagator mesPionRhoX2SeqSrc(const multi1d<LatticePropagator>& quark_propagators) 
  {
    return mesPionXSeqSrc(quark_propagators,
			  9);   // rho_x_2 = gamma_4 * gamma_1
  }


  //! Construct pion-source rho_y_2-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagators    array of quark propagators ( Read )
   *
   * \return gamma_5 * (gamma_4*gamma_2)^dag * gamma_5 * G * gamma_5
   */

  LatticePropagator mesPionRhoY2SeqSrc(const multi1d<LatticePropagator>& quark_propagators) 
  {
    return mesPionXSeqSrc(quark_propagators,
			  10);   // rho_y_2 = gamma_4 * gamma_2
  }


  //! Construct pion-source a1_z_1-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagators    array of quark propagators ( Read )
   *
   * \return gamma_5 * (gamma_1*gamma_2*gamma_4)^dag * gamma_5 * G * gamma_5
   */

  LatticePropagator mesPionA1Z1SeqSrc(const multi1d<LatticePropagator>& quark_propagators) 
  {
    return mesPionXSeqSrc(quark_propagators,
			  11);   // a1_z_1 = gamma_1 * gamma_2 * gamma_4
  }


  //! Construct pion-source rho_z_2-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagators    array of quark propagators ( Read )
   *
   * \return gamma_5 * (gamma_4*gamma_3)^dag * gamma_5 * G * gamma_5
   */

  LatticePropagator mesPionRhoZ2SeqSrc(const multi1d<LatticePropagator>& quark_propagators) 
  {
    return mesPionXSeqSrc(quark_propagators,
			  12);   // rho_z_2 = gamma_4 * gamma_3
  }


  //! Construct pion-source a1_y_1-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagators    array of quark propagators ( Read )
   *
   * \return gamma_5 * (gamma_1*gamma_3*gamma_4)^dag * gamma_5 * G * gamma_5
   */

  LatticePropagator mesPionA1Y1SeqSrc(const multi1d<LatticePropagator>& quark_propagators) 
  {
    return mesPionXSeqSrc(quark_propagators,
			  13);   // a1_y_1 = gamma_1 * gamma_3 * gamma_4
  }


  //! Construct pion-source a1_x_1-sink sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagators    array of quark propagators ( Read )
   *
   * \return gamma_5 * (gamma_2*gamma_3*gamma_4)^dag * gamma_5 * G * gamma_5
   */

  LatticePropagator mesPionA1X1SeqSrc(const multi1d<LatticePropagator>& quark_propagators) 
  {
    return mesPionXSeqSrc(quark_propagators,
			  14);   // a1_x_1 = gamma_2*gamma_3*gamma_4
  }


  //! Construct pion sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagators    array of quark propagators ( Read )
   *
   * \return gamma_5 * (gamma_5)^dag * gamma_5 * G * gamma_5
   */

  LatticePropagator mesPion1Pion1SeqSrc(const multi1d<LatticePropagator>& quark_propagators) 
  {
    return mesPionXSeqSrc(quark_propagators,
			  15);   // pion_1 = gamma_5
  }


  //! Meson sequential sources
  /*! \ingroup hadron */
  namespace MesonSeqSourceCallMapEnv
  { 
    bool registerAll(void) 
    {
      bool success = true;
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-A0_1"),
								  mesPionA01SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-RHO_X_1"),
								  mesPionRhoX1SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-RHO_Y_1"),
								  mesPionRhoY1SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-B1_Z_1"),
								  mesPionB1Z1SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-RHO_Z_1"),
								  mesPionRhoZ1SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-B1_Y_1"),
								  mesPionRhoZ1SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-B1_X_1"),
								  mesPionRhoZ1SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-PION_2"),
								  mesPion1Pion2SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-A0_2"),
								  mesPionA02SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-RHO_X_2"),
								  mesPionRhoX2SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-RHO_X_2"),
								  mesPionRhoX2SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-RHO_Y_2"),
								  mesPionRhoY2SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-A1_Z_1"),
								  mesPionA1Z1SeqSrc);
       
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-RHO_Z_2"),
								  mesPionRhoZ2SeqSrc);
       
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-A1_Y_1"),
								  mesPionA1Y1SeqSrc);
       
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-A1_X_1"),
								  mesPionA1X1SeqSrc);
       
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION"),
								  mesPion1Pion1SeqSrc);

      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION_1"), // same as PION
								  mesPion1Pion1SeqSrc);

      return success;
    }

    bool registered = registerAll();
  } // namespace MesonSeqSourceCallMapEnv


}  // end namespace Chroma
