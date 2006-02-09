// $Id: derivmesonseqsrc_w.cc,v 2.1 2006-02-09 02:25:24 edwards Exp $
/*! \file
 *  \brief Construct meson sequential sources.
 */

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/mesonseqsrc_w.h"
#include "meas/hadron/seqsrc_funcmap_w.h"

#include "meas/hadron/hadron_seqsource.h"

#include "util/ferm/symtensor.h"
#include "util/ferm/antisymtensor.h"

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



  //! Private namespace
  /*! * \ingroup hadron */
  namespace
  {
    //! Apply first deriv to the left onto source
    /*!
     * \ingroup hadron
     *
     * \f$\nabla_\mu f(x) = U_\mu(x)f(x+\mu) - U_{-\mu}(x)f(x-\mu)\f$
     *
     * \return $\f F(z,0)\nabla_\mu = F(x-\mu) U_\mu(x-\mu) - F(x+\mu) U_{-\mu}(x+\mu)\f$
     */
    LatticePropagator leftNabla(const LatticePropagator& F, 
				const multi1d<LatticeColorMatrix>& u,
				int mu)
    {
      return LatticePropagator(shift(F*u[mu], BACKWARD, mu) - shift(F, FORWARD, mu)*adj(u[mu]));
    }

    //! Apply "D_i" operator to the left onto source
    /*!
     * \ingroup hadron
     *
     * \f$D_i = s_{ijk}\nabla_j\nabla_k\f$
     *
     * where  \f$s_{ijk} = +1 \quad\forall i\ne j, j\ne k, i \ne k\f$
     * 
     * \return $\f F(z,0) D_\mu\f$
     */
    LatticePropagator leftD(const LatticePropagator& F,
			    const multi1d<LatticeColorMatrix>& u,
			    int mu)
    {
      LatticePropagator tmp = zero;

      // Slow implementation - to speed up could compute once the \nabla_j deriv
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (symTensor3d(mu,j,k) != 0)
	    tmp += leftNabla(leftNabla(F,u,j), u,k);
	}

      return tmp;
    }

    //! Apply "B_i" operator to the left onto source
    /*!
     * \ingroup hadron
     *
     * \f$B_i = \epsilon_{ijk}\nabla_j\nabla_k\f$
     *
     * \return $\f F(z,0) B_\mu\f$
     */
    LatticePropagator leftB(const LatticePropagator& F,
			    const multi1d<LatticeColorMatrix>& u,
			    int mu)
    {
      LatticePropagator tmp = zero;

      // Slow implementation - to speed up could compute once the \nabla_j deriv
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (antiSymTensor3d(mu,j,k) != 0)
	    tmp += Real(antiSymTensor3d(mu,j,k)) * leftNabla(leftNabla(F,u,j), u,k);
	}

      return tmp;
    }
  }


#if 0
  //! Construct pion_1-(A2=A1xD_T2) sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagators    array of quark propagators ( Read )
   *
   * Comes from a1xD_T2 in Manke's hep-lat/0210030 .
   * The sink interpolator is   \f$\Gamma_f \equiv \gamma_5\epsilon_{ijk}\gamma_j D_k\f$  
   * where   \f$D_i = s_{ijk}\nabla_j\nabla_k\f$
   * and     \f$\nabla_\mu f(x) = U_\mu(x)f(x+\mu) - U_{-\mu}(x)f(x-\mu)\f$
   *
   * \return \f$F^\dagger \gamma_5 \Gamma_f\f$
   */

  LatticePropagator mesPionA1xDT2SeqSrcMu(const multi1d<LatticePropagator>& quark_propagators,
					  const multi1d<LatticeColorMatrix>& u,
					  int mu)
  {
    START_CODE();

    LatticePropagator fin = zero;
    int G5 = Ns*Ns-1;

    if (quark_propagators.size() != 1)
    {
      QDPIO::cerr << __func__ << ": expect only 1 prop" << endl;
      QDP_abort(1);
    }

    LatticePropagator tmp = adj(quark_propagators[0]);

    // Slow implementation - to speed up could compute once the B_k deriv
    for(int j=0; j < 3; ++j)
      for(int k=0; k < 3; ++k)
      {
	if (antiSymTensor3d(mu,j,k) != 0)
	  tmp += Real(antiSymTensor3d(mu,j,k)) * (Gamma(1 << j) * leftD(tmp,u,k));
      }

    END_CODE();
    return fin;
  }

  //! Construct pion_1-(A2=A1xD_T2) sequential source
  /*!
   * \ingroup hadron
   *
   * \param quark_propagators    array of quark propagators ( Read )
   *
   * Comes from a1xD_T2 in Manke's hep-lat/0210030 .
   * The sink interpolator is   \f$\Gamma_f \equiv \gamma_5\epsilon_{ijk}\gamma_j D_k\f$  
   * where   \f$D_i = s_{ijk}\nabla_j\nabla_k\f$
   * and     \f$\nabla_\mu f(x) = U_\mu(x)f(x+\mu) - U_{-\mu}(x)f(x-\mu)\f$
   *
   * \return \f$F^\dagger \gamma_5 \Gamma_f\f$
   */

  LatticePropagator mesPionA1xDT2SeqSrcX(const multi1d<LatticePropagator>& quark_propagators,
					 const multi1d<LatticeColorMatrix>& u)
  {
    return mesPionA1xDT2SeqSrcMu(quark_propagators, 0);
  }

  LatticePropagator mesPionA1xDT2SeqSrcY(const multi1d<LatticePropagator>& quark_propagators,
					 const multi1d<LatticeColorMatrix>& u)
  {
    return mesPionA1xDT2SeqSrcMu(quark_propagators, 1);
  }

  LatticePropagator mesPionA1xDT2SeqSrcZ(const multi1d<LatticePropagator>& quark_propagators,
					 const multi1d<LatticeColorMatrix>& u)
  {
    return mesPionA1xDT2SeqSrcMu(quark_propagators, 2);
  }
#endif


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
								  mesPion1B1Y1SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-B1_X_1"),
								  mesPion1B1X1SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-PION_2"),
								  mesPion1Pion2SeqSrc);
      
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-A0_2"),
								  mesPionA02SeqSrc);
      
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

      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-PION_1"), // same as PION
								  mesPion1Pion1SeqSrc);

#if 0
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-A2=A1xD_T2_X"),
								  mesPionA1xDT2SeqSrcX);

      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-A2=A1xD_T2_Y"),
								  mesPionA1xDT2SeqSrcY);

      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-A2=A1xD_T2_Z"),
								  mesPionA1xDT2SeqSrcZ);
#endif

      return success;
    }

    bool registered = registerAll();
  } // namespace MesonSeqSourceCallMapEnv


}  // end namespace Chroma



namespace Chroma 
{

  // Read parameters
  void read(XMLReader& xml, const string& path, SimpleMesonSeqSourceEnv::Params& param)
  {
    SimpleMesonSeqSourceEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const SimpleMesonSeqSourceEnv::Params& param)
  {
    param.writeXML(xml, path);
  }



  //! Meson sequential sources
  /*! \ingroup hadron */
  namespace SimpleMesonSeqSourceEnv
  { 
    //! Initialize
    Params::Params() {}


    //! Read parameters
    Params::Params(XMLReader& xml, const string& path) {}


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);
      pop(xml);
    }

    //! Construct the source
    template<>
    LatticePropagator
    SimpleMesonSeqSource<LatticePropagator>::operator()(const multi1d<LatticeColorMatrix>& u,
							const multi1d<T>& forward_props,
							const multi1d<int>& sink_mom, 
							int t_sink, int j_decay) const
    {
      QDPIO::cout << "Simple meson sequential source: insertion = " << insertion << endl;

      return project(u,
		     mesPionXSeqSrc(quark_propagators, insertion),
		     sink_mom, t_sink, j_decay);
    }



    //! Callback function

    //! Construct pion-source a0_1-sink sequential source
    /*!
     * \ingroup hadron
     *
     * \return gamma_5 * 1 * gamma_5 * G * gamma_5
     */
    HadronSeqSource<LatticePropagator>* mesPionA01SeqSrc(XMLReader& xml_in,
							 const std::string& path)
    {
      return new SimpleMesonSeqSource<LatticePropagator>(0);   // a0_1 = 1
    }


    //! Construct pion-source rho_x_1-sink sequential source
    /*!
     * \ingroup hadron
     *
     * \return gamma_5 * gamma_1^dag * gamma_5 * G * gamma_5
     */
    HadronSeqSource<LatticePropagator>*  mesPionRhoX1SeqSrc(XMLReader& xml_in,
							    const std::string& path)
    {
      return new SimpleMesonSeqSource<LatticePropagator>(1);   // rho_x_1 = gamma_1
    }


    //! Construct pion-source rho_y_1-sink sequential source
    /*!
     * \ingroup hadron
     *
     * \return gamma_5 * gamma_2^dag * gamma_5 * G * gamma_5
     */
    HadronSeqSource<LatticePropagator>*  mesPionRhoY1SeqSrc(XMLReader& xml_in,
							    const std::string& path)
    {
      return new SimpleMesonSeqSource<LatticePropagator>(2);   // rho_y_1 = gamma_2
    }


    //! Construct pion-source b1_z_1-sink sequential source
    /*!
     * \ingroup hadron
     *
     * \return gamma_5 * (gamma_1*gamma_2)^dag * gamma_2  * gamma_5 * G * gamma_5
     */
    HadronSeqSource<LatticePropagator>* mesPionB1Z1SeqSrc(XMLReader& xml_in,
							  const std::string& path)
    {
      return new SimpleMesonSeqSource<LatticePropagator>(3);   // b1_z_1 = gamma_1*gamma_2
    }


    //! Construct pion-source rho_z_1-sink sequential source
    /*!
     * \ingroup hadron
     *
     * \return gamma_5 * gamma_3^dag * gamma_5 * G * gamma_5
     */
    HadronSeqSource<LatticePropagator>* mesPionRhoZ1SeqSrc(XMLReader& xml_in,
							   const std::string& path)
    {
      return new SimpleMesonSeqSource<LatticePropagator>(4);   // rho_z_1 = gamma_3
    }


    //! Construct pion_1-source b1_y_1-sink sequential source
    /*!
     * \ingroup hadron
     *
     * \return gamma_5 * (gamma_1*gamma_3)^dag * gamma_5 * G * gamma_5
     */
    HadronSeqSource<LatticePropagator>* mesPion1B1Y1SeqSrc(XMLReader& xml_in,
							   const std::string& path)
    {
      return new SimpleMesonSeqSource<LatticePropagator>(5);   // b1_y_1 =  gamma_1*gamma_3
    }


    //! Construct pion_1-source b1_x_1-sink sequential source
    /*!
     * \ingroup hadron
     *
     * \return gamma_5 * (gamma_2*gamma_3)^dag * gamma_5 * G * gamma_5
     */
    HadronSeqSource<LatticePropagator>* mesPion1B1X1SeqSrc(XMLReader& xml_in,
							   const std::string& path)
    {
      return new SimpleMesonSeqSource<LatticePropagator>(6);   // b1_x_1 =  gamma_2*gamma_3
    }


    //! Construct pion_1-source pion_2-sink sequential source
    /*!
     * \ingroup hadron
     *
     * \return gamma_5 * (gamma_4*gamma_5)^dag * gamma_5 * G * gamma_5
     */
    HadronSeqSource<LatticePropagator>* mesPion1Pion2SeqSrc(XMLReader& xml_in,
							    const std::string& path)
    {
      return new SimpleMesonSeqSource<LatticePropagator>(7);   // pion_2 =  gamma_4*gamma_5
    }


    //! Construct pion-source a0_2-sink sequential source
    /*!
     * \ingroup hadron
     *
     * \return gamma_5 * (gamma_4)^dag * gamma_5 * G * gamma_5
     */
    HadronSeqSource<LatticePropagator>* mesPionA02SeqSrc(XMLReader& xml_in,
							 const std::string& path)
    {
      return new SimpleMesonSeqSource<LatticePropagator>(8);   // a0_2 = gamma_4
    }


    //! Construct pion-source rho_x_2-sink sequential source
    /*!
     * \ingroup hadron
     *
     * \return gamma_5 * (gamma_4*gamma_1)^dag * gamma_5 * G * gamma_5
     */
    HadronSeqSource<LatticePropagator>* mesPionRhoX2SeqSrc(XMLReader& xml_in,
							   const std::string& path)
    {
      return new SimpleMesonSeqSource<LatticePropagator>(9);   // rho_x_2 = gamma_4 * gamma_1
    }


    //! Construct pion-source rho_y_2-sink sequential source
    /*!
     * \ingroup hadron
     *
     * \return gamma_5 * (gamma_4*gamma_2)^dag * gamma_5 * G * gamma_5
     */
    HadronSeqSource<LatticePropagator>* mesPionRhoY2SeqSrc(XMLReader& xml_in,
							 const std::string& path)
    {
      return new SimpleMesonSeqSource<LatticePropagator>(10);   // rho_y_2 = gamma_4 * gamma_2
    }


    //! Construct pion-source a1_z_1-sink sequential source
    /*!
     * \ingroup hadron
     *
     * \return gamma_5 * (gamma_1*gamma_2*gamma_4)^dag * gamma_5 * G * gamma_5
     */
    HadronSeqSource<LatticePropagator>* mesPionA1Z1SeqSrc(XMLReader& xml_in,
							  const std::string& path)
    {
      return new SimpleMesonSeqSource<LatticePropagator>(11);   // a1_z_1 = gamma_1 * gamma_2 * gamma_4
    }


    //! Construct pion-source rho_z_2-sink sequential source
    /*!
     * \ingroup hadron
     *
     * \return gamma_5 * (gamma_4*gamma_3)^dag * gamma_5 * G * gamma_5
     */
    HadronSeqSource<LatticePropagator>* mesPionRhoZ2SeqSrc(XMLReader& xml_in,
							   const std::string& path)
    {
      return new SimpleMesonSeqSource<LatticePropagator>(12);   // rho_z_2 = gamma_4 * gamma_3
    }


    //! Construct pion-source a1_y_1-sink sequential source
    /*!
     * \ingroup hadron
     *
     * \return gamma_5 * (gamma_1*gamma_3*gamma_4)^dag * gamma_5 * G * gamma_5
     */
    HadronSeqSource<LatticePropagator>* mesPionA1Y1SeqSrc(XMLReader& xml_in,
							  const std::string& path)
    {
      return new SimpleMesonSeqSource<LatticePropagator>(13);   // a1_y_1 = gamma_1 * gamma_3 * gamma_4
    }


    //! Construct pion-source a1_x_1-sink sequential source
    /*!
     * \ingroup hadron
     *
     * \return gamma_5 * (gamma_2*gamma_3*gamma_4)^dag * gamma_5 * G * gamma_5
     */
    HadronSeqSource<LatticePropagator>* mesPionA1X1SeqSrc(XMLReader& xml_in,
							  const std::string& path)
    {
      return new SimpleMesonSeqSource<LatticePropagator>(14);   // a1_x_1 = gamma_2*gamma_3*gamma_4
    }


    //! Construct pion sequential source
    /*!
     * \ingroup hadron
     *
     * \return gamma_5 * (gamma_5)^dag * gamma_5 * G * gamma_5
     */
    HadronSeqSource<LatticePropagator>* mesPion1Pion1SeqSrc(XMLReader& xml_in,
							    const std::string& path)
    {
      return new SimpleMesonSeqSource<LatticePropagator>(15);   // pion_1 = gamma_5
    }

    // Register all the possible simple mesons
    bool registerAll(void) 
    {
      bool success = true;

      //! Register all the factories
//      success &= LinkSmearingEnv::registered;
      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-A0_1"), 
										    mesPionA01SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-RHO_X_1"), 
										    mesPionRhoX1SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-RHO_Y_1"),
										    mesPionRhoY1SeqSrc);
      
      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-B1_Z_1"),
										    mesPionB1Z1SeqSrc);
      
      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-RHO_Z_1"),
										    mesPionRhoZ1SeqSrc);
      
      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-B1_Y_1"),
										    mesPion1B1Y1SeqSrc);
      
      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-B1_X_1"),
										    mesPion1B1X1SeqSrc);
      
      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-PION_2"),
										    mesPion1Pion2SeqSrc);
      
      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-A0_2"),
										    mesPionA02SeqSrc);
      
      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-RHO_X_2"),
										    mesPionRhoX2SeqSrc);
      
      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-RHO_Y_2"),
										    mesPionRhoY2SeqSrc);
      
      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-A1_Z_1"),
										    mesPionA1Z1SeqSrc);
       
      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-RHO_Z_2"),
										    mesPionRhoZ2SeqSrc);
       
      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-A1_Y_1"),
										    mesPionA1Y1SeqSrc);
       
      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-A1_X_1"),
										    mesPionA1X1SeqSrc);
       
      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION"),
										    mesPion1Pion1SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION_1"), // same as PION
										    mesPion1Pion1SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-PION_1"), // same as PION
										    mesPion1Pion1SeqSrc);

#if 0
      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-A2=A1xD_T2_X"),
								  mesPionA1xDT2SeqSrcX);

      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-A2=A1xD_T2_Y"),
								  mesPionA1xDT2SeqSrcY);

      success &= TheSeqSourceFuncMap::Instance().registerFunction(string("PION-A2=A1xD_T2_Z"),
								  mesPionA1xDT2SeqSrcZ);
#endif

      return success;
    }

    bool registered = registerAll();

  }  // end namespace SimpleMesonSeqSourceEnv

}  // end namespace Chroma


  
