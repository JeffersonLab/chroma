// $Id: mesonseqsrc_w.cc,v 3.0 2006-04-03 04:59:00 edwards Exp $
/*! \file
 *  \brief Construct meson sequential sources.
 */

#include "meas/hadron/mesonseqsrc_w.h"
#include "meas/hadron/seqsource_factory_w.h"

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
    //! Anonymous namespace
    namespace
    {
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


      //-------------------- callback functions ---------------------------------------

      //! Construct pion-source a0_1-sink sequential source
      /*!
       * \ingroup hadron
       *
       * \return gamma_5 * 1 * gamma_5 * G * gamma_5
       */
      HadronSeqSource<LatticePropagator>* mesPionA01SeqSrc(XMLReader& xml_in,
							   const std::string& path)
      {
	return new SimpleMesonSeqSource(Params(xml_in, path), 0);   // a0_1 = 1
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
	return new SimpleMesonSeqSource(Params(xml_in, path), 1);   // rho_x_1 = gam ma_1
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
	return new SimpleMesonSeqSource(Params(xml_in, path), 2);   // rho_y_1 = gamma_2
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
	return new SimpleMesonSeqSource(Params(xml_in, path), 3);   // b1_z_1 = gamma_1*gamma_2
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
	return new SimpleMesonSeqSource(Params(xml_in, path), 4);   // rho_z_1 = gamma_3
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
	return new SimpleMesonSeqSource(Params(xml_in, path), 5);   // b1_y_1 =  gamma_1*gamma_3
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
	return new SimpleMesonSeqSource(Params(xml_in, path), 6);   // b1_x_1 =  gamma_2*gamma_3
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
	return new SimpleMesonSeqSource(Params(xml_in, path), 7);   // pion_2 =  gamma_4*gamma_5
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
	return new SimpleMesonSeqSource(Params(xml_in, path), 8);   // a0_2 = gamma_4
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
	return new SimpleMesonSeqSource(Params(xml_in, path), 9);   // rho_x_2 = gamma_4 * gamma_1
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
	return new SimpleMesonSeqSource(Params(xml_in, path), 10);   // rho_y_2 = gamma_4 * gamma_2
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
	return new SimpleMesonSeqSource(Params(xml_in, path), 11);   // a1_z_1 = gamma_1 * gamma_2 * gamma_4
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
	return new SimpleMesonSeqSource(Params(xml_in, path), 12);   // rho_z_2 = gamma_4 * gamma_3
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
	return new SimpleMesonSeqSource(Params(xml_in, path), 13);   // a1_y_1 = gamma_1 * gamma_3 * gamma_4
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
	return new SimpleMesonSeqSource(Params(xml_in, path), 14);   // a1_x_1 = gamma_2*gamma_3*gamma_4
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
	return new SimpleMesonSeqSource(Params(xml_in, path), 15);   // pion_1 = gamma_5
      }

    } // end anonymous namespace


    //! Initialize
    Params::Params()
    {
      j_decay = -1;
      t_sink  = -1;
      sink_mom.resize(Nd-1);
      sink_mom = 0;
    }


    //! Read parameters
    Params::Params(XMLReader& xml, const string& path)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
	break;

      default:
	QDPIO::cerr << __func__ << ": parameter version " << version 
		    << " unsupported." << endl;
	QDP_abort(1);
      }

      read(paramtop, "j_decay", j_decay);
      read(paramtop, "t_sink", t_sink);
      read(paramtop, "sink_mom", sink_mom);
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);

      write(xml, "j_decay", j_decay);
      write(xml, "t_sink", t_sink);
      write(xml, "sink_mom", sink_mom);
      pop(xml);
    }


    //! Construct the source
    LatticePropagator
    SimpleMesonSeqSource::operator()(const multi1d<LatticeColorMatrix>& u,
				     const multi1d<LatticePropagator>& quark_propagators) const
    {
      QDPIO::cout << "Simple meson sequential source: gamma_sink = " << gamma_sink << endl;

      return project(mesPionXSeqSrc(quark_propagators, gamma_sink),
		     params.sink_mom, params.t_sink, params.j_decay);
    }


    // Register all the possible simple mesons
    bool registerAll(void) 
    {
      bool success = true;

      //! Register all the factories
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

      return success;
    }

    const bool registered = registerAll();

  }  // end namespace SimpleMesonSeqSourceEnv

}  // end namespace Chroma


  
