// $Id: simple_meson_seqsrc_w.cc,v 3.8 2006-12-04 20:38:29 edwards Exp $
/*! \file
 *  \brief Construct meson sequential sources.
 */

#include "meas/hadron/simple_meson_seqsrc_w.h"
#include "meas/hadron/seqsource_factory_w.h"
#include "util/ferm/gamma5_herm_w.h"
#include "util/ft/sftmom.h"

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
      //! Compute final gamma insertion
      /*!
       * \ingroup hadron
       *
       * \return \f I = $\gamma_5 * \Gamma_f^dag * \gamma_5\$
       */
      int sign_g5adjGfg5(int sink_insertion)
      {
	int sgn = 1;
	// The algorithm simply counts bits. If sink_insertion has 1 or 2 bits, there is a 
	// sign change, 0, 3 or 4 bits no sign change.
	int cnt = 0;
	if ((sink_insertion & 1) != 0) cnt++;
	if ((sink_insertion & 2) != 0) cnt++;
	if ((sink_insertion & 4) != 0) cnt++;
	if ((sink_insertion & 8) != 0) cnt++;

	switch (cnt)
	{
	case 0:
	case 3:
	case 4:
	  sgn = 1;
	  break;

	case 1:
	case 2:
	  sgn = -1;
	  break;

	default:
	  QDPIO::cerr << __func__ << ": internal error - number of bits not supported" << endl;
	  QDP_abort(1);
	}

	return sgn;
      }


      //! Construct the source
      /*!
       * \ingroup hadron
       *
       * \f$\gamma_5 * \Gamma(gamma_sink)^\dag * \gamma_5 * F\$
       */
      LatticePropagator mesA0XSeqSrc(const LatticePropagator& quark_prop, int gamma_sink)
      {
	START_CODE();
	
	QDPIO::cout << "Simple meson sequential source: gamma_sink = " << gamma_sink << endl;

	// \f$\gamma_5 * \Gamma_f^dag * \gamma_5 * F\$
//      LatticePropagator fin = Gamma(G5) * adj(Gamma(gamma_sink)) * Gamma(G5) * quark_prop;

	LatticePropagator fin = sign_g5adjGfg5(gamma_sink) * (Gamma(gamma_sink) * quark_prop);

	END_CODE();
	
	return fin;
      }


      //-------------------- callback functions ---------------------------------------

      //! Construct a0-source a0-sink sequential source
      /*!
       * \ingroup hadron
       *
       * \return \f$\gamma_5 * \Gamma(0)^dag * \gamma_5 * F\$
       */
      HadronSeqSource<LatticePropagator>* mesA0A01SeqSrc(XMLReader& xml_in,
							 const std::string& path)
      {
	return new SimpleMesonSeqSource(Params(xml_in, path), 0);   // a0 = 1
      }


      //! Construct a0-source rho_x_1-sink sequential source
      /*!
       * \ingroup hadron
       *
       * \return \f$\gamma_5 * \Gamma(1)^dag * \gamma_5 * F\$
       */
      HadronSeqSource<LatticePropagator>*  mesA0RhoX1SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new SimpleMesonSeqSource(Params(xml_in, path), 1);   // rho_x_1 = gam ma_1
      }


      //! Construct a0-source rho_y_1-sink sequential source
      /*!
       * \ingroup hadron
       *
       * \return \f$\gamma_5 * \Gamma(2)^dag * \gamma_5 * F\$
       */
      HadronSeqSource<LatticePropagator>*  mesA0RhoY1SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new SimpleMesonSeqSource(Params(xml_in, path), 2);   // rho_y_1 = gamma_2
      }


      //! Construct a0-source b1_z_1-sink sequential source
      /*!
       * \ingroup hadron
       *
       * \return \f$\gamma_5 * \Gamma(3)^dag * \gamma_5 * F\$
       */
      HadronSeqSource<LatticePropagator>* mesA0B1Z1SeqSrc(XMLReader& xml_in,
							  const std::string& path)
      {
	return new SimpleMesonSeqSource(Params(xml_in, path), 3);   // b1_z_1 = gamma_1*gamma_2
      }


      //! Construct a0-source rho_z_1-sink sequential source
      /*!
       * \ingroup hadron
       *
       * \return \f$\gamma_5 * \Gamma(4)^dag * \gamma_5 * F\$
       */
      HadronSeqSource<LatticePropagator>* mesA0RhoZ1SeqSrc(XMLReader& xml_in,
							   const std::string& path)
      {
	return new SimpleMesonSeqSource(Params(xml_in, path), 4);   // rho_z_1 = gamma_3
      }


      //! Construct a0-source b1_y_1-sink sequential source
      /*!
       * \ingroup hadron
       *
       * \return \f$\gamma_5 * \Gamma(5)^dag * \gamma_5 * F\$
       */
      HadronSeqSource<LatticePropagator>* mesA01B1Y1SeqSrc(XMLReader& xml_in,
							   const std::string& path)
      {
	return new SimpleMesonSeqSource(Params(xml_in, path), 5);   // b1_y_1 =  gamma_1*gamma_3
      }


      //! Construct a0-source b1_x_1-sink sequential source
      /*!
       * \ingroup hadron
       *
       * \return \f$\gamma_5 * \Gamma(6)^dag * \gamma_5 * F\$
       */
      HadronSeqSource<LatticePropagator>* mesA01B1X1SeqSrc(XMLReader& xml_in,
							   const std::string& path)
      {
	return new SimpleMesonSeqSource(Params(xml_in, path), 6);   // b1_x_1 =  gamma_2*gamma_3
      }


      //! Construct a0-source pion_2-sink sequential source
      /*!
       * \ingroup hadron
       *
       * \return \f$\gamma_5 * \Gamma(7)^dag * \gamma_5 * F\$
       */
      HadronSeqSource<LatticePropagator>* mesA01Pion2SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new SimpleMesonSeqSource(Params(xml_in, path), 7);   // pion_2 =  gamma_4*gamma_5
      }


      //! Construct a0-source a0_2-sink sequential source
      /*!
       * \ingroup hadron
       *
       * \return \f$\gamma_5 * \Gamma(8)^dag * \gamma_5 * F\$
       */
      HadronSeqSource<LatticePropagator>* mesA0A02SeqSrc(XMLReader& xml_in,
							 const std::string& path)
      {
	return new SimpleMesonSeqSource(Params(xml_in, path), 8);   // a0_2 = gamma_4
      }


      //! Construct a0-source rho_x_2-sink sequential source
      /*!
       * \ingroup hadron
       *
       * \return \f$\gamma_5 * \Gamma(9)^dag * \gamma_5 * F\$
       */
      HadronSeqSource<LatticePropagator>* mesA0RhoX2SeqSrc(XMLReader& xml_in,
							   const std::string& path)
      {
	return new SimpleMesonSeqSource(Params(xml_in, path), 9);   // rho_x_2 = gamma_4 * gamma_1
      }


      //! Construct a0-source rho_y_2-sink sequential source
      /*!
       * \ingroup hadron
       *
       * \return \f$\gamma_5 * \Gamma(10)^dag * \gamma_5 * F\$
       */
      HadronSeqSource<LatticePropagator>* mesA0RhoY2SeqSrc(XMLReader& xml_in,
							   const std::string& path)
      {
	return new SimpleMesonSeqSource(Params(xml_in, path), 10);   // rho_y_2 = gamma_4 * gamma_2
      }


      //! Construct a0-source a1_z_1-sink sequential source
      /*!
       * \ingroup hadron
       *
       * \return \f$\gamma_5 * \Gamma(11)^dag * \gamma_5 * F\$
       */
      HadronSeqSource<LatticePropagator>* mesA0A1Z1SeqSrc(XMLReader& xml_in,
							  const std::string& path)
      {
	return new SimpleMesonSeqSource(Params(xml_in, path), 11);   // a1_z_1 = gamma_1 * gamma_2 * gamma_4
      }


      //! Construct a0-source rho_z_2-sink sequential source
      /*!
       * \ingroup hadron
       *
       * \return \f$\gamma_5 * \Gamma(12)^dag * \gamma_5 * F\$
       */
      HadronSeqSource<LatticePropagator>* mesA0RhoZ2SeqSrc(XMLReader& xml_in,
							   const std::string& path)
      {
	return new SimpleMesonSeqSource(Params(xml_in, path), 12);   // rho_z_2 = gamma_4 * gamma_3
      }


      //! Construct a0-source a1_y_1-sink sequential source
      /*!
       * \ingroup hadron
       *
       * \return \f$\gamma_5 * \Gamma(13)^dag * \gamma_5 * F\$
       */
      HadronSeqSource<LatticePropagator>* mesA0A1Y1SeqSrc(XMLReader& xml_in,
							  const std::string& path)
      {
	return new SimpleMesonSeqSource(Params(xml_in, path), 13);   // a1_y_1 = gamma_1 * gamma_3 * gamma_4
      }


      //! Construct a0-source a1_x_1-sink sequential source
      /*!
       * \ingroup hadron
       *
       * \return \f$\gamma_5 * \Gamma(14)^dag * \gamma_5 * F\$
       */
      HadronSeqSource<LatticePropagator>* mesA0A1X1SeqSrc(XMLReader& xml_in,
							  const std::string& path)
      {
	return new SimpleMesonSeqSource(Params(xml_in, path), 14);   // a1_x_1 = gamma_2*gamma_3*gamma_4
      }


      //! Construct pion sequential source
      /*!
       * \ingroup hadron
       *
       * \return \f$\gamma_5 * \Gamma(15)^dag * \gamma_5 * F\$
       */
      HadronSeqSource<LatticePropagator>* mesA0Pion1SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new SimpleMesonSeqSource(Params(xml_in, path), 15);   // pion_1 = gamma_5
      }


      //! Construct pion-source pion-sink
      /*!
       * \ingroup hadron
       *
       * \return \f$\gamma_5 * \Gamma(15)^dag * \gamma_5 * F * \gamma_5\$
       */
      HadronSeqSource<LatticePropagator>* mesPion1Pion1SeqSrc(XMLReader& xml_in,
							      const std::string& path)
      {
	return new PionPionSeqSource(Params(xml_in, path));
      }


      //! Local registration flag
      bool registered = false;

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
				     const multi1d<ForwardProp_t>& forward_headers,
				     const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      QDPIO::cout << "Simple meson sequential source: gamma_sink = " << gamma_sink << endl;
      setTSrce(forward_headers);

      if (quark_propagators.size() != 1)
      {
	QDPIO::cerr << __func__ << ": expect only 1 prop" << endl;
	QDP_abort(1);
      }

      // \f$\gamma_5 * \Gamma(gamma_sink)^dag * \gamma_5 * F\$
//    LatticePropagator tmp = Gamma(G5) * adj(Gamma(gamma_sink)) * Gamma(G5) * quark_propagators[0];

//    int G5 = Ns*Ns-1;
      LatticePropagator tmp = mesA0XSeqSrc(quark_propagators[0], gamma_sink);
      LatticeComplex     ph = conj(phases());
      LatticePropagator fin = project(tmp * ph);

      END_CODE();

      return fin;
    }


    // Compute the 2-pt at the sink
    Complex 
    SimpleMesonSeqSource::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				    const multi1d<ForwardProp_t>& forward_headers,
				    const multi1d<LatticePropagator>& quark_propagators,
				    int gamma_insertion)
    {
      setTSrce(forward_headers);

      if (quark_propagators.size() != 2)
      {
	QDPIO::cerr << __func__ << ": expect 2 props" << endl;
	QDP_abort(1);
      }

      // Construct the meson correlation function
      LatticeComplex corr_fn = 
	trace(gamma5Herm(quark_propagators[1]) * (Gamma(gamma_sink) *
						  quark_propagators[0] * Gamma(gamma_insertion)));

      // Extract the sink at the appropriate momenta
      SftMom sft(0, getTSrce(), getSinkMom(), false, getDecayDir());
      multi2d<DComplex> hsum;
      hsum = sft.sft(corr_fn);

      return hsum[0][getTSink()];
    }


    //! Construct the source
    LatticePropagator
    PionPionSeqSource::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      QDPIO::cout << "Pion-Pion meson sequential source" << endl;
      setTSrce(forward_headers);

      if (quark_propagators.size() != 1)
      {
	QDPIO::cerr << __func__ << ": expect only 1 prop" << endl;
	QDP_abort(1);
      }

      // \f$\gamma_5 * \Gamma(gamma_sink)^dag * \gamma_5 * F * \gamma_5\$
//    LatticePropagator tmp = Gamma(G5) * adj(Gamma(gamma_sink)) * Gamma(G5) * quark_propagators[0];

      int G5 = Ns*Ns-1;
      LatticePropagator tmp = mesA0XSeqSrc(quark_propagators[0], G5) * Gamma(G5);
      LatticeComplex     ph = conj(phases());
      LatticePropagator fin = project(tmp * ph);

      END_CODE();

      return fin;
    }


    // Compute the 2-pt at the sink
    Complex 
    PionPionSeqSource::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<ForwardProp_t>& forward_headers,
				 const multi1d<LatticePropagator>& quark_propagators,
				 int gamma_insertion)
    {
      setTSrce(forward_headers);

      if (quark_propagators.size() != 2)
      {
	QDPIO::cerr << __func__ << ": expect 2 props" << endl;
	QDP_abort(1);
      }

      // Convert gamma_insertion to old convention
      int gamma_i = gamma_insertion ^ unsigned(Ns*Ns-1);
      QDPIO::cout << __func__ << ": gamma_i=" << gamma_i << endl;

      // Construct the meson correlation function
      LatticeComplex corr_fn = 
	trace(gamma5Herm(quark_propagators[1]) * (Gamma(Ns*Ns-1) *
						  quark_propagators[0] * Gamma(gamma_i)));

      // Extract the sink at the appropriate momenta
      SftMom sft(0, getTSrce(), getSinkMom(), false, getDecayDir());
      multi2d<DComplex> hsum;
      hsum = sft.sft(corr_fn);

      return hsum[0][getTSink()];
    }


    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	//! Register all the factories
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a0"), 
										      mesA0A01SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rho_x_1"), 
										      mesA0RhoX1SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rho_y_1"),
										      mesA0RhoY1SeqSrc);
      
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-b1_z"),
										      mesA0B1Z1SeqSrc);
      
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rho_z_1"),
										      mesA0RhoZ1SeqSrc);
      
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-b1_y"),
										      mesA01B1Y1SeqSrc);
      
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-b1_x"),
										      mesA01B1X1SeqSrc);
      
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-pion_2"),
										      mesA01Pion2SeqSrc);
      
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a0_2"),
										      mesA0A02SeqSrc);
      
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rho_x_2"),
										      mesA0RhoX2SeqSrc);
      
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rho_y_2"),
										      mesA0RhoY2SeqSrc);
      
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a1_z"),
										      mesA0A1Z1SeqSrc);
       
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rho_z_2"),
										      mesA0RhoZ2SeqSrc);
       
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a1_y"),
										      mesA0A1Y1SeqSrc);
       
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a1_x"),
										      mesA0A1X1SeqSrc);
       
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-pion_1"),
										      mesA0Pion1SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("pion_1-pion_1"),
										      mesPion1Pion1SeqSrc);

	// keep for historical purposes
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("pion"),
										      mesPion1Pion1SeqSrc);

	registered = true;
      }
      return success;
    }

  }  // end namespace SimpleMesonSeqSourceEnv

}  // end namespace Chroma


  
