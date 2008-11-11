// $Id: deriv_meson_seqsrc_w.cc,v 3.8 2008-11-11 21:27:41 edwards Exp $
/*! \file
 *  \brief Construct meson sequential sources.
 */

#include "meas/hadron/deriv_meson_seqsrc_w.h"
#include "meas/hadron/seqsource_factory_w.h"
#include "meas/smear/displace.h"

#include "util/ferm/gammasgn_w.h"
#include "util/ferm/gamma5_herm_w.h"
#include "util/ferm/symtensor.h"
#include "util/ferm/antisymtensor.h"
#include "util/ferm/etensor.h"
#include "util/ft/sftmom.h"

namespace Chroma 
{

  // Read parameters
  void read(XMLReader& xml, const string& path, DerivMesonSeqSourceEnv::Params& param)
  {
    DerivMesonSeqSourceEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const DerivMesonSeqSourceEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  // Read parameters
  void read(XMLReader& xml, const string& path, DerivMesonSeqSourceEnv::ParamsDir& param)
  {
    DerivMesonSeqSourceEnv::ParamsDir tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const DerivMesonSeqSourceEnv::ParamsDir& param)
  {
    param.writeXML(xml, path);
  }




  //---------------------------------------------------------------------------------
  // Anonymous namespace
  /*! \ingroup hadron */
  namespace
  {
    //! Check only 1 prop passed
    /*! \ingroup hadron */
    void check1Args(const char* name, const multi1d<LatticePropagator>& quark_propagators)
    {
      if (quark_propagators.size() != 1)
      {
	QDPIO::cerr << __func__ << ": expect only 1 prop" << endl;
	QDP_abort(1);
      }
    }

    //! Check only 2 props passed
    /*! \ingroup hadron */
    void check2Args(const char* name, const multi1d<LatticePropagator>& quark_propagators)
    {
      if (quark_propagators.size() != 2)
      {
	QDPIO::cerr << __func__ << ": expect only 2 prop" << endl;
	QDP_abort(1);
      }
    }
  }


  //---------------------------------------------------------------------------------
  //! Meson sequential sources
  /*! \ingroup hadron */
  namespace DerivMesonSeqSourceEnv
  { 
    //! Anonymous namespace
    namespace
    {

      //-------------------- callback functions ---------------------------------------

      //! Construct a0-(pionxNabla_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0PionxNablaT1SeqSrc(XMLReader& xml_in,
								  const std::string& path)
      {
	return new MesA0PionxNablaT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(a0xNabla_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A0xNablaT1SeqSrc(XMLReader& xml_in,
								const std::string& path)
      {
	return new MesA0A0xNablaT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(a0_2xNabla_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A02xNablaT1SeqSrc(XMLReader& xml_in,
								 const std::string& path)
      {
	return new MesA0A02xNablaT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(pion_2xNabla_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0Pion2xNablaT1SeqSrc(XMLReader& xml_in,
								   const std::string& path)
      {
	return new MesA0Pion2xNablaT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(rhoxNabla_A1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0RhoxNablaA1SeqSrc(XMLReader& xml_in,
								 const std::string& path)
      {
	return new MesA0RhoxNablaA1SeqSrc(Params(xml_in, path));
      }

      //! Construct a0-(rhoxNabla_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0RhoxNablaT1SeqSrc(XMLReader& xml_in,
								 const std::string& path)
      {
	return new MesA0RhoxNablaT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(rhoxNabla_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0RhoxNablaT2SeqSrc(XMLReader& xml_in,
								 const std::string& path)
      {
	return new MesA0RhoxNablaT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(rhoxNabla_E) sequential source
      HadronSeqSource<LatticePropagator>* mesA0RhoxNablaESeqSrc(XMLReader& xml_in,
								const std::string& path)
      {
	return new MesA0RhoxNablaESeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(rho_2xNabla_A1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0Rho2xNablaA1SeqSrc(XMLReader& xml_in,
								  const std::string& path)
      {
	return new MesA0Rho2xNablaA1SeqSrc(Params(xml_in, path));
      }

      //! Construct a0-(rho_2xNabla_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0Rho2xNablaT1SeqSrc(XMLReader& xml_in,
								  const std::string& path)
      {
	return new MesA0Rho2xNablaT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(rho_2xNabla_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0Rho2xNablaT2SeqSrc(XMLReader& xml_in,
								  const std::string& path)
      {
	return new MesA0Rho2xNablaT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(rho_2xNabla_E) sequential source
      HadronSeqSource<LatticePropagator>* mesA0Rho2xNablaESeqSrc(XMLReader& xml_in,
								 const std::string& path)
      {
	return new MesA0Rho2xNablaESeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(a1xNabla_A1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A1xNablaA1SeqSrc(XMLReader& xml_in,
								const std::string& path)
      {
	return new MesA0A1xNablaA1SeqSrc(Params(xml_in, path));
      }

      //! Construct a0-(a1xNabla_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A1xNablaT1SeqSrc(XMLReader& xml_in,
								const std::string& path)
      {
	return new MesA0A1xNablaT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(a1xNabla_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A1xNablaT2SeqSrc(XMLReader& xml_in,
								const std::string& path)
      {
	return new MesA0A1xNablaT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(a1xNabla_E) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A1xNablaESeqSrc(XMLReader& xml_in,
							       const std::string& path)
      {
	return new MesA0A1xNablaESeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(b1xNabla_A1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0B1xNablaA1SeqSrc(XMLReader& xml_in,
								const std::string& path)
      {
	return new MesA0B1xNablaA1SeqSrc(Params(xml_in, path));
      }

      //! Construct a0-(b1xNabla_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0B1xNablaT1SeqSrc(XMLReader& xml_in,
								const std::string& path)
      {
	return new MesA0B1xNablaT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(b1xNabla_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0B1xNablaT2SeqSrc(XMLReader& xml_in,
								const std::string& path)
      {
	return new MesA0B1xNablaT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(b1xNabla_E) sequential source
      HadronSeqSource<LatticePropagator>* mesA0B1xNablaESeqSrc(XMLReader& xml_in,
								const std::string& path)
      {
	return new MesA0B1xNablaESeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(pionxD_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0PionxDT2SeqSrc(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesA0PionxDT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(a0xD_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A0xDT2SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new MesA0A0xDT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(a0_2xD_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A02xDT2SeqSrc(XMLReader& xml_in,
							     const std::string& path)
      {
	return new MesA0A02xDT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(pion_2xD_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0Pion2xDT2SeqSrc(XMLReader& xml_in,
							       const std::string& path)
      {
	return new MesA0Pion2xDT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(rhoxD_A2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0RhoxDA2SeqSrc(XMLReader& xml_in,
							     const std::string& path)
      {
	return new MesA0RhoxDA2SeqSrc(Params(xml_in, path));
      }

      //! Construct a0-(rhoxD_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0RhoxDT1SeqSrc(XMLReader& xml_in,
							     const std::string& path)
      {
	return new MesA0RhoxDT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(rhoxD_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0RhoxDT2SeqSrc(XMLReader& xml_in,
							     const std::string& path)
      {
	return new MesA0RhoxDT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(rhoxD_E) sequential source
      HadronSeqSource<LatticePropagator>* mesA0RhoxDESeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new MesA0RhoxDESeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(rho_2xD_A2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0Rho2xDA2SeqSrc(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesA0Rho2xDA2SeqSrc(Params(xml_in, path));
      }

      //! Construct a0-(rho_2xD_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0Rho2xDT1SeqSrc(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesA0Rho2xDT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(rho_2xD_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0Rho2xDT2SeqSrc(XMLReader& xml_in,
							     const std::string& path)
      {
	return new MesA0Rho2xDT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(rho_2xD_E) sequential source
      HadronSeqSource<LatticePropagator>* mesA0Rho2xDESeqSrc(XMLReader& xml_in,
							     const std::string& path)
      {
	return new MesA0Rho2xDESeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(a1xD_A2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A1xDA2SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new MesA0A1xDA2SeqSrc(Params(xml_in, path));
      }

      //! Construct a0-(a1xD_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A1xDT1SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new MesA0A1xDT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(a1xD_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A1xDT2SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new MesA0A1xDT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(a1xD_E) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A1xDESeqSrc(XMLReader& xml_in,
							   const std::string& path)
      {
	return new MesA0A1xDESeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(b1xD_A2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0B1xDA2SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new MesA0B1xDA2SeqSrc(Params(xml_in, path));
      }

      //! Construct a0-(b1xD_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0B1xDT1SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new MesA0B1xDT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(b1xD_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0B1xDT2SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new MesA0B1xDT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(b1xD_E) sequential source
      HadronSeqSource<LatticePropagator>* mesA0B1xDESeqSrc(XMLReader& xml_in,
							   const std::string& path)
      {
	return new MesA0B1xDESeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(pionxB_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0PionxBT1SeqSrc(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesA0PionxBT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(a0xB_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A0xBT1SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new MesA0A0xBT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(a0_2xB_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A02xBT1SeqSrc(XMLReader& xml_in,
							     const std::string& path)
      {
	return new MesA0A02xBT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(pion_2xB_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0Pion2xBT1SeqSrc(XMLReader& xml_in,
							       const std::string& path)
      {
	return new MesA0Pion2xBT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(rhoxB_A1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0RhoxBA1SeqSrc(XMLReader& xml_in,
							     const std::string& path)
      {
	return new MesA0RhoxBA1SeqSrc(Params(xml_in, path));
      }

      //! Construct a0-(rhoxB_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0RhoxBT1SeqSrc(XMLReader& xml_in,
							     const std::string& path)
      {
	return new MesA0RhoxBT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(rhoxB_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0RhoxBT2SeqSrc(XMLReader& xml_in,
							     const std::string& path)
      {
	return new MesA0RhoxBT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(rhoxB_E) sequential source
      HadronSeqSource<LatticePropagator>* mesA0RhoxBESeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new MesA0RhoxBESeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(rho_2xB_A1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0Rho2xBA1SeqSrc(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesA0Rho2xBA1SeqSrc(Params(xml_in, path));
      }

      //! Construct a0-(rho_2xB_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0Rho2xBT1SeqSrc(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesA0Rho2xBT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(rho_2xB_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0Rho2xBT2SeqSrc(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesA0Rho2xBT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(rho_2xB_E) sequential source
      HadronSeqSource<LatticePropagator>* mesA0Rho2xBESeqSrc(XMLReader& xml_in,
							     const std::string& path)
      {
	return new MesA0Rho2xBESeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(a1xB_A1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A1xBA1SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new MesA0A1xBA1SeqSrc(Params(xml_in, path));
      }

      //! Construct a0-(a1xB_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A1xBT1SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new MesA0A1xBT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(a1xB_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A1xBT2SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new MesA0A1xBT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(a1xB_E) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A1xBESeqSrc(XMLReader& xml_in,
							   const std::string& path)
      {
	return new MesA0A1xBESeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(b1xB_A1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0B1xBA1SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new MesA0B1xBA1SeqSrc(Params(xml_in, path));
      }

      //! Construct a0-(b1xB_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0B1xBT1SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new MesA0B1xBT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(b1xB_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0B1xBT2SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new MesA0B1xBT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(b1xB_E) sequential source
      HadronSeqSource<LatticePropagator>* mesA0B1xBESeqSrc(XMLReader& xml_in,
							   const std::string& path)
      {
	return new MesA0B1xBESeqSrc(ParamsDir(xml_in, path));
      }

    } // end anonymous namespace


    //---------------------------------------------------------------------------------
    //! Initialize
    Params::Params()
    {
      deriv_length = 0;
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

      read(paramtop, "deriv_length", deriv_length);
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

      write(xml, "deriv_length", deriv_length);
      write(xml, "j_decay", j_decay);
      write(xml, "t_sink", t_sink);
      write(xml, "sink_mom", sink_mom);
      pop(xml);
    }


    //! Initialize
    ParamsDir::ParamsDir()
    {
      deriv_dir = -1;
      deriv_length = 0;
      j_decay = -1;
      t_sink  = -1;
      sink_mom.resize(Nd-1);
      sink_mom = 0;
    }


    //! Read parameters
    ParamsDir::ParamsDir(XMLReader& xml, const string& path)
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

      read(paramtop, "deriv_dir", deriv_dir);
      read(paramtop, "deriv_length", deriv_length);
      read(paramtop, "j_decay", j_decay);
      read(paramtop, "t_sink", t_sink);
      read(paramtop, "sink_mom", sink_mom);
    }


    // Writer
    void ParamsDir::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);

      write(xml, "deriv_dir", deriv_dir);
      write(xml, "deriv_length", deriv_length);
      write(xml, "j_decay", j_decay);
      write(xml, "t_sink", t_sink);
      write(xml, "sink_mom", sink_mom);
      pop(xml);
    }

 
    //---------------------------------------------------------------------------------
    //! Apply right nabla
    LatticePropagator 
    DerivMesonSeqSourceBase::nabla(const LatticePropagator& F, 
				   const multi1d<LatticeColorMatrix>& u,
				   int mu) const
    {
      return Chroma::rightNabla(F, u, mu, getDerivLength());
    }


    //! Apply right D
    LatticePropagator 
    DerivMesonSeqSourceBase::D(const LatticePropagator& F, 
			       const multi1d<LatticeColorMatrix>& u,
			       int mu) const
    {
      return Chroma::rightD(F, u, mu, getDerivLength());
    }


    //! Apply right nabla
    LatticePropagator 
    DerivMesonSeqSourceBase::B(const LatticePropagator& F, 
			       const multi1d<LatticeColorMatrix>& u,
			       int mu) const
    {
      return Chroma::rightB(F, u, mu, getDerivLength());
    }


    //---------------------------------------------------------------------------------
    //! Apply left and right "nabla_i" onto the source
    /*!
     * \f$\nabla_\mu f(x) = U_\mu(x)f(x+\mu) - U_{-\mu}(x)f(x-\mu)\f$
     *
     * \return $\f \nabla_\mu Fx,0) = U_\mu(x) F(x+\mu)  - U_{x-\mu}^\dag(x-\mu) F(x-\mu)\f$
     */
    LatticePropagator 
    DerivMesonSeqSourceBase::threePtNabla(const LatticePropagator& F, 
					  const multi1d<LatticeColorMatrix>& u,
					  int mu) const
    {
      LatticeComplex ph = conj(phases());
      return ph*nabla(F, u, mu) + nabla(ph*F, u, mu);
    }

    //! Apply left and right "D_i" operator onto source
    /*!
     * \f$D_i = s_{ijk}\nabla_j\nabla_k\f$
     *
     * where  \f$s_{ijk} = +1 \quad\forall i\ne j, j\ne k, i \ne k\f$
     * 
     * \return $\f p*D_\mu F(x,0) + D_\mu (p*F(x,0)) + 2 \sum_{j,k}s_{ijk} \nabla_j (p*\nabla_k F(x,0))\f$
     */
    LatticePropagator 
    DerivMesonSeqSourceBase::threePtD(const LatticePropagator& F,
				      const multi1d<LatticeColorMatrix>& u,
				      int mu) const
    {
      LatticePropagator tmp = zero;
      LatticeComplex ph = conj(phases());

      // Slow implementation - to speed up could compute once the \nabla_j deriv
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (symTensor3d(mu,j,k) != 0)
	  {
	    tmp += ph * nabla(nabla(F,u,j), u, k);
	    tmp += nabla(nabla(ph * F,u,j), u, k);
	    tmp += Real(2)*nabla(ph * nabla(F,u,j), u, k);
	  }
	}

      return tmp;
    }


    //! Apply left and right "B_i" operator onto source
    /*!
     * \f$B_i = \epsilon_{ijk}\nabla_j\nabla_k\f$
     *
     * \return $\f -B_\mu F(x,0)\f$
     */
    LatticePropagator 
    DerivMesonSeqSourceBase::threePtB(const LatticePropagator& F,
				      const multi1d<LatticeColorMatrix>& u,
				      int mu) const
    {
      LatticePropagator tmp = zero;
      LatticeComplex ph = conj(phases());

      // Slow implementation - to speed up could compute once the \nabla_j deriv
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (antiSymTensor3d(mu,j,k) != 0)
	    tmp -= Real(antiSymTensor3d(mu,j,k)) * nabla(nabla(ph*F,u,j), u, k);
	}

      return tmp;
    }


    //---------------------------------------------------------------------------------
    //! Vector from a Nabla operator
    multi1d<LatticePropagator> 
    DerivMesonSeqSourceBase::threePtNablaVector(const LatticePropagator& F,
						const multi1d<LatticeColorMatrix>& u) const
    {
      multi1d<LatticePropagator> out(3);
	
      for(int k=0; k < 3; ++k)
	out[k] = threePtNabla(F,u,k);
	
      return out;
    }


    //! Vector from a D operator
    multi1d<LatticePropagator> 
    DerivMesonSeqSourceBase::threePtDVector(const LatticePropagator& F,
					    const multi1d<LatticeColorMatrix>& u) const
    {
      multi1d<LatticePropagator> out(3);
	
      for(int k=0; k < 3; ++k)
	out[k] = threePtD(F,u,k);
	
      return out;
    }


    //! Vector from a B operator
    multi1d<LatticePropagator> 
    DerivMesonSeqSourceBase::threePtBVector(const LatticePropagator& F,
					    const multi1d<LatticeColorMatrix>& u) const
    {
      multi1d<LatticePropagator> out(3);
	
      for(int k=0; k < 3; ++k)
	out[k] = threePtB(F,u,k);
	
      return out;
    }


    //---------------------------------------------------------------------------------
    // Apply left and right "nabla_i" onto the source
    LatticeComplex
    DerivMesonSeqSourceBase::twoPtNabla(const LatticePropagator& F, 
					const multi1d<LatticeColorMatrix>& u,
					int mu, int g, int insertion) const
    {
      LatticePropagator drv = nabla(F,u,mu);
      LatticeComplex corr_fn = 
	  trace(gamma5Herm(drv) * Gamma(g) * F   * Gamma(insertion))
	- trace(gamma5Herm(F)   * Gamma(g) * drv * Gamma(insertion));

      return corr_fn;
    }


    // Apply left and right "D_i" operator onto source
    LatticeComplex
    DerivMesonSeqSourceBase::twoPtD(const LatticePropagator& F,
				    const multi1d<LatticeColorMatrix>& u,
				    int mu, int g, int insertion) const
    {
      LatticePropagator ddd = D(F,u,mu);

      LatticeComplex corr_fn = 
	  trace(gamma5Herm(ddd)  * Gamma(g) * F   * Gamma(insertion))
	+ trace(gamma5Herm(F)    * Gamma(g) * ddd * Gamma(insertion));

      // Slow implementation - to speed up could compute once the \nabla_j deriv
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (symTensor3d(mu,j,k) != 0)
	    corr_fn -= Real(2)*trace(gamma5Herm(nabla(F,u,j)) * Gamma(g) * 
				     nabla(F,u,k) * Gamma(insertion));
	}

      return corr_fn;
    }


    // Apply left and right "B_i" operator onto source
    LatticeComplex
    DerivMesonSeqSourceBase::twoPtB(const LatticePropagator& F,
				    const multi1d<LatticeColorMatrix>& u,
				    int mu, int g, int insertion) const
    {
      return trace(gamma5Herm(F) * Gamma(g) * B(F,u,mu) * Gamma(insertion));
    }


    // Project onto the fixed sink-momentum and return the 2-pt at the sink
    Complex 
    DerivMesonSeqSourceBase::momentumProject(const LatticeComplex& corr_fn) const
    {
      // Extract the sink at the appropriate momenta
      SftMom sft(0, getTSrce(), getSinkMom(), false, getDecayDir());
      multi2d<DComplex> hsum;
      hsum = sft.sft(corr_fn);

      END_CODE();

      return hsum[0][getTSink()];
    }


    //---------------------------------------------------------------------------------
    // All 3-pt seqsrc are variants on these basic objects
    namespace
    {
      int G5 = Ns*Ns-1;

      //---------------------------------------------------------------------------------
      //! a0xVector
      /*! \return \f$V_i\f$ */
      LatticePropagator a0xVector(const LatticePropagator& F)
      {
	return F;
      }


      //! pionxVector
      /*! \return \f$gamma_5 V_i\f$ */
      LatticePropagator pionxVector(const LatticePropagator& F)
      {
	return Gamma(G5) * F;
      }


      //! pion_2xVector
      /*! \return \f$\gamma_4\gamma_5 V_i\f$ */
      LatticePropagator pion_2xVector(const LatticePropagator& F)
      {
	return Gamma(1<<3) * (Gamma(G5) * F);
      }


      //! a0_2xVector
      /*! \return \f$\gamma_4 V_i\f$ */
      LatticePropagator a0_2xVector(const LatticePropagator& F)
      {
	return Gamma(1<<3) * F;
      }


      //---------------------------------------------------------------------------------
      //! rhoxVector_sum
      /*! \return \f$\gamma_i V_i\f$ */
      LatticePropagator rhoxVector_sum(const multi1d<LatticePropagator>& F)
      {
	LatticePropagator fin = zero;
	
	for(int k=0; k < 3; ++k)
	  fin += Gamma(1 << k) * F[k];
	
	return fin;
      }


      //! rhoxVector_antisym
      /*! \return \f$\epsilon_{ijk}\gamma_j V_k\f$ */
      LatticePropagator rhoxVector_antisym(const multi1d<LatticePropagator>& F, int dir)
      {
	LatticePropagator fin = zero;

	for(int j=0; j < 3; ++j)
	  for(int k=0; k < 3; ++k)
	  {
	    int a = antiSymTensor3d(dir,j,k);
	    if (a != 0)
	      fin += Real(a) * (Gamma(1 << j) * F[k]);
	  }

	return fin;
      }


      //! rhoxVector_sym
      /*! \return \f$s_{ijk}\gamma_j V_k\f$ */
      LatticePropagator rhoxVector_sym(const multi1d<LatticePropagator>& F, int dir)
      {
	LatticePropagator fin = zero;

	for(int j=0; j < 3; ++j)
	  for(int k=0; k < 3; ++k)
	  {
	    int s = symTensor3d(dir,j,k);
	    if (s != 0)
	      fin += Real(s) * (Gamma(1 << j) * F[k]);
	  }

	return fin;
      }


      //! rhoxVector_E
      /*! \return \f$Q_{\alpha jk}\gamma_j V_k\f$ */
      LatticePropagator rhoxVector_E(const multi1d<LatticePropagator>& F, int dir)
      {
	LatticePropagator fin = zero;

	for(int j=0; j < 3; ++j)
	  for(int k=0; k < 3; ++k)
	  {
	    Real e = ETensor3d(dir,j,k);
	    if (toBool(e != 0.0))
	      fin += e * (Gamma(1 << j) * F[k]);
	  }
      
	return fin;
      }


      //---------------------------------------------------------------------------------
      //! rho_2xVector_sum
      /*! \return \f$\gamma_4\gamma_i V_i\f$ */
      LatticePropagator rho_2xVector_sum(const multi1d<LatticePropagator>& F)
      {
	return Gamma(1<<3) * rhoxVector_sum(F);
      }


      //! rho_2xVector_antisym
      /*! \return \f$\epsilon_{ijk}\gamma_4\gamma_j V_k\f$ */
      LatticePropagator rho_2xVector_antisym(const multi1d<LatticePropagator>& F, int dir)
      {
	return Gamma(1<<3) * rhoxVector_antisym(F,dir);
      }


      //! rho_2xVector_sym
      /*! \return \f$s_{ijk}\gamma_4\gamma_j V_k\f$ */
      LatticePropagator rho_2xVector_sym(const multi1d<LatticePropagator>& F, int dir)
      {
	return Gamma(1<<3) * rhoxVector_sym(F,dir);
      }


      //! rho_2xVector_E
      /*! \return \f$Q_{\alpha jk}\gamma_4\gamma_j V_k\f$ */
      LatticePropagator rho_2xVector_E(const multi1d<LatticePropagator>& F, int dir)
      {
	return Gamma(1<<3) * rhoxVector_E(F,dir);
      }


      //---------------------------------------------------------------------------------
      //! a1xVector_sum
      /*! \return \f$\gamma_5\gamma_i V_i\f$ */
      LatticePropagator a1xVector_sum(const multi1d<LatticePropagator>& F)
      {
	return Gamma(G5) * rhoxVector_sum(F);
      }


      //! rho_2xVector_antisym
      /*! \return \f$\epsilon_{ijk}\gamma_5\gamma_j V_k\f$ */
      LatticePropagator a1xVector_antisym(const multi1d<LatticePropagator>& F,
						  int dir)
      {
	return Gamma(G5) * rhoxVector_antisym(F,dir);
      }


      //! rho_2xVector_sym
      /*! \return \f$s_{ijk}\gamma_5\gamma_j V_k\f$ */
      LatticePropagator a1xVector_sym(const multi1d<LatticePropagator>& F,
					      int dir)
      {
	return Gamma(G5) * rhoxVector_sym(F,dir);
      }


      //! rho_2xVector_E
      /*! \return \f$Q_{\alpha jk}\gamma_5\gamma_j V_k\f$ */
      LatticePropagator a1xVector_E(const multi1d<LatticePropagator>& F,
					    int dir)
      {
	return Gamma(G5) * rhoxVector_E(F,dir);
      }


      //---------------------------------------------------------------------------------
      //! b1xVector_sum
      /*! \return \f$\gamma_4\gamma_5\gamma_i V_i\f$ */
      LatticePropagator b1xVector_sum(const multi1d<LatticePropagator>& F)
      {
	return Gamma(1<<3) * (Gamma(G5) * rhoxVector_sum(F));
      }


      //! b1xVector_antisym
      /*! \return \f$\epsilon_{ijk}\gamma_4\gamma_5\gamma_j V_k\f$ */
      LatticePropagator b1xVector_antisym(const multi1d<LatticePropagator>& F,
						  int dir)
      {
	return Gamma(1<<3) * (Gamma(G5) * rhoxVector_antisym(F,dir));
      }


      //! b1xVector_sym
      /*! \return \f$s_{ijk}\gamma_4\gamma_5\gamma_j V_k\f$ */
      LatticePropagator b1xVector_sym(const multi1d<LatticePropagator>& F,
					      int dir)
      {
	return Gamma(1<<3) * (Gamma(G5) * rhoxVector_sym(F,dir));
      }


      //! b1xVector_E
      /*! \return \f$Q_{\alpha jk}\gamma_4\gamma_5\gamma_j V_k\f$ */
      LatticePropagator b1xVector_E(const multi1d<LatticePropagator>& F,
					    int dir)
      {
	return Gamma(1<<3) * (Gamma(G5) * rhoxVector_E(F,dir));
      }


    } // namespace anonymous


    //---------------------------------------------------------------------------------
    // Construct a0-(a0xNabla_T1) sequential source
    // See corresponding .h file for doxygen comments
    LatticePropagator
    MesA0A0xNablaT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<ForwardProp_t>& forward_headers,
				      const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0A0xNablaT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \nabla_i\f$
      return project(a0xVector(threePtNabla(quark_propagators[0],u,getDerivDir())));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0A0xNablaT1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				     const multi1d<ForwardProp_t>& forward_headers,
				     const multi1d<LatticePropagator>& quark_propagators,
				     int insertion)
    {
      START_CODE();

      check1Args("MesA0A0xNablaT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \nabla_i\f$
      LatticeComplex corr_fn = twoPtNabla(tmp,u,getDerivDir(),0,insertion);
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    // Construct a0-(pionxNabla_T1) sequential source
    // See corresponding .h file for doxygen comments
    LatticePropagator
    MesA0PionxNablaT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
					const multi1d<ForwardProp_t>& forward_headers,
					const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0PionxNablaT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_5\nabla_i\f$
      return project(pionxVector(threePtNabla(quark_propagators[0],u,getDerivDir())));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0PionxNablaT1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				       const multi1d<ForwardProp_t>& forward_headers,
				       const multi1d<LatticePropagator>& quark_propagators,
				       int insertion)
    {
      START_CODE();

      check1Args("MesA0PionxNablaT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator fin;
      const int G5 = Ns*Ns-1;

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5\nabla_i\f$
      LatticeComplex corr_fn = twoPtNabla(tmp,u,getDerivDir(),G5,insertion);
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    // Construct a0-(pion_2xNabla_T1) sequential source
    // See corresponding .h file for doxygen comments
    LatticePropagator
    MesA0Pion2xNablaT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
					 const multi1d<ForwardProp_t>& forward_headers,
					 const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0Pion2xNablaT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_4 \gamma_5 \nabla_i\f$
      return project(pion_2xVector(threePtNabla(quark_propagators[0],u,getDerivDir())));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0Pion2xNablaT1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
					const multi1d<ForwardProp_t>& forward_headers,
					const multi1d<LatticePropagator>& quark_propagators,
					int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn;
      const int G5 = Ns*Ns-1;

      check1Args("MesA0Pion2xNablaT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];
      
      // \f$\Gamma_f \equiv \gamma_4 \gamma_5 \nabla_i\f$
      corr_fn = gammaSgn(1<<3,G5) * twoPtNabla(tmp,u,getDerivDir(), (1<<3)^G5, insertion);
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    // Construct a0-(a0_2xNabla_T1) sequential source
    // See corresponding .h file for doxygen comments
    LatticePropagator
    MesA0A02xNablaT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				       const multi1d<ForwardProp_t>& forward_headers,
				       const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0A02xNablaT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);
      
      // \f$\Gamma_f \equiv \gamma_4 \nabla_i\f$
      return project(-a0_2xVector(threePtNabla(quark_propagators[0],u,getDerivDir())));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0A02xNablaT1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<ForwardProp_t>& forward_headers,
				      const multi1d<LatticePropagator>& quark_propagators,
				      int insertion)
    {
      START_CODE();

      check1Args("MesA0A02xNablaT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];
      
      // \f$\Gamma_f \equiv \gamma_4 \nabla_i\f$
      LatticeComplex corr_fn = twoPtNabla(tmp,u,getDerivDir(), 1<<3, insertion);
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //---------------------------------------------------------------------------------
    // Construct a0-(rhoxNabla_A1) sequential source
    LatticePropagator
    MesA0RhoxNablaA1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				       const multi1d<ForwardProp_t>& forward_headers,
				       const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0RhoxNablaA1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_i\nabla_i\f$
      return project(-rhoxVector_sum(threePtNablaVector(quark_propagators[0],u)));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0RhoxNablaA1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<ForwardProp_t>& forward_headers,
				      const multi1d<LatticePropagator>& quark_propagators,
				      int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;

      check1Args("MesA0RhoxNablaA1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_i\nabla_i\f$
      for(int k=0; k < 3; ++k)
	corr_fn += twoPtNabla(tmp,u,k, 1<<k, insertion);
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    // Construct a0-(rhoxNabla_T1) sequential source
    LatticePropagator
    MesA0RhoxNablaT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				       const multi1d<ForwardProp_t>& forward_headers,
				       const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0RhoxNablaT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j \nabla_k\f$  
      return project(-rhoxVector_antisym(threePtNablaVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0RhoxNablaT1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<ForwardProp_t>& forward_headers,
				      const multi1d<LatticePropagator>& quark_propagators,
				      int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;

      check1Args("MesA0RhoxNablaT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int a = antiSymTensor3d(getDerivDir(),j,k);
	  if (a != 0)
	    corr_fn += Real(a) * twoPtNabla(tmp,u,k, 1<<j, insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    // Construct a0-(rhoxNabla_T2) sequential source
    LatticePropagator
    MesA0RhoxNablaT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				       const multi1d<ForwardProp_t>& forward_headers,
				       const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0RhoxNablaT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv s_{ijk}\gamma_j \nabla_k\f$  
      return project(-rhoxVector_sym(threePtNablaVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0RhoxNablaT2SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<ForwardProp_t>& forward_headers,
				      const multi1d<LatticePropagator>& quark_propagators,
				      int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;

      check1Args("MesA0RhoxNablaT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv s_{ijk}\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int s = symTensor3d(getDerivDir(),j,k);
	  if (s != 0)
	    corr_fn += Real(s) * twoPtNabla(tmp,u,k, 1<<j, insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    // Construct a0-(rhoxNabla_E) sequential source
    LatticePropagator
    MesA0RhoxNablaESeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<ForwardProp_t>& forward_headers,
				      const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0RhoxNablaESeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv Q_{\alpha jk}\gamma_j \nabla_k\f$  
      return project(-rhoxVector_E(threePtNablaVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0RhoxNablaESeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				     const multi1d<ForwardProp_t>& forward_headers,
				     const multi1d<LatticePropagator>& quark_propagators,
				     int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;

      check1Args("MesA0RhoxNablaESeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv Q_{\alpha jk}\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  Real e = ETensor3d(getDerivDir(),j,k);
	  if (toBool(e != 0.0))
	    corr_fn += e * twoPtNabla(tmp,u,k, 1<<j, insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //---------------------------------------------------------------------------------
    // Construct a0-(rho_2xNabla_A1) sequential source
    LatticePropagator
    MesA0Rho2xNablaA1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				       const multi1d<ForwardProp_t>& forward_headers,
				       const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0Rho2xNablaA1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_4\gamma_i\nabla_i\f$
      return project(-rho_2xVector_sum(threePtNablaVector(quark_propagators[0],u)));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0Rho2xNablaA1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<ForwardProp_t>& forward_headers,
				      const multi1d<LatticePropagator>& quark_propagators,
				      int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;

      check1Args("MesA0Rho2xNablaA1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_4\gamma_i\nabla_i\f$
      for(int k=0; k < 3; ++k)
	corr_fn += gammaSgn(1<<3,1<<k) * twoPtNabla(tmp,u,k, (1<<3)^(1<<k), insertion);
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    // Construct a0-(rho_2xNabla_T1) sequential source
    LatticePropagator
    MesA0Rho2xNablaT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				       const multi1d<ForwardProp_t>& forward_headers,
				       const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0Rho2xNablaT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_4\gamma_j \nabla_k\f$  
      return project(-rho_2xVector_antisym(threePtNablaVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0Rho2xNablaT1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<ForwardProp_t>& forward_headers,
				      const multi1d<LatticePropagator>& quark_propagators,
				      int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;

      check1Args("MesA0Rho2xNablaT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_4\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int a = antiSymTensor3d(getDerivDir(),j,k);
	  if (a != 0)
	    corr_fn += Real(a) * gammaSgn(1<<3,1<<j) * twoPtNabla(tmp,u,k, (1<<3)^(1<<j), insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    // Construct a0-(rho_2xNabla_T2) sequential source
    LatticePropagator
    MesA0Rho2xNablaT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				       const multi1d<ForwardProp_t>& forward_headers,
				       const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0Rho2xNablaT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv s_{ijk}\gamma_4\gamma_j \nabla_k\f$  
      return project(-rho_2xVector_sym(threePtNablaVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0Rho2xNablaT2SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<ForwardProp_t>& forward_headers,
				      const multi1d<LatticePropagator>& quark_propagators,
				      int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;

      check1Args("MesA0Rho2xNablaT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv s_{ijk}\gamma_4\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int s = symTensor3d(getDerivDir(),j,k);
	  if (s != 0)
	    corr_fn += Real(s) * gammaSgn(1<<3,1<<j) * twoPtNabla(tmp,u,k, (1<<3)^(1<<j), insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    // Construct a0-(rho_2xNabla_E) sequential source
    LatticePropagator
    MesA0Rho2xNablaESeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<ForwardProp_t>& forward_headers,
				      const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0Rho2xNablaESeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv Q_{\alpha jk}\gamma_4\gamma_j \nabla_k\f$  
      return project(-rho_2xVector_E(threePtNablaVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0Rho2xNablaESeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<ForwardProp_t>& forward_headers,
				      const multi1d<LatticePropagator>& quark_propagators,
				      int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;

      check1Args("MesA0Rho2xNablaESeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv Q_{\alpha jk}\gamma_4\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  Real e = ETensor3d(getDerivDir(),j,k);
	  if (toBool(e != 0.0))
	    corr_fn += e * gammaSgn(1<<3,1<<j) * twoPtNabla(tmp,u,k, (1<<3)^(1<<j), insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //---------------------------------------------------------------------------------
    // Construct a0-(a1xNabla_A1) sequential source
    LatticePropagator
    MesA0A1xNablaA1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<ForwardProp_t>& forward_headers,
				      const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0A1xNablaA1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_5\gamma_i \nabla_i\f$  
      return project(a1xVector_sum(threePtNablaVector(quark_propagators[0],u)));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0A1xNablaA1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				     const multi1d<ForwardProp_t>& forward_headers,
				     const multi1d<LatticePropagator>& quark_propagators,
				     int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0A1xNablaA1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5\gamma_i \nabla_i\f$  
      for(int k=0; k < 3; ++k)
	corr_fn += gammaSgn(G5,1<<k) * twoPtNabla(tmp,u,k, G5^(1<<k), insertion);
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    // Construct a0-(a1xNabla_T1) sequential source
    LatticePropagator
    MesA0A1xNablaT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<ForwardProp_t>& forward_headers,
				      const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0A1xNablaT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_5 \epsilon_{ijk}\gamma_j \nabla_k\f$  
      return project(a1xVector_antisym(threePtNablaVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0A1xNablaT1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				     const multi1d<ForwardProp_t>& forward_headers,
				     const multi1d<LatticePropagator>& quark_propagators,
				     int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0A1xNablaT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5 \epsilon_{ijk}\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int a = antiSymTensor3d(getDerivDir(),j,k);
	  if (a != 0)
	    corr_fn += Real(a) * gammaSgn(G5,1<<j) * twoPtNabla(tmp,u,k,G5^(1<<j),insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    // Construct a0-(a1xNabla_T2) sequential source
    LatticePropagator
    MesA0A1xNablaT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<ForwardProp_t>& forward_headers,
				      const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0A1xNablaT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j \nabla_k\f$  
      return project(a1xVector_sym(threePtNablaVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0A1xNablaT2SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				     const multi1d<ForwardProp_t>& forward_headers,
				     const multi1d<LatticePropagator>& quark_propagators,
				     int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0A1xNablaT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int s = symTensor3d(getDerivDir(),j,k);
	  if (s != 0)
	    corr_fn += Real(s) * gammaSgn(G5,1<<j) * twoPtNabla(tmp,u,k,G5^(1<<j),insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    // Construct a0-(a1xNabla_E) sequential source
    LatticePropagator
    MesA0A1xNablaESeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				     const multi1d<ForwardProp_t>& forward_headers,
				     const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0A1xNablaESeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_5 Q_{\alpha jk}\gamma_j \nabla_k\f$  
      return project(a1xVector_E(threePtNablaVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0A1xNablaESeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				    const multi1d<ForwardProp_t>& forward_headers,
				    const multi1d<LatticePropagator>& quark_propagators,
				    int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0A1xNablaESeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5 Q_{\alpha jk}\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  Real e = ETensor3d(getDerivDir(),j,k);
	  if (toBool(e != 0.0))
	    corr_fn += e * gammaSgn(G5,1<<j) * twoPtNabla(tmp,u,k,G5^(1<<j),insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //---------------------------------------------------------------------------------
    // Construct a0-(b1xNabla_A1) sequential source
    LatticePropagator
    MesA0B1xNablaA1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<ForwardProp_t>& forward_headers,
				      const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0B1xNablaA1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 \gamma_i \nabla_i\f$  
      return project(-b1xVector_sum(threePtNablaVector(quark_propagators[0],u)));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0B1xNablaA1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				     const multi1d<ForwardProp_t>& forward_headers,
				     const multi1d<LatticePropagator>& quark_propagators,
				     int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0B1xNablaA1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 \gamma_i \nabla_i\f$  
      for(int k=0; k < 3; ++k)
	corr_fn += gammaSgn(1<<3,G5) * gammaSgn((1<<3)^G5,1<<k) * 
	  twoPtNabla(tmp,u,k,(1<<3)^G5^(1<<k),insertion);
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    // Construct a0-(b1xNabla_T1) sequential source
    LatticePropagator
    MesA0B1xNablaT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<ForwardProp_t>& forward_headers,
				      const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0B1xNablaT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_4\gamma_5\epsilon_{ijk}\gamma_j \nabla_k\f$  
      return project(-b1xVector_antisym(threePtNablaVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0B1xNablaT1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				     const multi1d<ForwardProp_t>& forward_headers,
				     const multi1d<LatticePropagator>& quark_propagators,
				     int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0B1xNablaT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_4\gamma_5\epsilon_{ijk}\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int a = antiSymTensor3d(getDerivDir(),j,k);
	  if (a != 0)
	    corr_fn += Real(a) * gammaSgn(1<<3,G5) * gammaSgn((1<<3)^G5,1<<j) * 
	      twoPtNabla(tmp,u,k,(1<<3)^G5^(1<<j),insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    // Construct a0-(b1xNabla_T2) sequential source
    LatticePropagator
    MesA0B1xNablaT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<ForwardProp_t>& forward_headers,
				      const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0B1xNablaT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 s_{ijk}\gamma_j \nabla_k\f$  
      return project(-b1xVector_sym(threePtNablaVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0B1xNablaT2SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				     const multi1d<ForwardProp_t>& forward_headers,
				     const multi1d<LatticePropagator>& quark_propagators,
				     int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0B1xNablaT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 s_{ijk}\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int s = symTensor3d(getDerivDir(),j,k);
	  if (s != 0)
	    corr_fn += Real(s) * gammaSgn(1<<3,G5) * gammaSgn((1<<3)^G5,1<<j) * 
	      twoPtNabla(tmp,u,k,(1<<3)^G5^(1<<j),insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    // Construct a0-(b1xNabla_E) sequential source
    LatticePropagator
    MesA0B1xNablaESeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				     const multi1d<ForwardProp_t>& forward_headers,
				     const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0B1xNablaESeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 Q_{\alpha jk}\gamma_j \nabla_k\f$  
      return project(-b1xVector_E(threePtNablaVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0B1xNablaESeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				     const multi1d<ForwardProp_t>& forward_headers,
				     const multi1d<LatticePropagator>& quark_propagators,
				     int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0B1xNablaESeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 Q_{\alpha jk}\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  Real e = ETensor3d(getDerivDir(),j,k);
	  if (toBool(e != 0.0))
	    corr_fn += e * gammaSgn(1<<3,G5) * gammaSgn((1<<3)^G5,1<<j) * 
	      twoPtNabla(tmp,u,k,(1<<3)^G5^(1<<j),insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //---------------------------------------------------------------------------------
    // Construct a0-(a0xD_T2) sequential source
    LatticePropagator
    MesA0A0xDT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0A0xDT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv D_i\f$  
      return project(a0xVector(threePtD(quark_propagators[0],u,getDerivDir())));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0A0xDT2SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<ForwardProp_t>& forward_headers,
				 const multi1d<LatticePropagator>& quark_propagators,
				 int insertion)
    {
      START_CODE();

      check1Args("MesA0A0xDT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv D_i\f$  
      LatticeComplex corr_fn = twoPtD(tmp,u,getDerivDir(),0,insertion);
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    // Construct a0-(pionxD_T2) sequential source
    LatticePropagator
    MesA0PionxDT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				    const multi1d<ForwardProp_t>& forward_headers,
				    const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0PionxDT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_5 D_i\f$  
      return project(pionxVector(threePtD(quark_propagators[0],u,getDerivDir())));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0PionxDT2SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators,
				  int insertion)
    {
      START_CODE();

      int G5 = Ns*Ns-1;

      check1Args("MesA0PionxDT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5 D_i\f$  
      LatticeComplex corr_fn = twoPtD(tmp,u,getDerivDir(),G5,insertion);
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    // Construct a0-(pion_2xD_T2) sequential source
    LatticePropagator
    MesA0Pion2xDT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				     const multi1d<ForwardProp_t>& forward_headers,
				     const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0Pion2xDT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 D_i\f$  
      return project(pion_2xVector(threePtD(quark_propagators[0],u,getDerivDir())));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0Pion2xDT2SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				    const multi1d<ForwardProp_t>& forward_headers,
				    const multi1d<LatticePropagator>& quark_propagators,
				    int insertion)
    {
      START_CODE();

      int G5 = Ns*Ns-1;

      check1Args("MesA0Pion2xDT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 D_i\f$  
      LatticeComplex corr_fn = gammaSgn(1<<3,G5) * twoPtD(tmp,u,getDerivDir(),(1<<3)^G5,insertion);
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    // Construct a0-(a0_2xD_T2) sequential source
    LatticePropagator
    MesA0A02xDT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0A02xDT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_4 D_i\f$  
      return project(-a0_2xVector(threePtD(quark_propagators[0],u,getDerivDir())));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0A02xDT2SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators,
				  int insertion)
    {
      START_CODE();

      check1Args("MesA0A02xDT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_4 D_i\f$  
      LatticeComplex corr_fn = twoPtD(tmp,u,getDerivDir(),1<<3,insertion);
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //---------------------------------------------------------------------------------
    // Construct a0-(rhoxD_A2) sequential source
    LatticePropagator
    MesA0RhoxDA2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0RhoxDA2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_i D_i\f$  
      return project(-rhoxVector_sum(threePtDVector(quark_propagators[0],u)));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0RhoxDA2SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators,
				  int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0RhoxDA2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_i D_i\f$  
      for(int k=0; k < 3; ++k)
	corr_fn += twoPtD(tmp,u,k, 1<<k, insertion);
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //! Construct a0-(rhoxD_T1) sequential source
    LatticePropagator
    MesA0RhoxDT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0RhoxDT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv s_{ijk}\gamma_j D_k\f$  
      return project(-rhoxVector_sym(threePtDVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0RhoxDT1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators,
				  int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0RhoxDT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv s_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int s = symTensor3d(getDerivDir(),j,k);
	  if (s != 0)
	    corr_fn += Real(s) * twoPtD(tmp,u,k, 1<<j, insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //! Construct a0-(rhoxD_T2) sequential source
    LatticePropagator
    MesA0RhoxDT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0RhoxDT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j D_k\f$  
      return project(-rhoxVector_antisym(threePtDVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0RhoxDT2SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators,
				  int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0RhoxDT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int a = antiSymTensor3d(getDerivDir(),j,k);
	  if (a != 0)
	    corr_fn += Real(a) * twoPtD(tmp,u,k, 1<<j, insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //! Construct a0-(rhoxD_E) sequential source
    LatticePropagator
    MesA0RhoxDESeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0RhoxDESeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv Q_{\alpha jk}\gamma_j D_k\f$  
      return project(-rhoxVector_E(threePtDVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0RhoxDESeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<ForwardProp_t>& forward_headers,
				 const multi1d<LatticePropagator>& quark_propagators,
				 int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0RhoxDESeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv Q_{\alpha jk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  Real e = ETensor3d(getDerivDir(),j,k);
	  if (toBool(e != 0.0))
	    corr_fn += e * twoPtD(tmp,u,k, 1<<j, insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //---------------------------------------------------------------------------------
    // Construct a0-(rho_2xD_A2) sequential source
    LatticePropagator
    MesA0Rho2xDA2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				    const multi1d<ForwardProp_t>& forward_headers,
				    const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0Rho2xDA2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_4\gamma_i D_i\f$  
      return project(-rho_2xVector_sum(threePtDVector(quark_propagators[0],u)));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0Rho2xDA2SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators,
				  int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0Rho2xDA2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_4 \gamma_i D_i\f$  
      for(int k=0; k < 3; ++k)
	corr_fn += gammaSgn(1<<3,1<<k) * twoPtD(tmp,u,k, (1<<3)^(1<<k), insertion);
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //! Construct a0-(rho_2xD_T1) sequential source
    LatticePropagator
    MesA0Rho2xDT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				    const multi1d<ForwardProp_t>& forward_headers,
				    const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0Rho2xDT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv s_{ijk}\gamma_4\gamma_j D_k\f$  
      return project(-rho_2xVector_sym(threePtDVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0Rho2xDT1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators,
				  int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0Rho2xDT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv s_{ijk}\gamma_4\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int s = symTensor3d(getDerivDir(),j,k);
	  if (s != 0)
	    corr_fn += Real(s) * gammaSgn(1<<3,1<<j) * twoPtD(tmp,u,k, (1<<3)^(1<<j), insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //! Construct a0-(rho_2xD_T2) sequential source
    LatticePropagator
    MesA0Rho2xDT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0Rho2xDT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_4\gamma_j D_k\f$  
      return project(-rho_2xVector_antisym(threePtDVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0Rho2xDT2SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& quark_propagators,
				   int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0Rho2xDT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_4\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int a = antiSymTensor3d(getDerivDir(),j,k);
	  if (a != 0)
	    corr_fn += Real(a) * gammaSgn(1<<3,1<<j) * twoPtD(tmp,u,k, (1<<3)^(1<<j), insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //! Construct a0-(rho_2xD_E) sequential source
    LatticePropagator
    MesA0Rho2xDESeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0Rho2xDESeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv Q_{\alpha jk}\gamma_4\gamma_j D_k\f$  
      return project(-rho_2xVector_E(threePtDVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0Rho2xDESeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators,
				  int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0Rho2xDESeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv Q_{\alpha jk}\gamma_4\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  Real e = ETensor3d(getDerivDir(),j,k);
	  if (toBool(e != 0.0))
	    corr_fn += e * gammaSgn(1<<3,1<<j) * twoPtD(tmp,u,k, (1<<3)^(1<<j), insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //---------------------------------------------------------------------------------
    // Construct a0-(a1xD_A2) sequential source
    LatticePropagator
    MesA0A1xDA2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0A1xDA2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_5\gamma_i D_i\f$  
      return project(a1xVector_sum(threePtDVector(quark_propagators[0],u)));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0A1xDA2SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<ForwardProp_t>& forward_headers,
				 const multi1d<LatticePropagator>& quark_propagators,
				 int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0A1xDA2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5\gamma_i D_i\f$  
      for(int k=0; k < 3; ++k)
	corr_fn += gammaSgn(G5,1<<k) * twoPtD(tmp,u,k,G5^(1<<k),insertion);
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    // Construct a0-(a1xD_T1) sequential source
    LatticePropagator
    MesA0A1xDT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0A1xDT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j D_k\f$  
      return project(a1xVector_sym(threePtDVector(quark_propagators[0],u), getDerivDir()));
    }

    
    // Compute the 2-pt at the sink
    Complex 
    MesA0A1xDT1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<ForwardProp_t>& forward_headers,
				 const multi1d<LatticePropagator>& quark_propagators,
				 int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0A1xDT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int s = symTensor3d(getDerivDir(),j,k);
	  if (s != 0)
	    corr_fn += Real(s) * gammaSgn(G5,1<<j) * twoPtD(tmp,u,k,G5^(1<<j),insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    // Construct a0-(a1xD_T2) sequential source
    LatticePropagator
    MesA0A1xDT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0A1xDT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_5 \epsilon_{ijk}\gamma_j D_k\f$  
      return project(a1xVector_antisym(threePtDVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0A1xDT2SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<ForwardProp_t>& forward_headers,
				 const multi1d<LatticePropagator>& quark_propagators,
				 int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0A1xDT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5 \epsilon_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int a = antiSymTensor3d(getDerivDir(),j,k);
	  if (a != 0)
	    corr_fn += Real(a) * gammaSgn(G5,1<<j) * twoPtD(tmp,u,k,G5^(1<<j),insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    // Construct a0-(a1xD_E) sequential source
    LatticePropagator
    MesA0A1xDESeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<ForwardProp_t>& forward_headers,
				 const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0A1xDESeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_5 Q_{\alpha jk}\gamma_j D_k\f$  
      return project(a1xVector_E(threePtDVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0A1xDESeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				const multi1d<ForwardProp_t>& forward_headers,
				const multi1d<LatticePropagator>& quark_propagators,
				int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0A1xDESeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5 Q_{\alpha jk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  Real e = ETensor3d(getDerivDir(),j,k);
	  if (toBool(e != 0.0))
	    corr_fn += e * gammaSgn(G5,1<<j) * twoPtD(tmp,u,k,G5^(1<<j),insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //---------------------------------------------------------------------------------
    // Construct a0-(b1xD_A2) sequential source
    LatticePropagator
    MesA0B1xDA2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0B1xDA2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 \gamma_i D_i\f$  
      return project(-b1xVector_sum(threePtDVector(quark_propagators[0],u)));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0B1xDA2SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<ForwardProp_t>& forward_headers,
				 const multi1d<LatticePropagator>& quark_propagators,
				 int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0B1xDA2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 \gamma_i D_i\f$  
      for(int k=0; k < 3; ++k)
	corr_fn += gammaSgn(1<<3,G5) * gammaSgn((1<<3)^G5,1<<k) * 
	  twoPtD(tmp,u,k,(1<<3)^G5^(1<<k),insertion);

      END_CODE();

      return momentumProject(corr_fn);
    }


    // Construct a0-(b1xD_T1) sequential source
    LatticePropagator
    MesA0B1xDT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0B1xDT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 s_{ijk}\gamma_j D_k\f$  
      return project(-b1xVector_sym(threePtDVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0B1xDT1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<ForwardProp_t>& forward_headers,
				 const multi1d<LatticePropagator>& quark_propagators,
				 int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0B1xDT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 s_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int s = symTensor3d(getDerivDir(),j,k);
	  if (s != 0)
	    corr_fn += Real(s) * gammaSgn(1<<3,G5) * gammaSgn((1<<3)^G5,1<<j) * 
	      twoPtD(tmp,u,k,(1<<3)^G5^(1<<j),insertion);
	}

      END_CODE();

      return momentumProject(corr_fn);
    }


    // Construct a0-(b1xD_T2) sequential source
    LatticePropagator
    MesA0B1xDT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0B1xDT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 s_{ijk}\gamma_j D_k\f$  
      return project(-b1xVector_antisym(threePtDVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0B1xDT2SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<ForwardProp_t>& forward_headers,
				 const multi1d<LatticePropagator>& quark_propagators,
				 int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0B1xDT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 \epsilon_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int a = antiSymTensor3d(getDerivDir(),j,k);
	  if (a != 0)
	    corr_fn += Real(a) * gammaSgn(1<<3,G5) * gammaSgn((1<<3)^G5,1<<j) * 
	      twoPtD(tmp,u,k,(1<<3)^G5^(1<<j),insertion);
	}

      END_CODE();

      return momentumProject(corr_fn);
    }


    // Construct a0-(b1xD_E) sequential source
    LatticePropagator
    MesA0B1xDESeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<ForwardProp_t>& forward_headers,
				 const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0B1xDESeqSrc", quark_propagators);
      setTSrce(forward_headers);
 
      // \f$\Gamma_f \equiv \gamma_4\gamma_5 Q_{\alpha jk}\gamma_j D_k\f$  
      return project(-b1xVector_E(threePtDVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0B1xDESeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				const multi1d<ForwardProp_t>& forward_headers,
				const multi1d<LatticePropagator>& quark_propagators,
				int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0B1xDESeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];
 
      // \f$\Gamma_f \equiv \gamma_4\gamma_5 Q_{\alpha jk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  Real e = ETensor3d(getDerivDir(),j,k);
	  if (toBool(e != 0.0))
	    corr_fn += e * gammaSgn(1<<3,G5) * gammaSgn((1<<3)^G5,1<<j) * 
	      twoPtD(tmp,u,k,(1<<3)^G5^(1<<j),insertion);
	}

      END_CODE();

      return momentumProject(corr_fn);
    }


    //---------------------------------------------------------------------------------
    //! Construct a0-(a0xB_T1) sequential source
    LatticePropagator
    MesA0A0xBT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0A0xBT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv B_i\f$  
      return project(a0xVector(threePtB(quark_propagators[0],u,getDerivDir())));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0A0xBT1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<ForwardProp_t>& forward_headers,
				 const multi1d<LatticePropagator>& quark_propagators,
				 int insertion)
    {
      START_CODE();

      int G5 = Ns*Ns-1;

      check1Args("MesA0A0xBT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv B_i\f$  
      LatticeComplex corr_fn = twoPtB(tmp,u,getDerivDir(),0,insertion);
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //! Construct a0-(pionxB_T1) sequential source
    LatticePropagator
    MesA0PionxBT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				    const multi1d<ForwardProp_t>& forward_headers,
				    const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0PionxBT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_5 B_i\f$  
      return project(pionxVector(threePtB(quark_propagators[0],u,getDerivDir())));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0PionxBT1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& quark_propagators,
				   int insertion)
    {
      START_CODE();

      int G5 = Ns*Ns-1;

      check1Args("MesA0PionxBT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5 B_i\f$  
      LatticeComplex corr_fn = twoPtB(tmp,u,getDerivDir(), G5, insertion);
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //! Construct a0-(pion_2xB_T1) sequential source
    LatticePropagator
    MesA0Pion2xBT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				     const multi1d<ForwardProp_t>& forward_headers,
				     const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0Pion2xBT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 B_i\f$  
      return project(pion_2xVector(threePtB(quark_propagators[0],u,getDerivDir())));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0Pion2xBT1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				    const multi1d<ForwardProp_t>& forward_headers,
				    const multi1d<LatticePropagator>& quark_propagators,
				    int insertion)
    {
      START_CODE();

      int G5 = Ns*Ns-1;

      check1Args("MesA0Pion2xBT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 B_i\f$  
      LatticeComplex corr_fn = gammaSgn(1<<3,G5) * twoPtB(tmp,u,getDerivDir(), (1<<3)^G5, insertion);
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //! Construct a0-(a0_2xB_T1) sequential source
    LatticePropagator
    MesA0A02xBT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0A02xBT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_4 B_i\f$  
      return project(-a0_2xVector(threePtB(quark_propagators[0],u,getDerivDir())));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0A02xBT1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators,
				  int insertion)
    {
      START_CODE();

      check1Args("MesA0A02xBT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_4 B_i\f$  
      LatticeComplex corr_fn = twoPtB(tmp,u,getDerivDir(), 1<<3, insertion);
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //---------------------------------------------------------------------------------
    //! Construct a0-(rhoxB_A1) sequential source
    LatticePropagator
    MesA0RhoxBA1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0RhoxBA1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_i B_i\f$  
      return project(-rhoxVector_sum(threePtBVector(quark_propagators[0],u)));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0RhoxBA1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators,
				  int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0RhoxBA1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_i B_i\f$  
      for(int k=0; k < 3; ++k)
	corr_fn += twoPtB(tmp,u,k, 1<<k, insertion);
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //! Construct a0-(rhoxB_T1) sequential source
    LatticePropagator
    MesA0RhoxBT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0RhoxBT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j B_k\f$  
      return project(-rhoxVector_antisym(threePtBVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0RhoxBT1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators,
				  int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0RhoxBT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int a = antiSymTensor3d(getDerivDir(),j,k);
	  if (a != 0)
	    corr_fn += Real(a) * twoPtB(tmp,u,k, 1<<j, insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //! Construct a0-(rhoxB_T2) sequential source
    LatticePropagator
    MesA0RhoxBT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0RhoxBT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv s_{ijk}\gamma_j B_k\f$  
      return project(-rhoxVector_sym(threePtBVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0RhoxBT2SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators,
				  int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0RhoxBT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv s_{ijk}\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int s = symTensor3d(getDerivDir(),j,k);
	  if (s != 0)
	    corr_fn += Real(s) * twoPtB(tmp,u,k, 1<<j, insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //! Construct a0-(rhoxB_E) sequential source
    LatticePropagator
    MesA0RhoxBESeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0RhoxBESeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv Q_{\alpha jk}\gamma_j B_k\f$  
      return project(-rhoxVector_E(threePtBVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0RhoxBESeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators,
				  int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0RhoxBESeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv Q_{\alpha jk}\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  Real e = ETensor3d(getDerivDir(),j,k);
	  if (toBool(e != 0.0))
	    corr_fn += e * twoPtB(tmp,u,k, 1<<j, insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }



    //---------------------------------------------------------------------------------
    //! Construct a0-(rho_2xB_A1) sequential source
    LatticePropagator
    MesA0Rho2xBA1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0Rho2xBA1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_4\gamma_i B_i\f$  
      return project(-rho_2xVector_sum(threePtBVector(quark_propagators[0],u)));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0Rho2xBA1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators,
				  int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0Rho2xBA1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_4\gamma_i B_i\f$  
      for(int k=0; k < 3; ++k)
	corr_fn += gammaSgn(1<<3,1<<k) * twoPtB(tmp,u,k, (1<<3)^(1<<k), insertion);
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //! Construct a0-(rho_2xB_T1) sequential source
    LatticePropagator
    MesA0Rho2xBT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0Rho2xBT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_4\gamma_j B_k\f$  
      return project(-rho_2xVector_antisym(threePtBVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0Rho2xBT1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& quark_propagators,
				   int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0Rho2xBT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_4\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int a = antiSymTensor3d(getDerivDir(),j,k);
	  if (a != 0)
	    corr_fn += Real(a) * gammaSgn(1<<3,1<<j) * twoPtB(tmp,u,k, (1<<3)^(1<<j), insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //! Construct a0-(rho_2xB_T2) sequential source
    LatticePropagator
    MesA0Rho2xBT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0Rho2xBT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv s_{ijk}\gamma_4\gamma_j B_k\f$  
      return project(-rho_2xVector_sym(threePtBVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0Rho2xBT2SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& quark_propagators,
				   int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0Rho2xBT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv s_{ijk}\gamma_4\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int s = symTensor3d(getDerivDir(),j,k);
	  if (s != 0)
	    corr_fn += Real(s) * gammaSgn(1<<3,1<<j) * twoPtB(tmp,u,k, (1<<3)^(1<<j), insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //! Construct a0-(rho_2xB_E) sequential source
    LatticePropagator
    MesA0Rho2xBESeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0Rho2xBESeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv Q_{\alpha jk}\gamma_4\gamma_j B_k\f$  
      return project(-rho_2xVector_E(threePtBVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0Rho2xBESeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators,
				  int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0Rho2xBESeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv Q_{\alpha jk}\gamma_4\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  Real e = ETensor3d(getDerivDir(),j,k);
	  if (toBool(e != 0.0))
	    corr_fn += e * gammaSgn(1<<3,1<<j) * twoPtB(tmp,u,k, (1<<3)^(1<<j), insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //---------------------------------------------------------------------------------
    //! Construct a0-(a1xB_A1) sequential source
    LatticePropagator
    MesA0A1xBA1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0A1xBA1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_5 \gamma_i B_i\f$  
      return project(a1xVector_sum(threePtBVector(quark_propagators[0],u)));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0A1xBA1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<ForwardProp_t>& forward_headers,
				 const multi1d<LatticePropagator>& quark_propagators,
				 int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0A1xBA1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5 \gamma_i B_i\f$  
      for(int k=0; k < 3; ++k)
	corr_fn += gammaSgn(G5,1<<k) * twoPtB(tmp,u,k, G5^(1<<k), insertion);
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //! Construct a0-(a1xB_T1) sequential source
    LatticePropagator
    MesA0A1xBT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0A1xBT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_5 \epsilon_{ijk}\gamma_j B_k\f$  
      return project(a1xVector_antisym(threePtBVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0A1xBT1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<ForwardProp_t>& forward_headers,
				 const multi1d<LatticePropagator>& quark_propagators,
				 int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0A1xBT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5 \epsilon_{ijk}\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int a = antiSymTensor3d(getDerivDir(),j,k);
	  if (a != 0)
	    corr_fn += Real(a) * gammaSgn(G5,1<<j) * twoPtB(tmp,u,k, G5^(1<<j), insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //! Construct a0-(a1xB_T2) sequential source
    LatticePropagator
    MesA0A1xBT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0A1xBT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j B_k\f$  
      return project(a1xVector_sym(threePtBVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0A1xBT2SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<ForwardProp_t>& forward_headers,
				 const multi1d<LatticePropagator>& quark_propagators,
				 int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0A1xBT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int s = symTensor3d(getDerivDir(),j,k);
	  if (s != 0)
	    corr_fn += Real(s) * gammaSgn(G5,1<<j) * twoPtB(tmp,u,k, G5^(1<<j), insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //! Construct a0-(a1xB_E) sequential source
    LatticePropagator
    MesA0A1xBESeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0A1xBESeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_5 Q_{\alpha jk}\gamma_j B_k\f$  
      return project(a1xVector_E(threePtBVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0A1xBESeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				const multi1d<ForwardProp_t>& forward_headers,
				const multi1d<LatticePropagator>& quark_propagators,
				int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0A1xBESeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5 Q_{\alpha jk}\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  Real e = ETensor3d(getDerivDir(),j,k);
	  if (toBool(e != 0.0))
	    corr_fn += e * gammaSgn(G5,1<<j) * twoPtB(tmp,u,k, G5^(1<<j), insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //---------------------------------------------------------------------------------
    //! Construct a0-(b1xB_A1) sequential source
    LatticePropagator
    MesA0B1xBA1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0B1xBA1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_4 \gamma_5 \gamma_i B_i\f$  
      return project(-b1xVector_sum(threePtBVector(quark_propagators[0],u)));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0B1xBA1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<ForwardProp_t>& forward_headers,
				 const multi1d<LatticePropagator>& quark_propagators,
				 int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0B1xBA1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_4 \gamma_5 \gamma_i B_i\f$  
      for(int k=0; k < 3; ++k)
	corr_fn += gammaSgn(1<<3,G5) * gammaSgn((1<<3)^G5,1<<k) * 
	  twoPtB(tmp,u,k, (1<<3)^G5^(1<<k), insertion);
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //! Construct a0-(b1xB_T1) sequential source
    LatticePropagator
    MesA0B1xBT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0B1xBT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 \epsilon_{ijk}\gamma_j B_k\f$  
      return project(-b1xVector_antisym(threePtBVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0B1xBT1SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<ForwardProp_t>& forward_headers,
				 const multi1d<LatticePropagator>& quark_propagators,
				 int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0B1xBT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_4 \gamma_5 \epsilon_{ijk}\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int a = antiSymTensor3d(getDerivDir(),j,k);
	  if (a != 0)
	    corr_fn += Real(a) * gammaSgn(1<<3,G5) * gammaSgn((1<<3)^G5,1<<j) * 
	      twoPtB(tmp,u,k, (1<<3)^G5^(1<<j), insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //! Construct a0-(b1xB_T2) sequential source
    LatticePropagator
    MesA0B1xBT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0B1xBT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_4 \gamma_5 s_{ijk}\gamma_j B_k\f$  
      return project(-b1xVector_sym(threePtBVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0B1xBT2SeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<ForwardProp_t>& forward_headers,
				 const multi1d<LatticePropagator>& quark_propagators,
				 int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0B1xBT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_4 \gamma_5 s_{ijk}\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int s = symTensor3d(getDerivDir(),j,k);
	  if (s != 0)
	    corr_fn += Real(s) * gammaSgn(1<<3,G5) * gammaSgn((1<<3)^G5,1<<j) * 
	      twoPtB(tmp,u,k, (1<<3)^G5^(1<<j), insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //! Construct a0-(b1xB_E) sequential source
    LatticePropagator
    MesA0B1xBESeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<ForwardProp_t>& forward_headers,
				 const multi1d<LatticePropagator>& quark_propagators)
    {
      check1Args("MesA0B1xBESeqSrc", quark_propagators);
      setTSrce(forward_headers);

      // \f$\Gamma_f \equiv \gamma_4 \gamma_5 Q_{\alpha jk}\gamma_j B_k\f$  
      return project(-b1xVector_E(threePtBVector(quark_propagators[0],u), getDerivDir()));
    }


    // Compute the 2-pt at the sink
    Complex 
    MesA0B1xBESeqSrc::twoPtSink(const multi1d<LatticeColorMatrix>& u,
				const multi1d<ForwardProp_t>& forward_headers,
				const multi1d<LatticePropagator>& quark_propagators,
				int insertion)
    {
      START_CODE();

      LatticeComplex corr_fn = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0B1xBESeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_4 \gamma_5 Q_{\alpha jk}\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  Real e = ETensor3d(getDerivDir(),j,k);
	  if (toBool(e != 0.0))
	    corr_fn += e * gammaSgn(1<<3,G5) * gammaSgn((1<<3)^G5,1<<j) * 
	      twoPtB(tmp,u,k, (1<<3)^G5^(1<<j), insertion);
	}
      
      END_CODE();

      return momentumProject(corr_fn);
    }


    //---------------------------------------------------------------------------------
    namespace
    {
      //! Local registration flag
      bool registered = false;
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	//! Register all the factories
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a0xNABLA_T1"),
										      mesA0A0xNablaT1SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-pionxNABLA_T1"),
										      mesA0PionxNablaT1SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-pion_2xNABLA_T1"),
										      mesA0Pion2xNablaT1SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a0_2xNABLA_T1"),
										      mesA0A02xNablaT1SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rhoxNABLA_A1"),
										      mesA0RhoxNablaA1SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rhoxNABLA_T1"),
										      mesA0RhoxNablaT1SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rhoxNABLA_T2"),
										      mesA0RhoxNablaT2SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rhoxNABLA_E"),
										      mesA0RhoxNablaESeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rho_2xNABLA_A1"),
										      mesA0Rho2xNablaA1SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rho_2xNABLA_T1"),
										      mesA0Rho2xNablaT1SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rho_2xNABLA_T2"),
										      mesA0Rho2xNablaT2SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rho_2xNABLA_E"),
										      mesA0Rho2xNablaESeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a1xNABLA_A1"),
										      mesA0A1xNablaA1SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a1xNABLA_T1"),
										      mesA0A1xNablaT1SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a1xNABLA_T2"),
										      mesA0A1xNablaT2SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a1xNABLA_E"),
										      mesA0A1xNablaESeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-b1xNABLA_A1"),
										      mesA0B1xNablaA1SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-b1xNABLA_T1"),
										      mesA0B1xNablaT1SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-b1xNABLA_T2"),
										      mesA0B1xNablaT2SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-b1xNABLA_E"),
										      mesA0B1xNablaESeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a0xD_T2"),
										      mesA0A0xDT2SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-pionxD_T2"),
										      mesA0PionxDT2SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-pion_2xD_T2"),
										      mesA0Pion2xDT2SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a0_2xD_T2"),
										      mesA0A02xDT2SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rhoxD_A2"),
										      mesA0RhoxDA2SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rhoxD_T1"),
										      mesA0RhoxDT1SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rhoxD_T2"),
										      mesA0RhoxDT2SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rhoxD_E"),
										      mesA0RhoxDESeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rho_2xD_A2"),
										      mesA0Rho2xDA2SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rho_2xD_T1"),
										      mesA0Rho2xDT1SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rho_2xD_T2"),
										      mesA0Rho2xDT2SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rho_2xD_E"),
										      mesA0Rho2xDESeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a1xD_A2"),
										      mesA0A1xDA2SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a1xD_T1"),
										      mesA0A1xDT1SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a1xD_T2"),
										      mesA0A1xDT2SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a1xD_E"),
										      mesA0A1xDESeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-b1xD_A2"),
										      mesA0B1xDA2SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-b1xD_T1"),
										      mesA0B1xDT1SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-b1xD_T2"),
										      mesA0B1xDT2SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-b1xD_E"),
										      mesA0B1xDESeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a0xB_T1"),
										      mesA0A0xBT1SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-pionxB_T1"),
										      mesA0PionxBT1SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-pion_2xB_T1"),
										      mesA0Pion2xBT1SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a0_2xB_T1"),
										      mesA0A02xBT1SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rhoxB_A1"),
										      mesA0RhoxBA1SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rhoxB_T1"),
										      mesA0RhoxBT1SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rhoxB_T2"),
										      mesA0RhoxBT2SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rhoxB_E"),
										      mesA0RhoxBESeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rho_2xB_A1"),
										      mesA0Rho2xBA1SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rho_2xB_T1"),
										      mesA0Rho2xBT1SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rho_2xB_T2"),
										      mesA0Rho2xBT2SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rho_2xB_E"),
										      mesA0Rho2xBESeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a1xB_A1"),
										      mesA0A1xBA1SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a1xB_T1"),
										      mesA0A1xBT1SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a1xB_T2"),
										      mesA0A1xBT2SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a1xB_E"),
										      mesA0A1xBESeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-b1xB_A1"),
										      mesA0B1xBA1SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-b1xB_T1"),
										      mesA0B1xBT1SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-b1xB_T2"),
										      mesA0B1xBT2SeqSrc);
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-b1xB_E"),
										      mesA0B1xBESeqSrc);

	registered = true;
      }
      return success;
    }

  }  // end namespace DerivMesonSeqSourceEnv

}  // end namespace Chroma
