// $Id: derivmesonseqsrc_w.cc,v 2.4 2006-02-10 02:51:57 edwards Exp $
/*! \file
 *  \brief Construct meson sequential sources.
 */

#include "meas/hadron/derivmesonseqsrc_w.h"
#include "meas/hadron/seqsource_factory_w.h"

#include "util/ferm/symtensor.h"
#include "util/ferm/antisymtensor.h"

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


  //! Meson sequential sources
  /*! \ingroup hadron */
  namespace DerivMesonSeqSourceEnv
  { 
    //! Anonymous namespace
    namespace
    {

      //! Private displacement
      LatticePropagator displacement(const multi1d<LatticeColorMatrix>& u, 
				     const LatticePropagator& psi, 
				     int length, int dir)
      {
	LatticePropagator chi = psi;

	if (length > 0)
	{
	  for(int n = 0; n < length; ++n)
	  {
	    LatticePropagator tmp = shift(chi, FORWARD, dir);
	    chi = u[dir] * tmp;
	  }
	}
	else // If length = or < 0.  If length == 0, does nothing.
	{
	  for(int n = 0; n > length; --n)
	  {
	    LatticePropagator tmp = shift(adj(u[dir])*chi, BACKWARD, dir);
	    chi = tmp;
	  }
	}
	return chi;
      }


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
				  int mu, int deriv_length)
      {
//	return LatticePropagator(shift(F*u[mu], BACKWARD, mu) - shift(F, FORWARD, mu)*adj(u[mu]));
	return displacement(u, F, -deriv_length, mu) - displacement(u, F, deriv_length, mu);
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
			      int mu, int deriv_length)
      {
	LatticePropagator tmp = zero;

	// Slow implementation - to speed up could compute once the \nabla_j deriv
	for(int j=0; j < 3; ++j)
	  for(int k=0; k < 3; ++k)
	  {
	    if (symTensor3d(mu,j,k) != 0)
	      tmp += leftNabla(leftNabla(F,u,j,deriv_length), u,k,deriv_length);
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
			      int mu, int deriv_length)
      {
	LatticePropagator tmp = zero;

	// Slow implementation - to speed up could compute once the \nabla_j deriv
	for(int j=0; j < 3; ++j)
	  for(int k=0; k < 3; ++k)
	  {
	    if (antiSymTensor3d(mu,j,k) != 0)
	      tmp += Real(antiSymTensor3d(mu,j,k)) * leftNabla(leftNabla(F,u,j,deriv_length), u,k,deriv_length);
	  }

	return tmp;
      }



      //-------------------- callback functions ---------------------------------------

      //! Construct pion_1-(A2=A1xD_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesPionA1xDT2SeqSrc(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesPionA1xDT2SeqSrc(ParamsDir(xml_in, path));
      }

#if 0
      //! Construct pion_1-(PionxNabla_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesPionPionxNablaT1SeqSrc(XMLReader& xml_in,
								    const std::string& path)
      {
	return new MesPionPionxNablaT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(A0_2xNabla_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesPionA02xNablaT1SeqSrc(XMLReader& xml_in,
								   const std::string& path)
      {
	return new MesPionA02xNablaT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(RhoxNabla_A1) sequential source
      HadronSeqSource<LatticePropagator>* mesPionRhoxNablaA1SeqSrc(XMLReader& xml_in,
								   const std::string& path)
      {
	return new MesPionRhoxNablaA1SeqSrc(Params(xml_in, path));
      }

      //! Construct pion_1-(RhoxNabla_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesPionRhoxNablaT1SeqSrc(XMLReader& xml_in,
								   const std::string& path)
      {
	return new MesPionRhoxNablaT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(RhoxNabla_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesPionRhoxNablaT2SeqSrc(XMLReader& xml_in,
								   const std::string& path)
      {
	return new MesPionRhoxNablaT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(A1xNabla_A1) sequential source
      HadronSeqSource<LatticePropagator>* mesPionA1xNablaA1SeqSrc(XMLReader& xml_in,
								  const std::string& path)
      {
	return new MesPionA1xNablaA1SeqSrc(Params(xml_in, path));
      }

      //! Construct pion_1-(A1xNabla_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesPionA1xNablaT2SeqSrc(XMLReader& xml_in,
								  const std::string& path)
      {
	return new MesPionA1xNablaT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(A1xNabla_E) sequential source
      HadronSeqSource<LatticePropagator>* mesPionA1xNablaESeqSrc(XMLReader& xml_in,
								 const std::string& path)
      {
	return new MesPionA1xNablaESeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(B1xNabla_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesPionB1xNablaT1SeqSrc(XMLReader& xml_in,
								  const std::string& path)
      {
	return new MesPionB1xNablaT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(A0_2xD_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesPionA02xDT2SeqSrc(XMLReader& xml_in,
							       const std::string& path)
      {
	return new MesPionA02xDT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(A1xD_A2) sequential source
      HadronSeqSource<LatticePropagator>* mesPionA1xDA2SeqSrc(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesPionA1xDA2SeqSrc(Params(xml_in, path));
      }

      //! Construct pion_1-(A1xD_E) sequential source
      HadronSeqSource<LatticePropagator>* mesPionA1xDESeqSrc(XMLReader& xml_in,
							     const std::string& path)
      {
	return new MesPionA1xDESeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(A1xD_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesPionA1xDT1SeqSrc(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesPionA1xDT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(A1xD_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesPionA1xDT2SeqSrc(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesPionA1xDT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(B1xD_A2) sequential source
      HadronSeqSource<LatticePropagator>* mesPionB1xDA2SeqSrc(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesPionB1xDA2SeqSrc(Params(xml_in, path));
      }

      //! Construct pion_1-(B1xD_E) sequential source
      HadronSeqSource<LatticePropagator>* mesPionB1xDESeqSrc(XMLReader& xml_in,
							     const std::string& path)
      {
	return new MesPionB1xDESeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(B1xD_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesPionB1xDT1SeqSrc(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesPionB1xDT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(B1xD_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesPionB1xDT2SeqSrc(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesPionB1xDT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(RhoxD_A2) sequential source
      HadronSeqSource<LatticePropagator>* mesPionRhoxDA2SeqSrc(XMLReader& xml_in,
							       const std::string& path)
      {
	return new MesPionRhoxDA2SeqSrc(Params(xml_in, path));
      }

      //! Construct pion_1-(RhoxD_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesPionRhoxDT1SeqSrc(XMLReader& xml_in,
							       const std::string& path)
      {
	return new MesPionRhoxDT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(RhoxD_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesPionRhoxDT2SeqSrc(XMLReader& xml_in,
							       const std::string& path)
      {
	return new MesPionRhoxDT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(RhoxD_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesPionPionxDT2SeqSrc(XMLReader& xml_in,
								const std::string& path)
      {
	return new MesPionPionxDT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(RhoxB_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesPionRhoxBT1SeqSrc(XMLReader& xml_in,
								const std::string& path)
      {
	return new MesPionRhoxBT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(RhoxB_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesPionRhoxBT2SeqSrc(XMLReader& xml_in,
								const std::string& path)
      {
	return new MesPionRhoxBT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(A1xB_A1) sequential source
      HadronSeqSource<LatticePropagator>* mesPionA1xBA1SeqSrc(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesPionA1xBA1SeqSrc(Params(xml_in, path));
      }

      //! Construct pion_1-(A1xB_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesPionA1xBT1SeqSrc(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesPionA1xBT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(A1xB_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesPionA1xBT2SeqSrc(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesPionA1xBT2SeqSrc(ParamsDir(xml_in, path));
      }
#endif

    } // end anonymous namespace


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


    // Construct pion_1-(A2=A1xD_T2) sequential source
    // See corresponding .h file for doxygen comments
    LatticePropagator
    MesPionA1xDT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				    const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesPionA1xDT2SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // Slow implementation - to speed up could compute once the B_k deriv
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (antiSymTensor3d(params.deriv_dir,j,k) != 0)
	    tmp += Real(antiSymTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * leftD(tmp,u,k,params.deriv_length));
	}
      
      END_CODE();

      return project(fin, params.sink_mom, params.t_sink, params.j_decay);
    }


    // Register all the possible deriv mesons
    bool registerAll(void) 
    {
      bool success = true;

      //! Register all the factories
      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-A2=A1xD_T2"),
										    mesPionA1xDT2SeqSrc);

      return success;
    }

    const bool registered = registerAll();

  }  // end namespace DerivMesonSeqSourceEnv

}  // end namespace Chroma



  
