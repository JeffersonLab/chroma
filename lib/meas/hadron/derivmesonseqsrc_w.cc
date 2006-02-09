// $Id: derivmesonseqsrc_w.cc,v 2.2 2006-02-09 03:50:06 edwards Exp $
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


  //! Meson sequential sources
  /*! \ingroup hadron */
  namespace DerivMesonSeqSourceEnv
  { 
    //! Anonymous namespace
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



      //-------------------- callback functions ---------------------------------------

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
      HadronSeqSource<LatticePropagator>* mesPionA1xDT2SeqSrc(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesPionA1xDT2SeqSrc(ParamsDir(xml_in, path));
      }

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
	  if (antiSymTensor3d(params.deriv_dir,j,k) != 0)
	    tmp += Real(antiSymTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * leftD(tmp,u,k));
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



  
