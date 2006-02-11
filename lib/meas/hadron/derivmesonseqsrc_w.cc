// $Id: derivmesonseqsrc_w.cc,v 2.8 2006-02-11 21:23:17 edwards Exp $
/*! \file
 *  \brief Construct meson sequential sources.
 */

#include "meas/hadron/derivmesonseqsrc_w.h"
#include "meas/hadron/seqsource_factory_w.h"

#include "util/ferm/symtensor.h"
#include "util/ferm/antisymtensor.h"
#include "util/ferm/etensor.h"

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

      //! Construct pion_1-(PionxNabla_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesPionPionxNablaT1SeqSrc(XMLReader& xml_in,
								    const std::string& path)
      {
	return new MesPionPionxNablaT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(A0xNabla_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesPionA0xNablaT1SeqSrc(XMLReader& xml_in,
								  const std::string& path)
      {
	return new MesPionA0xNablaT1SeqSrc(ParamsDir(xml_in, path));
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

      //! Construct pion_1-(PionxD_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesPionPionxDT2SeqSrc(XMLReader& xml_in,
								const std::string& path)
      {
	return new MesPionPionxDT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(PionxB_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesPionPionxBT1SeqSrc(XMLReader& xml_in,
								const std::string& path)
      {
	return new MesPionPionxBT1SeqSrc(ParamsDir(xml_in, path));
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


    // Construct pion_1-(PionxNabla_T1) sequential source
    // See corresponding .h file for doxygen comments
    LatticePropagator
    MesPionPionxNablaT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
					  const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin;
      const int G5 = Ns*Ns-1;

      check1Args("MesPionPionxNablaT1SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // \f$\Gamma_f \equiv \nabla_i\f$
      fin = Gamma(G5) * leftNabla(tmp,u,params.deriv_dir,params.deriv_length);
      
      END_CODE();

      return project(fin, params.sink_mom, params.t_sink, params.j_decay);
    }


    // Construct pion_1-(A0xNabla_T1) sequential source
    // See corresponding .h file for doxygen comments
    LatticePropagator
    MesPionA0xNablaT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
					const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin;

      check1Args("MesPionA0xNablaT1SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      fin = leftNabla(tmp,u,params.deriv_dir,params.deriv_length);
      
      END_CODE();

      return project(fin, params.sink_mom, params.t_sink, params.j_decay);
    }


    // Construct pion_1-(A0_2xNabla_T1) sequential source
    // See corresponding .h file for doxygen comments
    LatticePropagator
    MesPionA02xNablaT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
					 const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin;

      check1Args("MesPionA02xNablaT1SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);
      
      // \f$\Gamma_f \equiv \gamma_4 \nabla_i\f$
      fin = Gamma(1 << 3) * leftNabla(tmp,u,params.deriv_dir,params.deriv_length);
      
      END_CODE();

      return project(fin, params.sink_mom, params.t_sink, params.j_decay);
    }


    // Construct pion_1-(RhoxNabla_A1) sequential source
    LatticePropagator
    MesPionRhoxNablaA1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
					 const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin = zero;

      check1Args("MesPionRhoxNablaA1SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // \f$\Gamma_f \equiv \gamma_i\nabla_i\f$
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * leftNabla(tmp,u,k,params.deriv_length);
      
      END_CODE();

      return project(fin, params.sink_mom, params.t_sink, params.j_decay);
    }


    // Construct pion_1-(RhoxNabla_T1) sequential source
    LatticePropagator
    MesPionRhoxNablaT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
					 const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin = zero;

      check1Args("MesPionRhoxNablaT1SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (antiSymTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(antiSymTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * leftNabla(tmp,u,k,params.deriv_length));
	}
      
      END_CODE();

      return project(fin, params.sink_mom, params.t_sink, params.j_decay);
    }


    // Construct pion_1-(RhoxNabla_T2) sequential source
    LatticePropagator
    MesPionRhoxNablaT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
					 const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin = zero;

      check1Args("MesPionRhoxNablaT2SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // \f$\Gamma_f \equiv s_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (symTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(symTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * leftNabla(tmp,u,k,params.deriv_length));
	}
      
      END_CODE();

      return project(fin, params.sink_mom, params.t_sink, params.j_decay);
    }


    // Construct pion_1-(A1xNabla_A1) sequential source
    LatticePropagator
    MesPionA1xNablaA1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
					const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesPionA1xNablaA1SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // \f$\Gamma_f \equiv \gamma_5\gamma_i \nabla_i\f$  
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * leftNabla(tmp,u,k,params.deriv_length);
      
      END_CODE();

      return project(LatticePropagator(Gamma(G5) * fin), 
		     params.sink_mom, params.t_sink, params.j_decay);
    }


    // Construct pion_1-(A1xNabla_T2) sequential source
    LatticePropagator
    MesPionA1xNablaT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
					const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesPionA1xNablaT2SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (symTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(symTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * leftNabla(tmp,u,k,params.deriv_length));
	}
      
      END_CODE();

      return project(LatticePropagator(Gamma(G5) * fin), 
		     params.sink_mom, params.t_sink, params.j_decay);
    }


    // Construct pion_1-(A1xNabla_E) sequential source
    LatticePropagator
    MesPionA1xNablaESeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				       const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesPionA1xNablaESeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // \f$\Gamma_f \equiv \gamma_5 S_{\alpha jk}\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (ETensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(ETensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * leftNabla(tmp,u,k,params.deriv_length));
	}
      
      END_CODE();

      return project(LatticePropagator(Gamma(G5) * fin), 
		     params.sink_mom, params.t_sink, params.j_decay);
    }


    // Construct pion_1-(B1xNabla_T1) sequential source
    LatticePropagator
    MesPionB1xNablaT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
					const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesPionB1xNablaT1SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // \f$\Gamma_f \equiv \gamma_4\gamma_5\epsilon_{ijk}\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (antiSymTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(antiSymTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * leftD(tmp,u,k,params.deriv_length));
	}
      
      END_CODE();

      return project(LatticePropagator(Gamma(1 << 3) * (Gamma(G5) * fin)), 
		     params.sink_mom, params.t_sink, params.j_decay);
    }


    // Construct pion_1-(A0_2xD_T2) sequential source
    LatticePropagator
    MesPionA02xDT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				     const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin;

      check1Args("MesPionA02xDT2SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // \f$\Gamma_f \equiv \gamma_4 D_i\f$  
      fin = Gamma(1 << 3) * leftD(tmp,u,params.deriv_dir,params.deriv_length);
      
      END_CODE();

      return project(fin, params.sink_mom, params.t_sink, params.j_decay);
    }


    // Construct pion_1-(A1xD_A2) sequential source
    LatticePropagator
    MesPionA1xDA2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				    const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesPionA1xDA2SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // \f$\Gamma_f \equiv \gamma_5\gamma_i D_i\f$  
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * leftD(tmp,u,k,params.deriv_length);
      
      END_CODE();

      return project(LatticePropagator(Gamma(G5) * fin), 
		     params.sink_mom, params.t_sink, params.j_decay);
    }


    // Construct pion_1-(A1xD_E) sequential source
    LatticePropagator
    MesPionA1xDESeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesPionA1xDESeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // \f$\Gamma_f \equiv \gamma_5 S_{\alpha jk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (ETensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(ETensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * leftD(tmp,u,k,params.deriv_length));
	}
      
      END_CODE();

      return project(LatticePropagator(Gamma(G5) * fin), 
		     params.sink_mom, params.t_sink, params.j_decay);
    }


    // Construct pion_1-(A1xD_T1) sequential source
    LatticePropagator
    MesPionA1xDT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				    const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesPionA1xDT1SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (symTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(symTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * leftD(tmp,u,k,params.deriv_length));
	}
      
      END_CODE();

      return project(LatticePropagator(Gamma(G5) * fin), 
		     params.sink_mom, params.t_sink, params.j_decay);
    }

    
    // Construct pion_1-(A1xD_T2) sequential source
    LatticePropagator
    MesPionA1xDT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				    const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesPionA1xDT2SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // \f$\Gamma_f \equiv \gamma_5\epsilon_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (antiSymTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(antiSymTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * leftD(tmp,u,k,params.deriv_length));
	}
      
      END_CODE();

      return project(LatticePropagator(Gamma(1 << 3) * (Gamma(G5) * fin)), 
		     params.sink_mom, params.t_sink, params.j_decay);
    }


    // Construct pion_1-(B1xD_A2) sequential source
    LatticePropagator
    MesPionB1xDA2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				    const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesPionB1xDA2SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 \gamma_i D_i\f$  
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * leftD(tmp,u,k,params.deriv_length);

      END_CODE();

      return project(LatticePropagator(Gamma(1 << 3) * (Gamma(G5) * fin)), 
		     params.sink_mom, params.t_sink, params.j_decay);
    }


    // Construct pion_1-(B1xD_E) sequential source
    LatticePropagator
    MesPionB1xDESeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesPionB1xDESeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);
 
      // \f$\Gamma_f \equiv \gamma_4\gamma_5 S_{\alpha jk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (ETensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(ETensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * leftD(tmp,u,k,params.deriv_length));
	}

      END_CODE();

      return project(LatticePropagator(Gamma(1 << 3) * (Gamma(G5) * fin)), 
		     params.sink_mom, params.t_sink, params.j_decay);
    }


    // Construct pion_1-(B1xD_T1) sequential source
    LatticePropagator
    MesPionB1xDT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				    const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesPionB1xDT1SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 s_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (symTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(symTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * leftD(tmp,u,k,params.deriv_length));
	}

      END_CODE();

      return project(LatticePropagator(Gamma(1 << 3) * (Gamma(G5) * fin)), 
		     params.sink_mom, params.t_sink, params.j_decay);
    }


    // Construct pion_1-(B1xD_T2) sequential source
    LatticePropagator
    MesPionB1xDT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				    const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesPionB1xDT2SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 \epsilon_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (antiSymTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(antiSymTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * leftD(tmp,u,k,params.deriv_length));
	}

      END_CODE();

      return project(LatticePropagator(Gamma(1 << 3) * (Gamma(G5) * fin)), 
		     params.sink_mom, params.t_sink, params.j_decay);
    }


    // Construct pion_1-(RhoxD_A2) sequential source
    LatticePropagator
    MesPionRhoxDA2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				     const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesPionRhoxDA2SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // \f$\Gamma_f \equiv \gamma_i D_i\f$  
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * leftD(tmp,u,k,params.deriv_length);
      
      END_CODE();

      return project(fin, params.sink_mom, params.t_sink, params.j_decay);
    }


    //! Construct pion_1-(RhoxD_T1) sequential source
    LatticePropagator
    MesPionRhoxDT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				     const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesPionRhoxDT1SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // \f$\Gamma_f \equiv s_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (symTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(symTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * leftD(tmp,u,k,params.deriv_length));
	}
      
      END_CODE();

      return project(fin, params.sink_mom, params.t_sink, params.j_decay);
    }


    //! Construct pion_1-(RhoxD_T2) sequential source
    LatticePropagator
    MesPionRhoxDT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				     const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesPionRhoxDT2SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (antiSymTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(antiSymTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * leftD(tmp,u,k,params.deriv_length));
	}
      
      END_CODE();

      return project(fin, params.sink_mom, params.t_sink, params.j_decay);
    }


    // Construct pion_1-(PionxD_T2) sequential source
    LatticePropagator
    MesPionPionxDT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin;
      int G5 = Ns*Ns-1;

      check1Args("MesPionPionxDT2SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 D_i\f$  
      fin = Gamma(1 << 3) * (Gamma(G5) * leftD(tmp,u,params.deriv_dir,params.deriv_length));
      
      END_CODE();

      return project(fin, params.sink_mom, params.t_sink, params.j_decay);
    }

 
    //! Construct pion_1-(PionxB_T1) sequential source
    LatticePropagator
    MesPionPionxBT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin;
      int G5 = Ns*Ns-1;

      check1Args("MesPionPionxBT1SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // \f$\Gamma_f \equiv \gamma_5 B_i\f$  
      fin = Gamma(G5) * leftB(tmp,u,params.deriv_dir,params.deriv_length);
      
      END_CODE();

      return project(fin, params.sink_mom, params.t_sink, params.j_decay);
    }


    //! Construct pion_1-(RhoxB_T1) sequential source
    LatticePropagator
    MesPionRhoxBT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				     const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesPionRhoxBT1SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (antiSymTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(antiSymTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * leftB(tmp,u,k,params.deriv_length));
	}
      
      END_CODE();

      return project(fin, params.sink_mom, params.t_sink, params.j_decay);
    }


    //! Construct pion_1-(RhoxB_T2) sequential source
    LatticePropagator
    MesPionRhoxBT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				     const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesPionRhoxBT2SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // \f$\Gamma_f \equiv s_{ijk}\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (symTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(symTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * leftB(tmp,u,k,params.deriv_length));
	}
      
      END_CODE();

      return project(fin, params.sink_mom, params.t_sink, params.j_decay);
    }


    //! Construct pion_1-(A1xB_A1) sequential source
    LatticePropagator
    MesPionA1xBA1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				    const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesPionA1xBA1SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // \f$\Gamma_f \equiv \gamma_5 \gamma_i B_i\f$  
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * leftB(tmp,u,k,params.deriv_length);
      
      END_CODE();

      return project(LatticePropagator(Gamma(G5) * fin), 
		     params.sink_mom, params.t_sink, params.j_decay);
    }


    //! Construct pion_1-(RhoxB_T1) sequential source
    LatticePropagator
    MesPionA1xBT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				    const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesPionA1xBT1SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // \f$\Gamma_f \equiv \gamma_5 \epsilon_{ijk}\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (antiSymTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(antiSymTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * leftB(tmp,u,k,params.deriv_length));
	}
      
      END_CODE();

      return project(LatticePropagator(Gamma(G5) * fin), 
		     params.sink_mom, params.t_sink, params.j_decay);
    }


    //! Construct pion_1-(A1xB_T2) sequential source
    LatticePropagator
    MesPionA1xBT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				    const multi1d<LatticePropagator>& quark_propagators) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesPionA1xBT2SeqSrc", quark_propagators);

      LatticePropagator tmp = adj(quark_propagators[0]);

      // \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (symTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(symTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * leftB(tmp,u,k,params.deriv_length));
	}
      
      END_CODE();

      return project(LatticePropagator(Gamma(G5) * fin), 
		     params.sink_mom, params.t_sink, params.j_decay);
    }


    // Register all the possible deriv mesons
    bool registerAll(void) 
    {
      bool success = true;

      //! Register all the factories
      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-PIONxNABLA_T2"),
										    mesPionA1xDT2SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-A0xNABLA_T1"),
										    mesPionA0xNablaT1SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-A0_2xNABLA_T1"),
										    mesPionA02xNablaT1SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-RHOxNABLA_A1"),
										    mesPionRhoxNablaA1SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-RHOxNABLA_T1"),
										    mesPionRhoxNablaT1SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-RHOxNABLA_T2"),
										    mesPionRhoxNablaT2SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-A1xNABLA_A1"),
										    mesPionA1xNablaA1SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-A1xNABLA_T2"),
										    mesPionA1xNablaT2SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-A1xNABLA_E"),
										    mesPionA1xNablaESeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-B1xNABLA_T1"),
										    mesPionB1xNablaT1SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-A02xD_T2"),
										    mesPionA02xDT2SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-A1xD_E"),
										    mesPionA1xDESeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-A1xD_T1"),
										    mesPionA1xDT1SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-A1xD_T2"),
										    mesPionA1xDT2SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-B1xD_A2"),
										    mesPionB1xDA2SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-B1xD_E"),
										    mesPionB1xDESeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-B1xD_T1"),
										    mesPionB1xDT1SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-B1xD_T2"),
										    mesPionB1xDT2SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-RHOxD_A2"),
										    mesPionRhoxDA2SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-RHOxD_T1"),
										    mesPionRhoxDT1SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-RHOxD_T2"),
										    mesPionRhoxDT2SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-PIONxD_T2"),
										    mesPionPionxDT2SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-PIONxB_T1"),
										    mesPionPionxBT1SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-RHOxB_T1"),
										    mesPionRhoxBT1SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-RHOxB_T2"),
										    mesPionRhoxBT2SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-A1xB_A1"),
										    mesPionA1xBA1SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-A1xB_T1"),
										    mesPionA1xBT1SeqSrc);

      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-A1xB_T2"),
										    mesPionA1xBT2SeqSrc);

      return success;
    }

    const bool registered = registerAll();

  }  // end namespace DerivMesonSeqSourceEnv

}  // end namespace Chroma



  
