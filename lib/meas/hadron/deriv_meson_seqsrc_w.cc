// $Id: deriv_meson_seqsrc_w.cc,v 3.1 2006-11-27 04:33:35 edwards Exp $
/*! \file
 *  \brief Construct meson sequential sources.
 */

#include "meas/hadron/deriv_meson_seqsrc_w.h"
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

      //! Apply first deriv (nabla) to the right onto source
      /*!
       * \f$\nabla_\mu f(x) = U_\mu(x)f(x+\mu) - U_{-\mu}(x)f(x-\mu)\f$
       *
       * \return $\f \nabla_\mu Fx,0) = U_\mu(x) F(x+\mu)  - U_{x-\mu}^\dag(x-\mu) F(x-\mu)\f$
       */
      LatticePropagator 
      DerivMesonSeqSourceBase::rightNabla(const LatticePropagator& F, 
					  const multi1d<LatticeColorMatrix>& u,
					  int mu) const
      {
	return displacement(u, F, getDerivLength(), mu) - displacement(u, F, -getDerivLength(), mu);
      }

      //! Apply left and right "nabla_i" onto the source
      /*!
       * \f$\nabla_\mu f(x) = U_\mu(x)f(x+\mu) - U_{-\mu}(x)f(x-\mu)\f$
       *
       * \return $\f \nabla_\mu Fx,0) = U_\mu(x) F(x+\mu)  - U_{x-\mu}^\dag(x-\mu) F(x-\mu)\f$
       */
      LatticePropagator 
      DerivMesonSeqSourceBase::derivNabla(const LatticePropagator& F, 
					  const multi1d<LatticeColorMatrix>& u,
					  int mu) const
      {
	LatticeComplex ph = conj(phases());
	return ph*rightNabla(F, u, mu) + rightNabla(ph*F, u, mu);
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
      DerivMesonSeqSourceBase::derivD(const LatticePropagator& F,
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
	      tmp += ph * rightNabla(rightNabla(F,u,j), u, k);
	      tmp += rightNabla(rightNabla(ph * F,u,j), u, k);
	      tmp += Real(2)*rightNabla(ph * rightNabla(F,u,j), u, k);
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
      DerivMesonSeqSourceBase::derivB(const LatticePropagator& F,
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
	      tmp -= Real(antiSymTensor3d(mu,j,k)) * rightNabla(rightNabla(ph*F,u,j), u, k);
	  }

	return tmp;
      }



      //-------------------- callback functions ---------------------------------------

      //! Construct a0-(PionxNabla_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0PionxNablaT1SeqSrc(XMLReader& xml_in,
								  const std::string& path)
      {
	return new MesA0PionxNablaT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(A0xNabla_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A0xNablaT1SeqSrc(XMLReader& xml_in,
								const std::string& path)
      {
	return new MesA0A0xNablaT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(A0_2xNabla_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A02xNablaT1SeqSrc(XMLReader& xml_in,
								 const std::string& path)
      {
	return new MesA0A02xNablaT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(RhoxNabla_A1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0RhoxNablaA1SeqSrc(XMLReader& xml_in,
								 const std::string& path)
      {
	return new MesA0RhoxNablaA1SeqSrc(Params(xml_in, path));
      }

      //! Construct a0-(RhoxNabla_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0RhoxNablaT1SeqSrc(XMLReader& xml_in,
								 const std::string& path)
      {
	return new MesA0RhoxNablaT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(RhoxNabla_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0RhoxNablaT2SeqSrc(XMLReader& xml_in,
								 const std::string& path)
      {
	return new MesA0RhoxNablaT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(A1xNabla_A1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A1xNablaA1SeqSrc(XMLReader& xml_in,
								const std::string& path)
      {
	return new MesA0A1xNablaA1SeqSrc(Params(xml_in, path));
      }

      //! Construct a0-(A1xNabla_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A1xNablaT2SeqSrc(XMLReader& xml_in,
								const std::string& path)
      {
	return new MesA0A1xNablaT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(A1xNabla_E) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A1xNablaESeqSrc(XMLReader& xml_in,
							       const std::string& path)
      {
	return new MesA0A1xNablaESeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(B1xNabla_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0B1xNablaT1SeqSrc(XMLReader& xml_in,
								const std::string& path)
      {
	return new MesA0B1xNablaT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(A0_2xD_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A02xDT2SeqSrc(XMLReader& xml_in,
							     const std::string& path)
      {
	return new MesA0A02xDT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(A1xD_A2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A1xDA2SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new MesA0A1xDA2SeqSrc(Params(xml_in, path));
      }

      //! Construct a0-(A1xD_E) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A1xDESeqSrc(XMLReader& xml_in,
							   const std::string& path)
      {
	return new MesA0A1xDESeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(A1xD_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A1xDT1SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new MesA0A1xDT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(A1xD_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A1xDT2SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new MesA0A1xDT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(B1xD_A2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0B1xDA2SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new MesA0B1xDA2SeqSrc(Params(xml_in, path));
      }

      //! Construct a0-(B1xD_E) sequential source
      HadronSeqSource<LatticePropagator>* mesA0B1xDESeqSrc(XMLReader& xml_in,
							   const std::string& path)
      {
	return new MesA0B1xDESeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(B1xD_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0B1xDT1SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new MesA0B1xDT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(B1xD_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0B1xDT2SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new MesA0B1xDT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(RhoxD_A2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0RhoxDA2SeqSrc(XMLReader& xml_in,
							     const std::string& path)
      {
	return new MesA0RhoxDA2SeqSrc(Params(xml_in, path));
      }

      //! Construct a0-(RhoxD_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0RhoxDT1SeqSrc(XMLReader& xml_in,
							     const std::string& path)
      {
	return new MesA0RhoxDT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(RhoxD_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0RhoxDT2SeqSrc(XMLReader& xml_in,
							     const std::string& path)
      {
	return new MesA0RhoxDT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(PionxD_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0PionxDT2SeqSrc(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesA0PionxDT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(PionxB_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0PionxBT1SeqSrc(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesA0PionxBT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(RhoxB_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0RhoxBT1SeqSrc(XMLReader& xml_in,
							     const std::string& path)
      {
	return new MesA0RhoxBT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(RhoxB_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0RhoxBT2SeqSrc(XMLReader& xml_in,
							     const std::string& path)
      {
	return new MesA0RhoxBT2SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(A1xB_A1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A1xBA1SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new MesA0A1xBA1SeqSrc(Params(xml_in, path));
      }

      //! Construct a0-(A1xB_T1) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A1xBT1SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new MesA0A1xBT1SeqSrc(ParamsDir(xml_in, path));
      }

      //! Construct a0-(A1xB_T2) sequential source
      HadronSeqSource<LatticePropagator>* mesA0A1xBT2SeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new MesA0A1xBT2SeqSrc(ParamsDir(xml_in, path));
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


    // Construct a0-(PionxNabla_T1) sequential source
    // See corresponding .h file for doxygen comments
    LatticePropagator
    MesA0PionxNablaT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
					const multi1d<ForwardProp_t>& forward_headers,
					const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin;
      const int G5 = Ns*Ns-1;

      check1Args("MesA0PionxNablaT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5\nabla_i\f$
      // \f$\gamma_5 \Gamma_f^\dag \gamma_5\f$
      fin = Gamma(G5) * derivNabla(tmp,u,params.deriv_dir);
      
      END_CODE();

      return project(fin);
    }


    // Construct a0-(A0xNabla_T1) sequential source
    // See corresponding .h file for doxygen comments
    LatticePropagator
    MesA0A0xNablaT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<ForwardProp_t>& forward_headers,
				      const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin;

      check1Args("MesA0A0xNablaT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \nabla_i\f$
      fin = derivNabla(tmp,u,params.deriv_dir);
      
      END_CODE();

      return project(fin);
    }


    // Construct a0-(A0_2xNabla_T1) sequential source
    // See corresponding .h file for doxygen comments
    LatticePropagator
    MesA0A02xNablaT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				       const multi1d<ForwardProp_t>& forward_headers,
				       const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin;

      check1Args("MesA0A02xNablaT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];
      
      // \f$\Gamma_f \equiv \gamma_4 \nabla_i\f$
      fin = Gamma(1 << 3) * derivNabla(tmp,u,params.deriv_dir);
      
      END_CODE();

      return project(-fin);
    }


    // Construct a0-(RhoxNabla_A1) sequential source
    LatticePropagator
    MesA0RhoxNablaA1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				       const multi1d<ForwardProp_t>& forward_headers,
				       const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin = zero;

      check1Args("MesA0RhoxNablaA1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_i\nabla_i\f$
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * derivNabla(tmp,u,k);
      
      END_CODE();

      return project(-fin);
    }


    // Construct a0-(RhoxNabla_T1) sequential source
    LatticePropagator
    MesA0RhoxNablaT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				       const multi1d<ForwardProp_t>& forward_headers,
				       const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin = zero;

      check1Args("MesA0RhoxNablaT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int a = antiSymTensor3d(params.deriv_dir,j,k);
	  if (a != 0)
	    fin += Real(a) * (Gamma(1 << j) * derivNabla(tmp,u,k));
	}
      
      END_CODE();

      return project(-fin);
    }


    // Construct a0-(RhoxNabla_T2) sequential source
    LatticePropagator
    MesA0RhoxNablaT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				       const multi1d<ForwardProp_t>& forward_headers,
				       const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin = zero;

      check1Args("MesA0RhoxNablaT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv s_{ijk}\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int s = symTensor3d(params.deriv_dir,j,k);
	  if (s != 0)
	    fin += Real(s) * (Gamma(1 << j) * derivNabla(tmp,u,k));
	}
      
      END_CODE();

      return project(-fin);
    }


    // Construct a0-(A1xNabla_A1) sequential source
    LatticePropagator
    MesA0A1xNablaA1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<ForwardProp_t>& forward_headers,
				      const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0A1xNablaA1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5\gamma_i \nabla_i\f$  
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * derivNabla(tmp,u,k);
      
      END_CODE();

      return project(Gamma(G5) * fin);
    }


    // Construct a0-(A1xNabla_T2) sequential source
    LatticePropagator
    MesA0A1xNablaT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<ForwardProp_t>& forward_headers,
				      const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0A1xNablaT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int s = symTensor3d(params.deriv_dir,j,k);
	  if (s != 0)
	    fin += Real(s) * (Gamma(1 << j) * derivNabla(tmp,u,k));
	}
      
      END_CODE();

      return project(Gamma(G5) * fin);
    }


    // Construct a0-(A1xNabla_E) sequential source
    LatticePropagator
    MesA0A1xNablaESeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				     const multi1d<ForwardProp_t>& forward_headers,
				     const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0A1xNablaESeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5 S_{\alpha jk}\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int e = ETensor3d(params.deriv_dir,j,k);
	  if (e != 0)
	    fin += Real(e) * (Gamma(1 << j) * derivNabla(tmp,u,k));
	}
      
      END_CODE();

      return project(Gamma(G5) * fin);
    }


    // Construct a0-(B1xNabla_T1) sequential source
    LatticePropagator
    MesA0B1xNablaT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<ForwardProp_t>& forward_headers,
				      const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0B1xNablaT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_4\gamma_5\epsilon_{ijk}\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int a = antiSymTensor3d(params.deriv_dir,j,k);
	  if (a != 0)
	    fin += Real(a) * (Gamma(1 << j) * derivNabla(tmp,u,k));
	}
      
      END_CODE();

      return project(-(Gamma(1 << 3) * (Gamma(G5) * fin)));
    }


    // Construct a0-(A0_2xD_T2) sequential source
    LatticePropagator
    MesA0A02xDT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin;

      check1Args("MesA0A02xDT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_4 D_i\f$  
      fin = Gamma(1 << 3) * derivD(tmp,u,params.deriv_dir);
      
      END_CODE();

      return project(-fin);
    }


    // Construct a0-(A1xD_A2) sequential source
    LatticePropagator
    MesA0A1xDA2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0A1xDA2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5\gamma_i D_i\f$  
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * derivD(tmp,u,k);
      
      END_CODE();

      return project(Gamma(G5) * fin);
    }


    // Construct a0-(A1xD_E) sequential source
    LatticePropagator
    MesA0A1xDESeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<ForwardProp_t>& forward_headers,
				 const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0A1xDESeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5 S_{\alpha jk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int e = ETensor3d(params.deriv_dir,j,k);
	  if (e != 0)
	    fin += Real(e) * (Gamma(1 << j) * derivD(tmp,u,k));
	}
      
      END_CODE();

      return project(Gamma(G5) * fin);
    }


    // Construct a0-(A1xD_T1) sequential source
    LatticePropagator
    MesA0A1xDT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0A1xDT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int s = symTensor3d(params.deriv_dir,j,k);
	  if (s != 0)
	    fin += Real(s) * (Gamma(1 << j) * derivD(tmp,u,k));
	}
      
      END_CODE();

      return project(Gamma(G5) * fin);
    }

    
    // Construct a0-(A1xD_T2) sequential source
    LatticePropagator
    MesA0A1xDT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0A1xDT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5\epsilon_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int a = antiSymTensor3d(params.deriv_dir,j,k);
	  if (a != 0)
	    fin += Real(a) * (Gamma(1 << j) * derivD(tmp,u,k));
	}
      
      END_CODE();

      return project(Gamma(G5) * fin);
    }


    // Construct a0-(B1xD_A2) sequential source
    LatticePropagator
    MesA0B1xDA2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0B1xDA2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 \gamma_i D_i\f$  
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * derivD(tmp,u,k);

      END_CODE();

      return project(-(Gamma(1 << 3) * (Gamma(G5) * fin)));
    }


    // Construct a0-(B1xD_E) sequential source
    LatticePropagator
    MesA0B1xDESeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<ForwardProp_t>& forward_headers,
				 const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0B1xDESeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];
 
      // \f$\Gamma_f \equiv \gamma_4\gamma_5 S_{\alpha jk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int e = ETensor3d(params.deriv_dir,j,k);
	  if (e != 0)
	    fin += Real(e) * (Gamma(1 << j) * derivD(tmp,u,k));
	}

      END_CODE();

      return project(-(Gamma(1 << 3) * (Gamma(G5) * fin)));
    }


    // Construct a0-(B1xD_T1) sequential source
    LatticePropagator
    MesA0B1xDT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0B1xDT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 s_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int s = symTensor3d(params.deriv_dir,j,k);
	  if (s != 0)
	    fin += Real(s) * (Gamma(1 << j) * derivD(tmp,u,k));
	}

      END_CODE();

      return project(-(Gamma(1 << 3) * (Gamma(G5) * fin)));
    }


    // Construct a0-(B1xD_T2) sequential source
    LatticePropagator
    MesA0B1xDT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0B1xDT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 \epsilon_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int a = antiSymTensor3d(params.deriv_dir,j,k);
	  if (a != 0)
	    fin += Real(a) * (Gamma(1 << j) * derivD(tmp,u,k));
	}

      END_CODE();

      return project(-(Gamma(1 << 3) * (Gamma(G5) * fin)));
    }


    // Construct a0-(RhoxD_A2) sequential source
    LatticePropagator
    MesA0RhoxDA2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0RhoxDA2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_i D_i\f$  
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * derivD(tmp,u,k);
      
      END_CODE();

      return project(-fin);
    }


    //! Construct a0-(RhoxD_T1) sequential source
    LatticePropagator
    MesA0RhoxDT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0RhoxDT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv s_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int s = symTensor3d(params.deriv_dir,j,k);
	  if (s != 0)
	    fin += Real(s) * (Gamma(1 << j) * derivD(tmp,u,k));
	}
      
      END_CODE();

      return project(-fin);
    }


    //! Construct a0-(RhoxD_T2) sequential source
    LatticePropagator
    MesA0RhoxDT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0RhoxDT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int a = antiSymTensor3d(params.deriv_dir,j,k);
	  if (a != 0)
	    fin += Real(a) * (Gamma(1 << j) * derivD(tmp,u,k));
	}
      
      END_CODE();

      return project(-fin);
    }


    // Construct a0-(PionxD_T2) sequential source
    LatticePropagator
    MesA0PionxDT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				    const multi1d<ForwardProp_t>& forward_headers,
				    const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin;
      int G5 = Ns*Ns-1;

      check1Args("MesA0PionxDT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 D_i\f$  
      fin = Gamma(1 << 3) * (Gamma(G5) * derivD(tmp,u,params.deriv_dir));
      
      END_CODE();

      return project(fin);
    }

 
    //! Construct a0-(PionxB_T1) sequential source
    LatticePropagator
    MesA0PionxBT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				    const multi1d<ForwardProp_t>& forward_headers,
				    const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin;
      int G5 = Ns*Ns-1;

      check1Args("MesA0PionxBT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5 B_i\f$  
      fin = Gamma(G5) * derivB(tmp,u,params.deriv_dir);
      
      END_CODE();

      return project(fin);
    }


    //! Construct a0-(RhoxB_T1) sequential source
    LatticePropagator
    MesA0RhoxBT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0RhoxBT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int a = antiSymTensor3d(params.deriv_dir,j,k);
	  if (a != 0)
	    fin += Real(a) * (Gamma(1 << j) * derivB(tmp,u,k));
	}
      
      END_CODE();

      return project(-fin);
    }


    //! Construct a0-(RhoxB_T2) sequential source
    LatticePropagator
    MesA0RhoxBT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<ForwardProp_t>& forward_headers,
				   const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0RhoxBT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv s_{ijk}\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int s = symTensor3d(params.deriv_dir,j,k);
	  if (s != 0)
	    fin += Real(s) * (Gamma(1 << j) * derivB(tmp,u,k));
	}
      
      END_CODE();

      return project(-fin);
    }


    //! Construct a0-(A1xB_A1) sequential source
    LatticePropagator
    MesA0A1xBA1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0A1xBA1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5 \gamma_i B_i\f$  
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * derivB(tmp,u,k);
      
      END_CODE();

      return project(Gamma(G5) * fin);
    }


    //! Construct a0-(A1xB_T1) sequential source
    LatticePropagator
    MesA0A1xBT1SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0A1xBT1SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5 \epsilon_{ijk}\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int a = antiSymTensor3d(params.deriv_dir,j,k);
	  if (a != 0)
	    fin += Real(a) * (Gamma(1 << j) * derivB(tmp,u,k));
	}
      
      END_CODE();

      return project(Gamma(G5) * fin);
    }


    //! Construct a0-(A1xB_T2) sequential source
    LatticePropagator
    MesA0A1xBT2SeqSrc::operator()(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<ForwardProp_t>& forward_headers,
				  const multi1d<LatticePropagator>& quark_propagators)
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      check1Args("MesA0A1xBT2SeqSrc", quark_propagators);
      setTSrce(forward_headers);

      LatticePropagator tmp = quark_propagators[0];

      // \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  int s = symTensor3d(params.deriv_dir,j,k);
	  if (s != 0)
	    fin += Real(s) * (Gamma(1 << j) * derivB(tmp,u,k));
	}
      
      END_CODE();

      return project(Gamma(G5) * fin);
    }


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
	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-pionxNABLA_T1"),
										      mesA0PionxNablaT1SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a0xNABLA_T1"),
										      mesA0A0xNablaT1SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a0_2xNABLA_T1"),
										      mesA0A02xNablaT1SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rhoxNABLA_A1"),
										      mesA0RhoxNablaA1SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rhoxNABLA_T1"),
										      mesA0RhoxNablaT1SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rhoxNABLA_T2"),
										      mesA0RhoxNablaT2SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a1xNABLA_A1"),
										      mesA0A1xNablaA1SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a1xNABLA_T2"),
										      mesA0A1xNablaT2SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a1xNABLA_E"),
										      mesA0A1xNablaESeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-b1xNABLA_T1"),
										      mesA0B1xNablaT1SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a0_2xD_T2"),
										      mesA0A02xDT2SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a1xD_A2"),
										      mesA0A1xDA2SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a1xD_E"),
										      mesA0A1xDESeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a1xD_T1"),
										      mesA0A1xDT1SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a1xD_T2"),
										      mesA0A1xDT2SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-b1xD_A2"),
										      mesA0B1xDA2SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-b1xD_E"),
										      mesA0B1xDESeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-b1xD_T1"),
										      mesA0B1xDT1SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-b1xD_T2"),
										      mesA0B1xDT2SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rhoxD_A2"),
										      mesA0RhoxDA2SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rhoxD_T1"),
										      mesA0RhoxDT1SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rhoxD_T2"),
										      mesA0RhoxDT2SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-pionxD_T2"),
										      mesA0PionxDT2SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-pionxB_T1"),
										      mesA0PionxBT1SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rhoxB_T1"),
										      mesA0RhoxBT1SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-rhoxB_T2"),
										      mesA0RhoxBT2SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a1xB_A1"),
										      mesA0A1xBA1SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a1xB_T1"),
										      mesA0A1xBT1SeqSrc);

	success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("a0-a1xB_T2"),
										      mesA0A1xBT2SeqSrc);

	registered = true;
      }
      return success;
    }

  }  // end namespace DerivMesonSeqSourceEnv

}  // end namespace Chroma



  
