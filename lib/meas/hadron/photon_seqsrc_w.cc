// $Id: photon_seqsrc_w.cc,v 3.0 2006-04-03 04:59:00 edwards Exp $
/*! \file
 *  \brief Construct a photon sequential sources via LSZ reduction
 */

#include "meas/hadron/photon_seqsrc_w.h"
#include "meas/hadron/seqsource_factory_w.h"

namespace Chroma 
{

  // Read parameters
  void read(XMLReader& xml, const string& path, PhotonRhoSeqSourceEnv::Params& param)
  {
    PhotonRhoSeqSourceEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const PhotonRhoSeqSourceEnv::Params& param)
  {
    param.writeXML(xml, path);
  }





  //! Meson sequential sources
  /*! \ingroup hadron */
  namespace PhotonRhoSeqSourceEnv
  { 
    //! Anonymous namespace
    namespace
    {
      //! Construct a0-photon sequential source
      /*!
       * \ingroup hadron
       */
      HadronSeqSource<LatticePropagator>* mesA0PhotonSeqSrc(XMLReader& xml_in,
							    const std::string& path)
      {
	return new PhotonRhoSeqSource(Params(xml_in, path));
      }

    } // end anonymous namespace


    //! Initialize
    Params::Params()
    {
      Q_sq     = zero;
      c_sq     = 1.0;
      xi       = 1.0;
      pol_dir  = -1;
      j_decay  = -1;
      t_sink   = -1;
      t_sink_start = -1;
      t_sink_end   = -1;
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

      read(paramtop, "Q_sq", Q_sq);
      read(paramtop, "c_sq", c_sq);
      read(paramtop, "xi", xi);
      read(paramtop, "pol_dir", pol_dir);
      read(paramtop, "j_decay", j_decay);
      read(paramtop, "sink_mom", sink_mom);
      read(paramtop, "t_sink_start", t_sink_start);
      read(paramtop, "t_sink_end", t_sink_end);
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);

      write(xml, "Q_sq", Q_sq);
      write(xml, "c_sq", c_sq);
      write(xml, "xi", xi);
      write(xml, "pol_dir", pol_dir);
      write(xml, "j_decay", j_decay);
      write(xml, "t_sink_start", t_sink_start);
      write(xml, "t_sink_end", t_sink_end);
      write(xml, "t_sink", t_sink);
      write(xml, "sink_mom", sink_mom);
      pop(xml);
    }


    //! Construct the source
    LatticePropagator
    PhotonRhoSeqSource::operator()(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& quark_propagators) const
    {
      QDPIO::cout << "Photon sequential source " << endl;

      START_CODE();

      if (quark_propagators.size() != 1)
      {
	QDPIO::cerr << __func__ << ": expect only 1 prop" << endl;
	QDP_abort(1);
      }

      if (params.j_decay < 0 || params.j_decay >= Nd)
      {
	QDPIO::cerr << __func__ << ": j_decay out of bounds: j_decay = " << params.j_decay << endl;
	QDP_abort(1);
      }

      if (params.pol_dir < 0 || params.pol_dir >= Nd)
      {
	QDPIO::cerr << __func__ << ": pol_dir out of bounds: pol_dir = " << params.pol_dir << endl;
	QDP_abort(1);
      }

      // Initial sequential source
      LatticePropagator seq_src_tmp;
      {
	int G5 = Ns*Ns-1;
	LatticePropagator tmp = adj(quark_propagators[0]) * Gamma(G5) * Gamma(1 << params.pol_dir);
      
	// Now take hermitian conjugate and multiply on both sides with gamma_5 = Gamma(15)
	seq_src_tmp = Gamma(15) * adj(tmp) * Gamma(15);
      }
        
      // Compute 4-vector correction of sink phase and energy exp
      LatticeComplex exp_p_dot_x;
      LatticeBoolean mask;
      {
	LatticeReal p_dot_x = zero;
	multi1d<Real> pp_f(Nd-1);
	for(int mu=0, j=0; mu < Nd; mu++)
	{
	  if (mu != params.j_decay)
	  {
	    pp_f[j] = params.sink_mom[j] * twopi / Real(Layout::lattSize()[mu]);
	    
	    if (params.sink_mom[j] != 0)
	      p_dot_x += Layout::latticeCoordinate(mu) * pp_f[j];

	    j++;
	  }
	}
	
	// Spatial contribution is solely a phase
	exp_p_dot_x = cmplx(cos(p_dot_x),sin(p_dot_x));

	// Photon energy determined from virtuality, \f$Q_f^2 = c^2|\vec{p_f}|^2 - \omega_f^2\f$
	Real omega;
	{
	  Real norm2_ppf = zero;
	  for(int i=0; i < pp_f.size(); ++i)
	    norm2_ppf += pp_f[i] * pp_f[i];

	  Real c_sq      =  params.c_sq;
	  Real xi        =  params.xi;
	  Real xi_sq     =  xi * xi;

	  omega = sqrt( Real(c_sq/xi_sq)*norm2_ppf - params.Q_sq );
	}
	QDPIO::cout << __func__ << ": omega= " << omega << endl;

	// Multiply in exp from time dependence of 4-vector
	// Note positive sign of omega
	exp_p_dot_x *= exp(omega * Layout::latticeCoordinate(params.j_decay));

	// Zap any regions in time outside integration region
	mask = where((Layout::latticeCoordinate(params.j_decay) >= params.t_sink_start) &&
		     (Layout::latticeCoordinate(params.j_decay) <= params.t_sink_end),
		     true, false);
      }
      
      
      // Multiply in by exp(4-vector)
      LatticePropagator fin = where(mask, 
				    exp_p_dot_x * seq_src_tmp, 
				    LatticePropagator(zero));

      END_CODE();

      return fin;
    }


    // Register all the possible simple mesons
    bool registerAll(void) 
    {
      bool success = true;

      //! Register all the factories
      success &= Chroma::TheWilsonHadronSeqSourceFactory::Instance().registerObject(string("PION-PHOTON"), 
										    mesA0PhotonSeqSrc);

      return success;
    }

    const bool registered = registerAll();

  }  // end namespace PhotonRhoSeqSourceEnv

}  // end namespace Chroma


  
