// $Id: deriv_quark_displacement_s.cc,v 3.1 2006-11-17 02:55:11 edwards Exp $
/*! \file
 *  \brief Derivative displacements
 */

#include "chromabase.h"

#include "meas/smear/deriv_quark_displacement_s.h"
#include "meas/smear/quark_displacement_factory.h"
#include "meas/smear/quark_displacement_aggregate.h"

namespace Chroma 
{

  // Read parameters
  void read(XMLReader& xml, const string& path, StaggeredDerivQuarkDisplacementEnv::Params& param)
  {
    StaggeredDerivQuarkDisplacementEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const StaggeredDerivQuarkDisplacementEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  // Read parameters
  void read(XMLReader& xml, const string& path, StaggeredDerivQuarkDisplacementEnv::ParamsDir& param)
  {
    StaggeredDerivQuarkDisplacementEnv::ParamsDir tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const StaggeredDerivQuarkDisplacementEnv::ParamsDir& param)
  {
    param.writeXML(xml, path);
  }



  //! Meson sources
  /*! \ingroup sources */
  namespace StaggeredDerivQuarkDisplacementEnv
  { 
    //! Anonymous namespace
    namespace
    {

      //! Determine sign of plusminus
      /*!
       * \ingroup sources
       */
      int plusMinus(enum PlusMinus isign)
      {
	int is = 0;
	switch (isign)
	{
	case PLUS:
	  is = +1;
	  break;

	case MINUS:
	  is = -1;
	  break;

	default:
	  QDP_error_exit("illegal isign in plusminus");
	}
	return is;
      }


      //! Private displacement
      LatticeStaggeredPropagator displacement(const multi1d<LatticeColorMatrix>& u, 
					      const LatticeStaggeredPropagator& psi, 
					      int length, int dir)
      {
	if (dir < 0 || dir >= Nd)
	{
	  QDPIO::cerr << __func__ << ": invalid direction: dir=" << dir << endl;
	  QDP_abort(1);
	}

	LatticeStaggeredPropagator chi = psi;

	if (length > 0)
	{
	  for(int n = 0; n < length; ++n)
	  {
	    LatticeStaggeredPropagator tmp = shift(chi, FORWARD, dir);
	    chi = u[dir] * tmp;
	  }
	}
	else // If length = or < 0.  If length == 0, does nothing.
	{
	  for(int n = 0; n > length; --n)
	  {
	    LatticeStaggeredPropagator tmp = shift(adj(u[dir])*chi, BACKWARD, dir);
	    chi = tmp;
	  }
	}
	return chi;
      }


      //! Apply first deriv to the right onto source
      /*!
       * \ingroup sources
       *
       * \f$\nabla_\mu f(x) = U_\mu(x)f(x+\mu) - U_{-\mu}(x)f(x-\mu)\f$
       *
       * \return $\f \nabla_\mu F(x)\f$
       */
      LatticeStaggeredPropagator rightNabla(const LatticeStaggeredPropagator& F, 
					    const multi1d<LatticeColorMatrix>& u,
					    int mu, int length)
      {
	return displacement(u, F, length, mu) - displacement(u, F, -length, mu);
      }


      //-------------------- callback functions ---------------------------------------

      //! Construct (right Nabla) source
      QuarkDisplacement<LatticeStaggeredPropagator>* rightNablaDisplace(XMLReader& xml_in,
									const std::string& path)
      {
	return new RightNablaDisplace<LatticeStaggeredPropagator>(ParamsDir(xml_in, path));
      }

    } // end anonymous namespace


    //! Initialize
    Params::Params()
    {
      deriv_length = 0;
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

      read(paramtop, "DisplacementType",  displacement_type);
      read(paramtop, "deriv_length", deriv_length);
    }

    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      QDP::write(xml, "version", version);

      write(xml, "DisplacementType", displacement_type);
      write(xml, "deriv_length", deriv_length);

      pop(xml);
    }


    //! Initialize
    ParamsDir::ParamsDir()
    {
      deriv_dir = -1;
      deriv_length = 0;
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

      read(paramtop, "DisplacementType",  displacement_type);

      read(paramtop, "deriv_dir", deriv_dir);
      read(paramtop, "deriv_length", deriv_length);
    }


    // Writer
    void ParamsDir::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      QDP::write(xml, "version", version);

      write(xml, "DisplacementType", displacement_type);
      write(xml, "deriv_dir", deriv_dir);
      write(xml, "deriv_length", deriv_length);

      pop(xml);
    }


    // Construct (right Nabla) source
    // See corresponding .h file for doxygen comments
    template<>
    void
    RightNablaDisplace<LatticeStaggeredPropagator>::operator()(LatticeStaggeredPropagator& tmp,
							       const multi1d<LatticeColorMatrix>& u,
							       enum PlusMinus isign) const
    {
      START_CODE();

      LatticeStaggeredPropagator fin;
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \nabla_i\f$
      fin = rightNabla(tmp,u,params.deriv_dir,length);
      tmp = fin;

      END_CODE();
    }


    //! Local registration flag
    static bool registered = false;

    //! Register all the possible deriv mesons
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	//! Register all the factories
	success &= Chroma::TheStagPropDisplacementFactory::Instance().registerObject(string("NABLA-DERIV"),
										     rightNablaDisplace);

	registered = true;
      }
      return success;
    }
  }  // end namespace StaggeredDerivQuarkDisplacementEnv

}  // end namespace Chroma



  
