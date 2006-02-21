// $Id: deriv_quark_displacement_w.cc,v 1.1 2006-02-21 06:45:40 edwards Exp $
/*! \file
 *  \brief Derivative displacements
 */

#include "chromabase.h"

#include "meas/smear/deriv_quark_displacement_w.h"
#include "meas/smear/quark_displacement_factory.h"
#include "meas/smear/quark_displacement_aggregate.h"

#include "util/ferm/symtensor.h"
#include "util/ferm/antisymtensor.h"
#include "util/ferm/etensor.h"

namespace Chroma 
{

  // Read parameters
  void read(XMLReader& xml, const string& path, DerivQuarkDisplacementEnv::Params& param)
  {
    DerivQuarkDisplacementEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const DerivQuarkDisplacementEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  // Read parameters
  void read(XMLReader& xml, const string& path, DerivQuarkDisplacementEnv::ParamsDir& param)
  {
    DerivQuarkDisplacementEnv::ParamsDir tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const DerivQuarkDisplacementEnv::ParamsDir& param)
  {
    param.writeXML(xml, path);
  }




  //! Meson sources
  /*! \ingroup sources */
  namespace DerivQuarkDisplacementEnv
  { 
    //! Anonymous namespace
    namespace
    {

      //! Private displacement
      LatticePropagator displacement(const multi1d<LatticeColorMatrix>& u, 
				     const LatticePropagator& psi, 
				     int length, int dir)
      {
	if (dir < 0 || dir >= Nd)
	{
	  QDPIO::cerr << __func__ << ": invalid direction: dir=" << dir << endl;
	  QDP_abort(1);
	}

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


      //! Apply first deriv to the right onto source
      /*!
       * \ingroup sources
       *
       * \f$\nabla_\mu f(x) = U_\mu(x)f(x+\mu) - U_{-\mu}(x)f(x-\mu)\f$
       *
       * \return $\f \nabla_\mu F(x)\f$
       */
      LatticePropagator rightNabla(const LatticePropagator& F, 
				   const multi1d<LatticeColorMatrix>& u,
				   int mu, int length)
      {
	return displacement(u, F, -length, mu) - displacement(u, F, length, mu);
      }

      //! Apply "D_i" operator to the right onto source
      /*!
       * \ingroup sources
       *
       * \f$D_i = s_{ijk}\nabla_j\nabla_k\f$
       *
       * where  \f$s_{ijk} = +1 \quad\forall i\ne j, j\ne k, i \ne k\f$
       * 
       * \return $\f F(z,0) D_\mu\f$
       */
      LatticePropagator rightD(const LatticePropagator& F,
			       const multi1d<LatticeColorMatrix>& u,
			       int mu, int length)
      {
	LatticePropagator tmp = zero;

	// Slow implementation - to speed up could compute once the \nabla_j deriv
	for(int j=0; j < 3; ++j)
	  for(int k=0; k < 3; ++k)
	  {
	    if (symTensor3d(mu,j,k) != 0)
	      tmp += rightNabla(rightNabla(F,u,j,length), u,k,length);
	  }

	return tmp;
      }

      //! Apply "B_i" operator to the right onto source
      /*!
       * \ingroup sources
       *
       * \f$B_i = \epsilon_{ijk}\nabla_j\nabla_k\f$
       *
       * \return $\f F(z,0) B_\mu\f$
       */
      LatticePropagator rightB(const LatticePropagator& F,
			       const multi1d<LatticeColorMatrix>& u,
			       int mu, int length)
      {
	LatticePropagator tmp = zero;

	// Slow implementation - to speed up could compute once the \nabla_j deriv
	for(int j=0; j < 3; ++j)
	  for(int k=0; k < 3; ++k)
	  {
	    if (antiSymTensor3d(mu,j,k) != 0)
	      tmp += Real(antiSymTensor3d(mu,j,k)) * rightNabla(rightNabla(F,u,j,length), u,k,length);
	  }

	return tmp;
      }



      //-------------------- callback functions ---------------------------------------

      //! Construct (PionxNabla_T1) source
      QuarkDisplacement<LatticePropagator>* mesPionxNablaT1Displace(XMLReader& xml_in,
								    const std::string& path)
      {
	return new MesPionxNablaT1Displace<LatticePropagator>(ParamsDir(xml_in, path));
      }

      //! Construct (A0xNabla_T1) source
      QuarkDisplacement<LatticePropagator>* mesA0xNablaT1Displace(XMLReader& xml_in,
								  const std::string& path)
      {
	return new MesA0xNablaT1Displace<LatticePropagator>(ParamsDir(xml_in, path));
      }

      //! Construct (A0_2xNabla_T1) source
      QuarkDisplacement<LatticePropagator>* mesA02xNablaT1Displace(XMLReader& xml_in,
								   const std::string& path)
      {
	return new MesA02xNablaT1Displace<LatticePropagator>(ParamsDir(xml_in, path));
      }

      //! Construct (RhoxNabla_A1) source
      QuarkDisplacement<LatticePropagator>* mesRhoxNablaA1Displace(XMLReader& xml_in,
								   const std::string& path)
      {
	return new MesRhoxNablaA1Displace<LatticePropagator>(Params(xml_in, path));
      }

      //! Construct (RhoxNabla_T1) source
      QuarkDisplacement<LatticePropagator>* mesRhoxNablaT1Displace(XMLReader& xml_in,
								   const std::string& path)
      {
	return new MesRhoxNablaT1Displace<LatticePropagator>(ParamsDir(xml_in, path));
      }

      //! Construct (RhoxNabla_T2) source
      QuarkDisplacement<LatticePropagator>* mesRhoxNablaT2Displace(XMLReader& xml_in,
								   const std::string& path)
      {
	return new MesRhoxNablaT2Displace<LatticePropagator>(ParamsDir(xml_in, path));
      }

      //! Construct (A1xNabla_A1) source
      QuarkDisplacement<LatticePropagator>* mesA1xNablaA1Displace(XMLReader& xml_in,
								  const std::string& path)
      {
	return new MesA1xNablaA1Displace<LatticePropagator>(Params(xml_in, path));
      }

      //! Construct (A1xNabla_T2) source
      QuarkDisplacement<LatticePropagator>* mesA1xNablaT2Displace(XMLReader& xml_in,
								  const std::string& path)
      {
	return new MesA1xNablaT2Displace<LatticePropagator>(ParamsDir(xml_in, path));
      }

      //! Construct (A1xNabla_E) source
      QuarkDisplacement<LatticePropagator>* mesA1xNablaEDisplace(XMLReader& xml_in,
								 const std::string& path)
      {
	return new MesA1xNablaEDisplace<LatticePropagator>(ParamsDir(xml_in, path));
      }

      //! Construct (B1xNabla_T1) source
      QuarkDisplacement<LatticePropagator>* mesB1xNablaT1Displace(XMLReader& xml_in,
								  const std::string& path)
      {
	return new MesB1xNablaT1Displace<LatticePropagator>(ParamsDir(xml_in, path));
      }

      //! Construct (A0_2xD_T2) source
      QuarkDisplacement<LatticePropagator>* mesA02xDT2Displace(XMLReader& xml_in,
							       const std::string& path)
      {
	return new MesA02xDT2Displace<LatticePropagator>(ParamsDir(xml_in, path));
      }

      //! Construct (A1xD_A2) source
      QuarkDisplacement<LatticePropagator>* mesA1xDA2Displace(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesA1xDA2Displace<LatticePropagator>(Params(xml_in, path));
      }

      //! Construct (A1xD_E) source
      QuarkDisplacement<LatticePropagator>* mesA1xDEDisplace(XMLReader& xml_in,
							     const std::string& path)
      {
	return new MesA1xDEDisplace<LatticePropagator>(ParamsDir(xml_in, path));
      }

      //! Construct (A1xD_T1) source
      QuarkDisplacement<LatticePropagator>* mesA1xDT1Displace(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesA1xDT1Displace<LatticePropagator>(ParamsDir(xml_in, path));
      }

      //! Construct (A1xD_T2) source
      QuarkDisplacement<LatticePropagator>* mesA1xDT2Displace(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesA1xDT2Displace<LatticePropagator>(ParamsDir(xml_in, path));
      }

      //! Construct (B1xD_A2) source
      QuarkDisplacement<LatticePropagator>* mesB1xDA2Displace(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesB1xDA2Displace<LatticePropagator>(Params(xml_in, path));
      }

      //! Construct (B1xD_E) source
      QuarkDisplacement<LatticePropagator>* mesB1xDEDisplace(XMLReader& xml_in,
							     const std::string& path)
      {
	return new MesB1xDEDisplace<LatticePropagator>(ParamsDir(xml_in, path));
      }

      //! Construct (B1xD_T1) source
      QuarkDisplacement<LatticePropagator>* mesB1xDT1Displace(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesB1xDT1Displace<LatticePropagator>(ParamsDir(xml_in, path));
      }

      //! Construct (B1xD_T2) source
      QuarkDisplacement<LatticePropagator>* mesB1xDT2Displace(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesB1xDT2Displace<LatticePropagator>(ParamsDir(xml_in, path));
      }

      //! Construct (RhoxD_A2) source
      QuarkDisplacement<LatticePropagator>* mesRhoxDA2Displace(XMLReader& xml_in,
							       const std::string& path)
      {
	return new MesRhoxDA2Displace<LatticePropagator>(Params(xml_in, path));
      }

      //! Construct (RhoxD_T1) source
      QuarkDisplacement<LatticePropagator>* mesRhoxDT1Displace(XMLReader& xml_in,
							       const std::string& path)
      {
	return new MesRhoxDT1Displace<LatticePropagator>(ParamsDir(xml_in, path));
      }

      //! Construct (RhoxD_T2) source
      QuarkDisplacement<LatticePropagator>* mesRhoxDT2Displace(XMLReader& xml_in,
							       const std::string& path)
      {
	return new MesRhoxDT2Displace<LatticePropagator>(ParamsDir(xml_in, path));
      }

      //! Construct (PionxD_T2) source
      QuarkDisplacement<LatticePropagator>* mesPionxDT2Displace(XMLReader& xml_in,
								const std::string& path)
      {
	return new MesPionxDT2Displace<LatticePropagator>(ParamsDir(xml_in, path));
      }

      //! Construct (PionxB_T1) source
      QuarkDisplacement<LatticePropagator>* mesPionxBT1Displace(XMLReader& xml_in,
								const std::string& path)
      {
	return new MesPionxBT1Displace<LatticePropagator>(ParamsDir(xml_in, path));
      }

      //! Construct (RhoxB_T1) source
      QuarkDisplacement<LatticePropagator>* mesRhoxBT1Displace(XMLReader& xml_in,
							       const std::string& path)
      {
	return new MesRhoxBT1Displace<LatticePropagator>(ParamsDir(xml_in, path));
      }

      //! Construct (RhoxB_T2) source
      QuarkDisplacement<LatticePropagator>* mesRhoxBT2Displace(XMLReader& xml_in,
							       const std::string& path)
      {
	return new MesRhoxBT2Displace<LatticePropagator>(ParamsDir(xml_in, path));
      }

      //! Construct (A1xB_A1) source
      QuarkDisplacement<LatticePropagator>* mesA1xBA1Displace(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesA1xBA1Displace<LatticePropagator>(Params(xml_in, path));
      }

      //! Construct (A1xB_T1) source
      QuarkDisplacement<LatticePropagator>* mesA1xBT1Displace(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesA1xBT1Displace<LatticePropagator>(ParamsDir(xml_in, path));
      }

      //! Construct (A1xB_T2) source
      QuarkDisplacement<LatticePropagator>* mesA1xBT2Displace(XMLReader& xml_in,
							      const std::string& path)
      {
	return new MesA1xBT2Displace<LatticePropagator>(ParamsDir(xml_in, path));
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


    // Construct (PionxNabla_T1) source
    // See corresponding .h file for doxygen comments
    template<>
    void
    MesPionxNablaT1Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
							   const multi1d<LatticeColorMatrix>& u,
							   enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin;
      const int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \nabla_i\f$
      fin = Gamma(G5) * rightNabla(tmp,u,params.deriv_dir,params.deriv_length);
      tmp = fin;

      END_CODE();
    }


    // Construct (A0xNabla_T1) source
    // See corresponding .h file for doxygen comments
    template<>
    void
    MesA0xNablaT1Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
							 const multi1d<LatticeColorMatrix>& u,
							 enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin;

      // \f$\Gamma_f \equiv \nabla_i\f$
      fin = rightNabla(tmp,u,params.deriv_dir,params.deriv_length);
      tmp = fin;
      
      END_CODE();
    }


    // Construct (A0_2xNabla_T1) source
    // See corresponding .h file for doxygen comments
    template<>
    void
    MesA02xNablaT1Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
							  const multi1d<LatticeColorMatrix>& u,
							  enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin;

      // \f$\Gamma_f \equiv \gamma_4 \nabla_i\f$
      fin = Gamma(1 << 3) * rightNabla(tmp,u,params.deriv_dir,params.deriv_length);
      tmp = fin;
      
      END_CODE();
    }


    // Construct (RhoxNabla_A1) source
    template<>
    void
    MesRhoxNablaA1Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
							  const multi1d<LatticeColorMatrix>& u,
							  enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin = zero;

      // \f$\Gamma_f \equiv \gamma_i\nabla_i\f$
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * rightNabla(tmp,u,k,params.deriv_length);
      
      tmp = fin;

      END_CODE();
    }


    // Construct (RhoxNabla_T1) source
    template<>
    void
    MesRhoxNablaT1Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
							  const multi1d<LatticeColorMatrix>& u,
							  enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin = zero;

      // \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (antiSymTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(antiSymTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightNabla(tmp,u,k,params.deriv_length));
	}
      
      tmp = fin;

      END_CODE();
    }


    // Construct (RhoxNabla_T2) source
    template<>
    void
    MesRhoxNablaT2Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
							  const multi1d<LatticeColorMatrix>& u,
							  enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin = zero;

      // \f$\Gamma_f \equiv s_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (symTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(symTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightNabla(tmp,u,k,params.deriv_length));
	}
      
      tmp = fin;

      END_CODE();
    }


    // Construct (A1xNabla_A1) source
    template<>
    void
    MesA1xNablaA1Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
							 const multi1d<LatticeColorMatrix>& u,
							 enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_5\gamma_i \nabla_i\f$  
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * rightNabla(tmp,u,k,params.deriv_length);
      
      tmp = fin;

      END_CODE();
    }


    // Construct (A1xNabla_T2) source
    template<>
    void
    MesA1xNablaT2Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
							 const multi1d<LatticeColorMatrix>& u,
							 enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (symTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(symTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightNabla(tmp,u,k,params.deriv_length));
	}
      
      tmp = fin;

      END_CODE();
    }


    // Construct (A1xNabla_E) source
    template<>
    void
    MesA1xNablaEDisplace<LatticePropagator>::operator()(LatticePropagator& tmp,
							const multi1d<LatticeColorMatrix>& u,
							enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_5 S_{\alpha jk}\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (ETensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(ETensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightNabla(tmp,u,k,params.deriv_length));
	}
      
      tmp = fin;

      END_CODE();
    }


    // Construct (B1xNabla_T1) source
    template<>
    void
    MesB1xNablaT1Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
							 const multi1d<LatticeColorMatrix>& u,
							 enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_4\gamma_5\epsilon_{ijk}\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (antiSymTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(antiSymTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightD(tmp,u,k,params.deriv_length));
	}
      
      tmp = fin;

      END_CODE();
    }


    // Construct (A0_2xD_T2) source
    template<>
    void
    MesA02xDT2Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
						      const multi1d<LatticeColorMatrix>& u,
						      enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin;

      // \f$\Gamma_f \equiv \gamma_4 D_i\f$  
      fin = Gamma(1 << 3) * rightD(tmp,u,params.deriv_dir,params.deriv_length);
      
      tmp = fin;

      END_CODE();
    }


    // Construct (A1xD_A2) source
    template<>
    void
    MesA1xDA2Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
						     const multi1d<LatticeColorMatrix>& u,
						     enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_5\gamma_i D_i\f$  
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * rightD(tmp,u,k,params.deriv_length);
      
      tmp = fin;

      END_CODE();
    }


    // Construct (A1xD_E) source
    template<>
    void
    MesA1xDEDisplace<LatticePropagator>::operator()(LatticePropagator& tmp,
						    const multi1d<LatticeColorMatrix>& u,
						    enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_5 S_{\alpha jk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (ETensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(ETensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightD(tmp,u,k,params.deriv_length));
	}
      
      tmp = fin;

      END_CODE();
    }


    // Construct (A1xD_T1) source
    template<>
    void
    MesA1xDT1Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
						     const multi1d<LatticeColorMatrix>& u,
						     enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (symTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(symTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightD(tmp,u,k,params.deriv_length));
	}
      
      tmp = fin;

      END_CODE();
    }

    
    // Construct (A1xD_T2) source
    template<>
    void
    MesA1xDT2Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
						     const multi1d<LatticeColorMatrix>& u,
						     enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_5\epsilon_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (antiSymTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(antiSymTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightD(tmp,u,k,params.deriv_length));
	}
      
      tmp = fin;

      END_CODE();
    }


    // Construct (B1xD_A2) source
    template<>
    void
    MesB1xDA2Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
						     const multi1d<LatticeColorMatrix>& u,
						     enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 \gamma_i D_i\f$  
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * rightD(tmp,u,k,params.deriv_length);

      tmp = fin;

      END_CODE();
    }


    // Construct (B1xD_E) source
    template<>
    void
    MesB1xDEDisplace<LatticePropagator>::operator()(LatticePropagator& tmp,
						    const multi1d<LatticeColorMatrix>& u,
						    enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 S_{\alpha jk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (ETensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(ETensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightD(tmp,u,k,params.deriv_length));
	}

      tmp = fin;

      END_CODE();
    }


    // Construct (B1xD_T1) source
    template<>
    void
    MesB1xDT1Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
						     const multi1d<LatticeColorMatrix>& u,
						     enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 s_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (symTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(symTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightD(tmp,u,k,params.deriv_length));
	}

      tmp = fin;

      END_CODE();
    }


    // Construct (B1xD_T2) source
    template<>
    void
    MesB1xDT2Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
						     const multi1d<LatticeColorMatrix>& u,
						     enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 \epsilon_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (antiSymTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(antiSymTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightD(tmp,u,k,params.deriv_length));
	}

      tmp = fin;

      END_CODE();
    }


    // Construct (RhoxD_A2) source
    template<>
    void
    MesRhoxDA2Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
						      const multi1d<LatticeColorMatrix>& u,
						      enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_i D_i\f$  
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * rightD(tmp,u,k,params.deriv_length);
      
      tmp = fin;

      END_CODE();
    }


    //! Construct (RhoxD_T1) source
    template<>
    void
    MesRhoxDT1Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
						      const multi1d<LatticeColorMatrix>& u,
						      enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv s_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (symTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(symTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightD(tmp,u,k,params.deriv_length));
	}
      
      tmp = fin;

      END_CODE();
    }


    //! Construct (RhoxD_T2) source
    template<>
    void
    MesRhoxDT2Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
						      const multi1d<LatticeColorMatrix>& u,
						      enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (antiSymTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(antiSymTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightD(tmp,u,k,params.deriv_length));
	}
      
      tmp = fin;

      END_CODE();
    }


    // Construct (PionxD_T2) source
    template<>
    void
    MesPionxDT2Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
						       const multi1d<LatticeColorMatrix>& u,
						       enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 D_i\f$  
      fin = Gamma(1 << 3) * (Gamma(G5) * rightD(tmp,u,params.deriv_dir,params.deriv_length));
      
      tmp = fin;

      END_CODE();
    }

 
    //! Construct (PionxB_T1) source
    template<>
    void
    MesPionxBT1Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
						       const multi1d<LatticeColorMatrix>& u,
						       enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_5 B_i\f$  
      fin = Gamma(G5) * rightB(tmp,u,params.deriv_dir,params.deriv_length);
      
      tmp = fin;

      END_CODE();
    }


    //! Construct (RhoxB_T1) source
    template<>
    void
    MesRhoxBT1Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
						      const multi1d<LatticeColorMatrix>& u,
						      enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (antiSymTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(antiSymTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightB(tmp,u,k,params.deriv_length));
	}
      
      tmp = fin;

      END_CODE();
    }


    //! Construct (RhoxB_T2) source
    template<>
    void
    MesRhoxBT2Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
						      const multi1d<LatticeColorMatrix>& u,
						      enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv s_{ijk}\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (symTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(symTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightB(tmp,u,k,params.deriv_length));
	}
      
      tmp = fin;

      END_CODE();
    }


    //! Construct (A1xB_A1) source
    template<>
    void
    MesA1xBA1Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
						     const multi1d<LatticeColorMatrix>& u,
						     enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_5 \gamma_i B_i\f$  
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * rightB(tmp,u,k,params.deriv_length);
      
      tmp = fin;

      END_CODE();
    }


    //! Construct (RhoxB_T1) source
    template<>
    void
    MesA1xBT1Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
						     const multi1d<LatticeColorMatrix>& u,
						     enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_5 \epsilon_{ijk}\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (antiSymTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(antiSymTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightB(tmp,u,k,params.deriv_length));
	}
      
      tmp = fin;

      END_CODE();
    }


    //! Construct (A1xB_T2) source
    template<>
    void
    MesA1xBT2Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
						     const multi1d<LatticeColorMatrix>& u,
						     enum PlusMinus isign) const
    {
      START_CODE();
      
      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (symTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(symTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightB(tmp,u,k,params.deriv_length));
	}
      
      tmp = fin;

      END_CODE();
    }


    // Register all the possible deriv mesons
    bool registerAll(void) 
    {
      bool foo = true;

      //! Register all the factories
      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-PIONxNABLA_T1"),
									   mesPionxNablaT1Displace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-A0xNABLA_T1"),
									   mesA0xNablaT1Displace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-A0_2xNABLA_T1"),
									   mesA02xNablaT1Displace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-RHOxNABLA_A1"),
									   mesRhoxNablaA1Displace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-RHOxNABLA_T1"),
									   mesRhoxNablaT1Displace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-RHOxNABLA_T2"),
									   mesRhoxNablaT2Displace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-A1xNABLA_A1"),
									   mesA1xNablaA1Displace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-A1xNABLA_T2"),
									   mesA1xNablaT2Displace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-A1xNABLA_E"),
									   mesA1xNablaEDisplace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-B1xNABLA_T1"),
									   mesB1xNablaT1Displace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-A0_2xD_T2"),
									   mesA02xDT2Displace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-A1xD_A2"),
									   mesA1xDA2Displace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-A1xD_E"),
									   mesA1xDEDisplace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-A1xD_T1"),
									   mesA1xDT1Displace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-A1xD_T2"),
									   mesA1xDT2Displace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-B1xD_A2"),
									   mesB1xDA2Displace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-B1xD_E"),
									   mesB1xDEDisplace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-B1xD_T1"),
									   mesB1xDT1Displace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-B1xD_T2"),
									   mesB1xDT2Displace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-RHOxD_A2"),
									   mesRhoxDA2Displace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-RHOxD_T1"),
									   mesRhoxDT1Displace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-RHOxD_T2"),
									   mesRhoxDT2Displace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-PIONxD_T2"),
									   mesPionxDT2Displace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-PIONxB_T1"),
									   mesPionxBT1Displace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-RHOxB_T1"),
									   mesRhoxBT1Displace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-RHOxB_T2"),
									   mesRhoxBT2Displace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-A1xB_A1"),
									   mesA1xBA1Displace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-A1xB_T1"),
									   mesA1xBT1Displace);

      foo &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("DISPLACE-A1xB_T2"),
									   mesA1xBT2Displace);

      return foo;
    }

    const bool registered = registerAll();

  }  // end namespace DerivQuarkDisplacementEnv

}  // end namespace Chroma



  
