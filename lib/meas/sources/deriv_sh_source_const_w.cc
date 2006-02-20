// $Id: deriv_sh_source_const_w.cc,v 1.6 2006-02-20 22:51:17 edwards Exp $
/*! \file
 *  \brief Construct derivative source construction
 */

#include "chromabase.h"
#include "handle.h"

#include "meas/sources/source_const_factory.h"
#include "meas/sources/sh_source_const.h"
#include "meas/sources/srcfil.h"
#include "util/ferm/transf.h"

#include "meas/sources/deriv_sh_source_const_w.h"
#include "meas/sources/source_const_factory.h"

#include "meas/smear/quark_smearing_factory.h"
#include "meas/smear/quark_smearing_aggregate.h"

#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/link_smearing_factory.h"

#include "util/ferm/symtensor.h"
#include "util/ferm/antisymtensor.h"
#include "util/ferm/etensor.h"

namespace Chroma 
{

  // Read parameters
  void read(XMLReader& xml, const string& path, DerivShellSourceConstEnv::Params& param)
  {
    DerivShellSourceConstEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const DerivShellSourceConstEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  // Read parameters
  void read(XMLReader& xml, const string& path, DerivShellSourceConstEnv::ParamsDir& param)
  {
    DerivShellSourceConstEnv::ParamsDir tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const DerivShellSourceConstEnv::ParamsDir& param)
  {
    param.writeXML(xml, path);
  }




  //! Meson sources
  /*! \ingroup sources */
  namespace DerivShellSourceConstEnv
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
      QuarkSourceConstruction<LatticePropagator>* mesPionxNablaT1SrcConst(XMLReader& xml_in,
									  const std::string& path)
      {
	return new MesPionxNablaT1SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct (A0xNabla_T1) source
      QuarkSourceConstruction<LatticePropagator>* mesA0xNablaT1SrcConst(XMLReader& xml_in,
									const std::string& path)
      {
	return new MesA0xNablaT1SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct (A0_2xNabla_T1) source
      QuarkSourceConstruction<LatticePropagator>* mesA02xNablaT1SrcConst(XMLReader& xml_in,
									 const std::string& path)
      {
	return new MesA02xNablaT1SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct (RhoxNabla_A1) source
      QuarkSourceConstruction<LatticePropagator>* mesRhoxNablaA1SrcConst(XMLReader& xml_in,
									 const std::string& path)
      {
	return new MesRhoxNablaA1SrcConst(Params(xml_in, path));
      }

      //! Construct (RhoxNabla_T1) source
      QuarkSourceConstruction<LatticePropagator>* mesRhoxNablaT1SrcConst(XMLReader& xml_in,
									 const std::string& path)
      {
	return new MesRhoxNablaT1SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct (RhoxNabla_T2) source
      QuarkSourceConstruction<LatticePropagator>* mesRhoxNablaT2SrcConst(XMLReader& xml_in,
									 const std::string& path)
      {
	return new MesRhoxNablaT2SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct (A1xNabla_A1) source
      QuarkSourceConstruction<LatticePropagator>* mesA1xNablaA1SrcConst(XMLReader& xml_in,
									const std::string& path)
      {
	return new MesA1xNablaA1SrcConst(Params(xml_in, path));
      }

      //! Construct (A1xNabla_T2) source
      QuarkSourceConstruction<LatticePropagator>* mesA1xNablaT2SrcConst(XMLReader& xml_in,
									const std::string& path)
      {
	return new MesA1xNablaT2SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct (A1xNabla_E) source
      QuarkSourceConstruction<LatticePropagator>* mesA1xNablaESrcConst(XMLReader& xml_in,
								       const std::string& path)
      {
	return new MesA1xNablaESrcConst(ParamsDir(xml_in, path));
      }

      //! Construct (B1xNabla_T1) source
      QuarkSourceConstruction<LatticePropagator>* mesB1xNablaT1SrcConst(XMLReader& xml_in,
									const std::string& path)
      {
	return new MesB1xNablaT1SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct (A0_2xD_T2) source
      QuarkSourceConstruction<LatticePropagator>* mesA02xDT2SrcConst(XMLReader& xml_in,
								     const std::string& path)
      {
	return new MesA02xDT2SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct (A1xD_A2) source
      QuarkSourceConstruction<LatticePropagator>* mesA1xDA2SrcConst(XMLReader& xml_in,
								    const std::string& path)
      {
	return new MesA1xDA2SrcConst(Params(xml_in, path));
      }

      //! Construct (A1xD_E) source
      QuarkSourceConstruction<LatticePropagator>* mesA1xDESrcConst(XMLReader& xml_in,
								   const std::string& path)
      {
	return new MesA1xDESrcConst(ParamsDir(xml_in, path));
      }

      //! Construct (A1xD_T1) source
      QuarkSourceConstruction<LatticePropagator>* mesA1xDT1SrcConst(XMLReader& xml_in,
								    const std::string& path)
      {
	return new MesA1xDT1SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct (A1xD_T2) source
      QuarkSourceConstruction<LatticePropagator>* mesA1xDT2SrcConst(XMLReader& xml_in,
								    const std::string& path)
      {
	return new MesA1xDT2SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct (B1xD_A2) source
      QuarkSourceConstruction<LatticePropagator>* mesB1xDA2SrcConst(XMLReader& xml_in,
								    const std::string& path)
      {
	return new MesB1xDA2SrcConst(Params(xml_in, path));
      }

      //! Construct (B1xD_E) source
      QuarkSourceConstruction<LatticePropagator>* mesB1xDESrcConst(XMLReader& xml_in,
								   const std::string& path)
      {
	return new MesB1xDESrcConst(ParamsDir(xml_in, path));
      }

      //! Construct (B1xD_T1) source
      QuarkSourceConstruction<LatticePropagator>* mesB1xDT1SrcConst(XMLReader& xml_in,
								    const std::string& path)
      {
	return new MesB1xDT1SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct (B1xD_T2) source
      QuarkSourceConstruction<LatticePropagator>* mesB1xDT2SrcConst(XMLReader& xml_in,
								    const std::string& path)
      {
	return new MesB1xDT2SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct (RhoxD_A2) source
      QuarkSourceConstruction<LatticePropagator>* mesRhoxDA2SrcConst(XMLReader& xml_in,
								     const std::string& path)
      {
	return new MesRhoxDA2SrcConst(Params(xml_in, path));
      }

      //! Construct (RhoxD_T1) source
      QuarkSourceConstruction<LatticePropagator>* mesRhoxDT1SrcConst(XMLReader& xml_in,
								     const std::string& path)
      {
	return new MesRhoxDT1SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct (RhoxD_T2) source
      QuarkSourceConstruction<LatticePropagator>* mesRhoxDT2SrcConst(XMLReader& xml_in,
								     const std::string& path)
      {
	return new MesRhoxDT2SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct (PionxD_T2) source
      QuarkSourceConstruction<LatticePropagator>* mesPionxDT2SrcConst(XMLReader& xml_in,
								      const std::string& path)
      {
	return new MesPionxDT2SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct (PionxB_T1) source
      QuarkSourceConstruction<LatticePropagator>* mesPionxBT1SrcConst(XMLReader& xml_in,
								      const std::string& path)
      {
	return new MesPionxBT1SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct (RhoxB_T1) source
      QuarkSourceConstruction<LatticePropagator>* mesRhoxBT1SrcConst(XMLReader& xml_in,
								     const std::string& path)
      {
	return new MesRhoxBT1SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct (RhoxB_T2) source
      QuarkSourceConstruction<LatticePropagator>* mesRhoxBT2SrcConst(XMLReader& xml_in,
								     const std::string& path)
      {
	return new MesRhoxBT2SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct (A1xB_A1) source
      QuarkSourceConstruction<LatticePropagator>* mesA1xBA1SrcConst(XMLReader& xml_in,
								    const std::string& path)
      {
	return new MesA1xBA1SrcConst(Params(xml_in, path));
      }

      //! Construct (A1xB_T1) source
      QuarkSourceConstruction<LatticePropagator>* mesA1xBT1SrcConst(XMLReader& xml_in,
								    const std::string& path)
      {
	return new MesA1xBT1SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct (A1xB_T2) source
      QuarkSourceConstruction<LatticePropagator>* mesA1xBT2SrcConst(XMLReader& xml_in,
								    const std::string& path)
      {
	return new MesA1xBT2SrcConst(ParamsDir(xml_in, path));
      }

    } // end anonymous namespace


    //! Initialize
    Params::Params()
    {
      j_decay = -1;
      deriv_length = 0;
      t_srce.resize(Nd);
      t_srce = 0;
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

      read(paramtop, "SourceType",  source_type);

      read(paramtop, "deriv_length", deriv_length);

      {
	XMLReader xml_tmp(paramtop, "SmearingParam");
	std::ostringstream os;
	xml_tmp.print(os);
	read(xml_tmp, "wvf_kind", quark_smearing_type);
	quark_smearing = os.str();
      }

      if (paramtop.count("LinkSmearing") != 0)
      {
	XMLReader xml_tmp(paramtop, "LinkSmearing");
	std::ostringstream os;
	xml_tmp.print(os);
	read(xml_tmp, "LinkSmearingType", link_smearing_type);
	link_smearing = os.str();
      }

      read(paramtop, "t_srce", t_srce);
      read(paramtop, "j_decay",  j_decay);
    }

    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      QDP::write(xml, "version", version);

      write(xml, "SourceType", source_type);
      xml << quark_smearing;
      write(xml, "deriv_length", deriv_length);
      xml << link_smearing;
      pop(xml);

      write(xml, "t_srce",  t_srce);
      write(xml, "j_decay",  j_decay);

      pop(xml);
    }


    //! Initialize
    ParamsDir::ParamsDir()
    {
      deriv_dir = -1;
      deriv_length = 0;
      j_decay = -1;
      t_srce.resize(Nd);
      t_srce = 0;
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

      read(paramtop, "SourceType",  source_type);

      read(paramtop, "deriv_dir", deriv_dir);
      read(paramtop, "deriv_length", deriv_length);

      {
	XMLReader xml_tmp(paramtop, "SmearingParam");
	std::ostringstream os;
	xml_tmp.print(os);
	read(xml_tmp, "wvf_kind", quark_smearing_type);
	quark_smearing = os.str();
      }

      if (paramtop.count("LinkSmearing") != 0)
      {
	XMLReader xml_tmp(paramtop, "LinkSmearing");
	std::ostringstream os;
	xml_tmp.print(os);
	read(xml_tmp, "LinkSmearingType", link_smearing_type);
	link_smearing = os.str();
      }

      read(paramtop, "t_srce", t_srce);
      read(paramtop, "j_decay",  j_decay);
    }


    // Writer
    void ParamsDir::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      QDP::write(xml, "version", version);

      write(xml, "SourceType", source_type);
      xml << quark_smearing;
      write(xml, "deriv_dir", deriv_dir);
      write(xml, "deriv_length", deriv_length);
      xml << link_smearing;
      pop(xml);

      write(xml, "t_srce",  t_srce);
      write(xml, "j_decay",  j_decay);

      pop(xml);
    }


    // Construct the source, do the spin actions, and smear it
    template<>
    LatticePropagator
    DerivSourceConst<LatticePropagator>::operator()(const multi1d<LatticeColorMatrix>& u) const
    {
      QDPIO::cout << "Deriv Shell source" << endl;

      LatticePropagator quark_source;
      const DerivShellSourceConstEnv::ParamsDir& params = getParamsDir();

      try
      {
	//
	// Smear the gauge field if needed
	//
	multi1d<LatticeColorMatrix> u_smr = u;
	linkSmear(u_smr, std::string("/LinkSmearing"), params.link_smearing, params.link_smearing_type);


	//
	// Create the quark smearing object
	//
	std::istringstream  xml_s(params.quark_smearing);
	XMLReader  smeartop(xml_s);
	const string smear_path = "/SmearingParam";
	
	Handle< QuarkSmearing<LatticePropagator> >
	  quarkSmearing(ThePropSmearingFactory::Instance().createObject(params.quark_smearing_type,
									smeartop,
									smear_path));


	//
	// Create quark source
	//
	LatticePropagator pt_source;

	for(int color_source = 0; color_source < Nc; ++color_source)
	{
	  QDPIO::cout << "color = " << color_source << endl; 

	  LatticeColorVector src_color_vec = zero;

	  // Make a point source at coordinates t_srce
	  srcfil(src_color_vec, params.t_srce, color_source);

	  for(int spin_source = 0; spin_source < Ns; ++spin_source)
	  {
	    QDPIO::cout << "spin = " << spin_source << endl; 

	    // Insert a ColorVector into spin index spin_source
	    // This only overwrites sections, so need to initialize first
	    LatticeFermion chi = zero;

	    CvToFerm(src_color_vec, chi, spin_source);
      
	    /*
	     *  Move the source to the appropriate components
	     *  of quark source.
	     */
	    FermToProp(chi, pt_source, color_source, spin_source);
	  }
	}

	// Do the spin operations
	quark_source = deriv(u, pt_source);

	// do the smearing
	(*quarkSmearing)(quark_source, u_smr);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << "DerivSourceConst: Caught Exception in deriv source construction: " << e << endl;
	QDP_abort(1);
      }


      return quark_source;
    }




    // Construct (PionxNabla_T1) source
    // See corresponding .h file for doxygen comments
    LatticePropagator
    MesPionxNablaT1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				   const LatticePropagator& tmp) const
    {
      START_CODE();

      LatticePropagator fin;
      const int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \nabla_i\f$
      fin = Gamma(G5) * rightNabla(tmp,u,params.deriv_dir,params.deriv_length);
      
      END_CODE();

      return fin;
    }


    // Construct (A0xNabla_T1) source
    // See corresponding .h file for doxygen comments
    LatticePropagator
    MesA0xNablaT1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				 const LatticePropagator& tmp) const
    {
      START_CODE();

      LatticePropagator fin;

      // \f$\Gamma_f \equiv \nabla_i\f$
      fin = rightNabla(tmp,u,params.deriv_dir,params.deriv_length);
      
      END_CODE();

      return fin;
    }


    // Construct (A0_2xNabla_T1) source
    // See corresponding .h file for doxygen comments
    LatticePropagator
    MesA02xNablaT1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				  const LatticePropagator& tmp) const
    {
      START_CODE();

      LatticePropagator fin;

      // \f$\Gamma_f \equiv \gamma_4 \nabla_i\f$
      fin = Gamma(1 << 3) * rightNabla(tmp,u,params.deriv_dir,params.deriv_length);
      
      END_CODE();

      return fin;
    }


    // Construct (RhoxNabla_A1) source
    LatticePropagator
    MesRhoxNablaA1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				  const LatticePropagator& tmp) const
    {
      START_CODE();

      LatticePropagator fin = zero;

      // \f$\Gamma_f \equiv \gamma_i\nabla_i\f$
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * rightNabla(tmp,u,k,params.deriv_length);
      
      END_CODE();

      return fin;
    }


    // Construct (RhoxNabla_T1) source
    LatticePropagator
    MesRhoxNablaT1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				  const LatticePropagator& tmp) const
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
      
      END_CODE();

      return fin;
    }


    // Construct (RhoxNabla_T2) source
    LatticePropagator
    MesRhoxNablaT2SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				  const LatticePropagator& tmp) const
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
      
      END_CODE();

      return fin;
    }


    // Construct (A1xNabla_A1) source
    LatticePropagator
    MesA1xNablaA1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				 const LatticePropagator& tmp) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_5\gamma_i \nabla_i\f$  
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * rightNabla(tmp,u,k,params.deriv_length);
      
      END_CODE();

      return fin;
    }


    // Construct (A1xNabla_T2) source
    LatticePropagator
    MesA1xNablaT2SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				 const LatticePropagator& tmp) const
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
      
      END_CODE();

      return fin;
    }


    // Construct (A1xNabla_E) source
    LatticePropagator
    MesA1xNablaESrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				const LatticePropagator& tmp) const
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
      
      END_CODE();

      return fin;
    }


    // Construct (B1xNabla_T1) source
    LatticePropagator
    MesB1xNablaT1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				 const LatticePropagator& tmp) const
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
      
      END_CODE();

      return fin;
    }


    // Construct (A0_2xD_T2) source
    LatticePropagator
    MesA02xDT2SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
			      const LatticePropagator& tmp) const
    {
      START_CODE();

      LatticePropagator fin;

      // \f$\Gamma_f \equiv \gamma_4 D_i\f$  
      fin = Gamma(1 << 3) * rightD(tmp,u,params.deriv_dir,params.deriv_length);
      
      END_CODE();

      return fin;
    }


    // Construct (A1xD_A2) source
    LatticePropagator
    MesA1xDA2SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
			     const LatticePropagator& tmp) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_5\gamma_i D_i\f$  
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * rightD(tmp,u,k,params.deriv_length);
      
      END_CODE();

      return fin;
    }


    // Construct (A1xD_E) source
    LatticePropagator
    MesA1xDESrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
			    const LatticePropagator& tmp) const
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
      
      END_CODE();

      return fin;
    }


    // Construct (A1xD_T1) source
    LatticePropagator
    MesA1xDT1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
			     const LatticePropagator& tmp) const
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
      
      END_CODE();

      return fin;
    }

    
    // Construct (A1xD_T2) source
    LatticePropagator
    MesA1xDT2SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
			     const LatticePropagator& tmp) const
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
      
      END_CODE();

      return fin;
    }


    // Construct (B1xD_A2) source
    LatticePropagator
    MesB1xDA2SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
			     const LatticePropagator& tmp) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 \gamma_i D_i\f$  
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * rightD(tmp,u,k,params.deriv_length);

      END_CODE();

      return fin;
    }


    // Construct (B1xD_E) source
    LatticePropagator
    MesB1xDESrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
			    const LatticePropagator& tmp) const
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

      END_CODE();

      return fin;
    }


    // Construct (B1xD_T1) source
    LatticePropagator
    MesB1xDT1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
			     const LatticePropagator& tmp) const
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

      END_CODE();

      return fin;
    }


    // Construct (B1xD_T2) source
    LatticePropagator
    MesB1xDT2SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
			     const LatticePropagator& tmp) const
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

      END_CODE();

      return fin;
    }


    // Construct (RhoxD_A2) source
    LatticePropagator
    MesRhoxDA2SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
			      const LatticePropagator& tmp) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_i D_i\f$  
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * rightD(tmp,u,k,params.deriv_length);
      
      END_CODE();

      return fin;
    }


    //! Construct (RhoxD_T1) source
    LatticePropagator
    MesRhoxDT1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
			      const LatticePropagator& tmp) const
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
      
      END_CODE();

      return fin;
    }


    //! Construct (RhoxD_T2) source
    LatticePropagator
    MesRhoxDT2SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
			      const LatticePropagator& tmp) const
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
      
      END_CODE();

      return fin;
    }


    // Construct (PionxD_T2) source
    LatticePropagator
    MesPionxDT2SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
			       const LatticePropagator& tmp) const
    {
      START_CODE();

      LatticePropagator fin;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 D_i\f$  
      fin = Gamma(1 << 3) * (Gamma(G5) * rightD(tmp,u,params.deriv_dir,params.deriv_length));
      
      END_CODE();

      return fin;
    }

 
    //! Construct (PionxB_T1) source
    LatticePropagator
    MesPionxBT1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
			       const LatticePropagator& tmp) const
    {
      START_CODE();

      LatticePropagator fin;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_5 B_i\f$  
      fin = Gamma(G5) * rightB(tmp,u,params.deriv_dir,params.deriv_length);
      
      END_CODE();

      return fin;
    }


    //! Construct (RhoxB_T1) source
    LatticePropagator
    MesRhoxBT1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
			      const LatticePropagator& tmp) const
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
      
      END_CODE();

      return fin;
    }


    //! Construct (RhoxB_T2) source
    LatticePropagator
    MesRhoxBT2SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
			      const LatticePropagator& tmp) const
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
      
      END_CODE();

      return fin;
    }


    //! Construct (A1xB_A1) source
    LatticePropagator
    MesA1xBA1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
			     const LatticePropagator& tmp) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_5 \gamma_i B_i\f$  
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * rightB(tmp,u,k,params.deriv_length);
      
      END_CODE();

      return fin;
    }


    //! Construct (RhoxB_T1) source
    LatticePropagator
    MesA1xBT1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
			     const LatticePropagator& tmp) const
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
      
      END_CODE();

      return fin;
    }


    //! Construct (A1xB_T2) source
    LatticePropagator
    MesA1xBT2SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
			     const LatticePropagator& tmp) const
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
      
      END_CODE();

      return fin;
    }


    // Register all the possible deriv mesons
    bool registerAll(void) 
    {
      bool foo = true;

      foo &= LinkSmearingEnv::registered;
      foo &= QuarkSmearingEnv::registered;

      //! Register all the factories
      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-PIONxNABLA_T1"),
										 mesPionxNablaT1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-A0xNABLA_T1"),
										 mesA0xNablaT1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-A0_2xNABLA_T1"),
										 mesA02xNablaT1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-RHOxNABLA_A1"),
										 mesRhoxNablaA1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-RHOxNABLA_T1"),
										 mesRhoxNablaT1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-RHOxNABLA_T2"),
										 mesRhoxNablaT2SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-A1xNABLA_A1"),
										 mesA1xNablaA1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-A1xNABLA_T2"),
										 mesA1xNablaT2SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-A1xNABLA_E"),
										 mesA1xNablaESrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-B1xNABLA_T1"),
										 mesB1xNablaT1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-A0_2xD_T2"),
										 mesA02xDT2SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-A1xD_A2"),
										 mesA1xDA2SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-A1xD_E"),
										 mesA1xDESrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-A1xD_T1"),
										 mesA1xDT1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-A1xD_T2"),
										 mesA1xDT2SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-B1xD_A2"),
										 mesB1xDA2SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-B1xD_E"),
										 mesB1xDESrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-B1xD_T1"),
										 mesB1xDT1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-B1xD_T2"),
										 mesB1xDT2SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-RHOxD_A2"),
										 mesRhoxDA2SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-RHOxD_T1"),
										 mesRhoxDT1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-RHOxD_T2"),
										 mesRhoxDT2SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-PIONxD_T2"),
										 mesPionxDT2SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-PIONxB_T1"),
										 mesPionxBT1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-RHOxB_T1"),
										 mesRhoxBT1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-RHOxB_T2"),
										 mesRhoxBT2SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-A1xB_A1"),
										 mesA1xBA1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-A1xB_T1"),
										 mesA1xBT1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("SOURCE-A1xB_T2"),
										 mesA1xBT2SrcConst);

      return foo;
    }

    const bool registered = registerAll();

  }  // end namespace DerivShellSourceConstEnv

}  // end namespace Chroma



  
