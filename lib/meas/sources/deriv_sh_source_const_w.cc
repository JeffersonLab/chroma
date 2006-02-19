// $Id: deriv_sh_source_const_w.cc,v 1.1 2006-02-19 05:22:08 edwards Exp $
/*! \file
 *  \brief Construct derivative source construction
 */

#include "meas/hadron/deriv_sh_source_const_w.h"
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
  /*! \ingroup hadron */
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
       * \ingroup hadron
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
       * \ingroup hadron
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
       * \ingroup hadron
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

      //! Construct pion_1-(PionxNabla_T1) source
      QuarkSourceConstruction<LatticePropagator>* mesPionPionxNablaT1SrcConst(XMLReader& xml_in,
									      const std::string& path)
      {
	return new MesPionPionxNablaT1SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(A0xNabla_T1) source
      QuarkSourceConstruction<LatticePropagator>* mesPionA0xNablaT1SrcConst(XMLReader& xml_in,
									    const std::string& path)
      {
	return new MesPionA0xNablaT1SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(A0_2xNabla_T1) source
      QuarkSourceConstruction<LatticePropagator>* mesPionA02xNablaT1SrcConst(XMLReader& xml_in,
									     const std::string& path)
      {
	return new MesPionA02xNablaT1SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(RhoxNabla_A1) source
      QuarkSourceConstruction<LatticePropagator>* mesPionRhoxNablaA1SrcConst(XMLReader& xml_in,
									     const std::string& path)
      {
	return new MesPionRhoxNablaA1SrcConst(Params(xml_in, path));
      }

      //! Construct pion_1-(RhoxNabla_T1) source
      QuarkSourceConstruction<LatticePropagator>* mesPionRhoxNablaT1SrcConst(XMLReader& xml_in,
									     const std::string& path)
      {
	return new MesPionRhoxNablaT1SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(RhoxNabla_T2) source
      QuarkSourceConstruction<LatticePropagator>* mesPionRhoxNablaT2SrcConst(XMLReader& xml_in,
									     const std::string& path)
      {
	return new MesPionRhoxNablaT2SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(A1xNabla_A1) source
      QuarkSourceConstruction<LatticePropagator>* mesPionA1xNablaA1SrcConst(XMLReader& xml_in,
									    const std::string& path)
      {
	return new MesPionA1xNablaA1SrcConst(Params(xml_in, path));
      }

      //! Construct pion_1-(A1xNabla_T2) source
      QuarkSourceConstruction<LatticePropagator>* mesPionA1xNablaT2SrcConst(XMLReader& xml_in,
									    const std::string& path)
      {
	return new MesPionA1xNablaT2SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(A1xNabla_E) source
      QuarkSourceConstruction<LatticePropagator>* mesPionA1xNablaESrcConst(XMLReader& xml_in,
									   const std::string& path)
      {
	return new MesPionA1xNablaESrcConst(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(B1xNabla_T1) source
      QuarkSourceConstruction<LatticePropagator>* mesPionB1xNablaT1SrcConst(XMLReader& xml_in,
									    const std::string& path)
      {
	return new MesPionB1xNablaT1SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(A0_2xD_T2) source
      QuarkSourceConstruction<LatticePropagator>* mesPionA02xDT2SrcConst(XMLReader& xml_in,
									 const std::string& path)
      {
	return new MesPionA02xDT2SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(A1xD_A2) source
      QuarkSourceConstruction<LatticePropagator>* mesPionA1xDA2SrcConst(XMLReader& xml_in,
									const std::string& path)
      {
	return new MesPionA1xDA2SrcConst(Params(xml_in, path));
      }

      //! Construct pion_1-(A1xD_E) source
      QuarkSourceConstruction<LatticePropagator>* mesPionA1xDESrcConst(XMLReader& xml_in,
								       const std::string& path)
      {
	return new MesPionA1xDESrcConst(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(A1xD_T1) source
      QuarkSourceConstruction<LatticePropagator>* mesPionA1xDT1SrcConst(XMLReader& xml_in,
									const std::string& path)
      {
	return new MesPionA1xDT1SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(A1xD_T2) source
      QuarkSourceConstruction<LatticePropagator>* mesPionA1xDT2SrcConst(XMLReader& xml_in,
									const std::string& path)
      {
	return new MesPionA1xDT2SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(B1xD_A2) source
      QuarkSourceConstruction<LatticePropagator>* mesPionB1xDA2SrcConst(XMLReader& xml_in,
									const std::string& path)
      {
	return new MesPionB1xDA2SrcConst(Params(xml_in, path));
      }

      //! Construct pion_1-(B1xD_E) source
      QuarkSourceConstruction<LatticePropagator>* mesPionB1xDESrcConst(XMLReader& xml_in,
								       const std::string& path)
      {
	return new MesPionB1xDESrcConst(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(B1xD_T1) source
      QuarkSourceConstruction<LatticePropagator>* mesPionB1xDT1SrcConst(XMLReader& xml_in,
									const std::string& path)
      {
	return new MesPionB1xDT1SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(B1xD_T2) source
      QuarkSourceConstruction<LatticePropagator>* mesPionB1xDT2SrcConst(XMLReader& xml_in,
									const std::string& path)
      {
	return new MesPionB1xDT2SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(RhoxD_A2) source
      QuarkSourceConstruction<LatticePropagator>* mesPionRhoxDA2SrcConst(XMLReader& xml_in,
									 const std::string& path)
      {
	return new MesPionRhoxDA2SrcConst(Params(xml_in, path));
      }

      //! Construct pion_1-(RhoxD_T1) source
      QuarkSourceConstruction<LatticePropagator>* mesPionRhoxDT1SrcConst(XMLReader& xml_in,
									 const std::string& path)
      {
	return new MesPionRhoxDT1SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(RhoxD_T2) source
      QuarkSourceConstruction<LatticePropagator>* mesPionRhoxDT2SrcConst(XMLReader& xml_in,
									 const std::string& path)
      {
	return new MesPionRhoxDT2SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(PionxD_T2) source
      QuarkSourceConstruction<LatticePropagator>* mesPionPionxDT2SrcConst(XMLReader& xml_in,
									  const std::string& path)
      {
	return new MesPionPionxDT2SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(PionxB_T1) source
      QuarkSourceConstruction<LatticePropagator>* mesPionPionxBT1SrcConst(XMLReader& xml_in,
									  const std::string& path)
      {
	return new MesPionPionxBT1SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(RhoxB_T1) source
      QuarkSourceConstruction<LatticePropagator>* mesPionRhoxBT1SrcConst(XMLReader& xml_in,
									 const std::string& path)
      {
	return new MesPionRhoxBT1SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(RhoxB_T2) source
      QuarkSourceConstruction<LatticePropagator>* mesPionRhoxBT2SrcConst(XMLReader& xml_in,
									 const std::string& path)
      {
	return new MesPionRhoxBT2SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(A1xB_A1) source
      QuarkSourceConstruction<LatticePropagator>* mesPionA1xBA1SrcConst(XMLReader& xml_in,
									const std::string& path)
      {
	return new MesPionA1xBA1SrcConst(Params(xml_in, path));
      }

      //! Construct pion_1-(A1xB_T1) source
      QuarkSourceConstruction<LatticePropagator>* mesPionA1xBT1SrcConst(XMLReader& xml_in,
									const std::string& path)
      {
	return new MesPionA1xBT1SrcConst(ParamsDir(xml_in, path));
      }

      //! Construct pion_1-(A1xB_T2) source
      QuarkSourceConstruction<LatticePropagator>* mesPionA1xBT2SrcConst(XMLReader& xml_in,
									const std::string& path)
      {
	return new MesPionA1xBT2SrcConst(ParamsDir(xml_in, path));
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
    LatticePropagator
    DerivSourceConst::operator()(const multi1d<LatticeColorMatrix>& u) const
    {
      QDPIO::cout << "Deriv Shell source" << endl;

      LatticePropagator quark_source;
      const DerivShellSourceConstEnv::Params& = getParams();

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
	QDPIO::cerr << name << ": Caught Exception in deriv source construction: " << e << endl;
	QDP_abort(1);
      }


      return quark_source;
    }




    // Construct pion_1-(PionxNabla_T1) source
    // See corresponding .h file for doxygen comments
    LatticePropagator
    MesPionPionxNablaT1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				       const multi1d<LatticePropagator>& tmp) const
    {
      START_CODE();

      LatticePropagator fin;
      const int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \nabla_i\f$
      fin = Gamma(G5) * rightNabla(tmp,u,params.deriv_dir,params.deriv_length);
      
      END_CODE();

      return fin;
    }


    // Construct pion_1-(A0xNabla_T1) source
    // See corresponding .h file for doxygen comments
    LatticePropagator
    MesPionA0xNablaT1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				     const multi1d<LatticePropagator>& tmp) const
    {
      START_CODE();

      LatticePropagator fin;

      // \f$\Gamma_f \equiv \nabla_i\f$
      fin = rightNabla(tmp,u,params.deriv_dir,params.deriv_length);
      
      END_CODE();

      return fin;
    }


    // Construct pion_1-(A0_2xNabla_T1) source
    // See corresponding .h file for doxygen comments
    LatticePropagator
    MesPionA02xNablaT1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<LatticePropagator>& tmp) const
    {
      START_CODE();

      LatticePropagator fin;

      // \f$\Gamma_f \equiv \gamma_4 \nabla_i\f$
      fin = Gamma(1 << 3) * rightNabla(tmp,u,params.deriv_dir,params.deriv_length);
      
      END_CODE();
    }


    // Construct pion_1-(RhoxNabla_A1) source
    LatticePropagator
    MesPionRhoxNablaA1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<LatticePropagator>& tmp) const
    {
      START_CODE();

      LatticePropagator fin = zero;

      // \f$\Gamma_f \equiv \gamma_i\nabla_i\f$
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * rightNabla(tmp,u,k,params.deriv_length);
      
      END_CODE();

      return fin;
    }


    // Construct pion_1-(RhoxNabla_T1) source
    LatticePropagator
    MesPionRhoxNablaT1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<LatticePropagator>& tmp) const
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


    // Construct pion_1-(RhoxNabla_T2) source
    LatticePropagator
    MesPionRhoxNablaT2SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				      const multi1d<LatticePropagator>& tmp) const
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


    // Construct pion_1-(A1xNabla_A1) source
    LatticePropagator
    MesPionA1xNablaA1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				     const multi1d<LatticePropagator>& tmp) const
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


    // Construct pion_1-(A1xNabla_T2) source
    LatticePropagator
    MesPionA1xNablaT2SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				     const multi1d<LatticePropagator>& tmp) const
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


    // Construct pion_1-(A1xNabla_E) source
    LatticePropagator
    MesPionA1xNablaESrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				    const multi1d<LatticePropagator>& tmp) const
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


    // Construct pion_1-(B1xNabla_T1) source
    LatticePropagator
    MesPionB1xNablaT1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				     const multi1d<LatticePropagator>& tmp) const
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


    // Construct pion_1-(A0_2xD_T2) source
    LatticePropagator
    MesPionA02xDT2SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<LatticePropagator>& tmp) const
    {
      START_CODE();

      LatticePropagator fin;

      // \f$\Gamma_f \equiv \gamma_4 D_i\f$  
      fin = Gamma(1 << 3) * rightD(tmp,u,params.deriv_dir,params.deriv_length);
      
      END_CODE();

      return fin;
    }


    // Construct pion_1-(A1xD_A2) source
    LatticePropagator
    MesPionA1xDA2SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<LatticePropagator>& tmp) const
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


    // Construct pion_1-(A1xD_E) source
    LatticePropagator
    MesPionA1xDESrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				const multi1d<LatticePropagator>& tmp) const
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


    // Construct pion_1-(A1xD_T1) source
    LatticePropagator
    MesPionA1xDT1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<LatticePropagator>& tmp) const
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

    
    // Construct pion_1-(A1xD_T2) source
    LatticePropagator
    MesPionA1xDT2SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<LatticePropagator>& tmp) const
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


    // Construct pion_1-(B1xD_A2) source
    LatticePropagator
    MesPionB1xDA2SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				       const multi1d<LatticePropagator>& tmp) const
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


    // Construct pion_1-(B1xD_E) source
    LatticePropagator
    MesPionB1xDESrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				const multi1d<LatticePropagator>& tmp) const
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


    // Construct pion_1-(B1xD_T1) source
    LatticePropagator
    MesPionB1xDT1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<LatticePropagator>& tmp) const
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


    // Construct pion_1-(B1xD_T2) source
    LatticePropagator
    MesPionB1xDT2SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<LatticePropagator>& tmp) const
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


    // Construct pion_1-(RhoxD_A2) source
    LatticePropagator
    MesPionRhoxDA2SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<LatticePropagator>& tmp) const
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


    //! Construct pion_1-(RhoxD_T1) source
    LatticePropagator
    MesPionRhoxDT1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<LatticePropagator>& tmp) const
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


    //! Construct pion_1-(RhoxD_T2) source
    LatticePropagator
    MesPionRhoxDT2SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<LatticePropagator>& tmp) const
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


    // Construct pion_1-(PionxD_T2) source
    LatticePropagator
    MesPionPionxDT2SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& tmp) const
    {
      START_CODE();

      LatticePropagator fin;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 D_i\f$  
      fin = Gamma(1 << 3) * (Gamma(G5) * rightD(tmp,u,params.deriv_dir,params.deriv_length));
      
      END_CODE();

      return fin;
    }

 
    //! Construct pion_1-(PionxB_T1) source
    LatticePropagator
    MesPionPionxBT1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				   const multi1d<LatticePropagator>& tmp) const
    {
      START_CODE();

      LatticePropagator fin;
      int G5 = Ns*Ns-1;

      // \f$\Gamma_f \equiv \gamma_5 B_i\f$  
      fin = Gamma(G5) * rightB(tmp,u,params.deriv_dir,params.deriv_length);
      
      END_CODE();

      return fin;
    }


    //! Construct pion_1-(RhoxB_T1) source
    LatticePropagator
    MesPionRhoxBT1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<LatticePropagator>& tmp) const
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


    //! Construct pion_1-(RhoxB_T2) source
    LatticePropagator
    MesPionRhoxBT2SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				  const multi1d<LatticePropagator>& tmp) const
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


    //! Construct pion_1-(A1xB_A1) source
    LatticePropagator
    MesPionA1xBA1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<LatticePropagator>& tmp) const
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


    //! Construct pion_1-(RhoxB_T1) source
    LatticePropagator
    MesPionA1xBT1SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<LatticePropagator>& tmp) const
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


    //! Construct pion_1-(A1xB_T2) source
    LatticePropagator
    MesPionA1xBT2SrcConst::deriv(const multi1d<LatticeColorMatrix>& u,
				 const multi1d<LatticePropagator>& tmp) const
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
      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-PIONxNABLA_T2"),
										 mesPionA1xDT2SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-A0xNABLA_T1"),
										 mesPionA0xNablaT1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-A0_2xNABLA_T1"),
										 mesPionA02xNablaT1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-RHOxNABLA_A1"),
										 mesPionRhoxNablaA1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-RHOxNABLA_T1"),
										 mesPionRhoxNablaT1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-RHOxNABLA_T2"),
										 mesPionRhoxNablaT2SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-A1xNABLA_A1"),
										 mesPionA1xNablaA1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-A1xNABLA_T2"),
										 mesPionA1xNablaT2SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-A1xNABLA_E"),
										 mesPionA1xNablaESrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-B1xNABLA_T1"),
										 mesPionB1xNablaT1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-A0_2xD_T2"),
										 mesPionA02xDT2SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-A1xD_A2"),
										 mesPionA1xDA2SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-A1xD_E"),
										 mesPionA1xDESrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-A1xD_T1"),
										 mesPionA1xDT1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-A1xD_T2"),
										 mesPionA1xDT2SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-B1xD_A2"),
										 mesPionB1xDA2SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-B1xD_E"),
										 mesPionB1xDESrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-B1xD_T1"),
										 mesPionB1xDT1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-B1xD_T2"),
										 mesPionB1xDT2SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-RHOxD_A2"),
										 mesPionRhoxDA2SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-RHOxD_T1"),
										 mesPionRhoxDT1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-RHOxD_T2"),
										 mesPionRhoxDT2SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-PIONxD_T2"),
										 mesPionPionxDT2SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-PIONxB_T1"),
										 mesPionPionxBT1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-RHOxB_T1"),
										 mesPionRhoxBT1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-RHOxB_T2"),
										 mesPionRhoxBT2SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-A1xB_A1"),
										 mesPionA1xBA1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-A1xB_T1"),
										 mesPionA1xBT1SrcConst);

      foo &= Chroma::ThePropSourceConstructionFactory::Instance().registerObject(string("PION-A1xB_T2"),
										 mesPionA1xBT2SrcConst);

      return foo;
    }

    const bool registered = registerAll();

  }  // end namespace DerivShellSourceConstEnv

}  // end namespace Chroma



  
