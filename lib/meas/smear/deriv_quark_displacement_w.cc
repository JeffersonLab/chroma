// $Id: deriv_quark_displacement_w.cc,v 3.5 2008-11-11 21:27:42 edwards Exp $
/*! \file
 *  \brief Derivative displacements
 */

#include "chromabase.h"

#include "meas/smear/deriv_quark_displacement_w.h"
#include "meas/smear/quark_displacement_factory.h"
#include "meas/smear/quark_displacement_aggregate.h"
#include "meas/smear/displace.h"

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




      //-------------------- callback functions ---------------------------------------

      //! Construct (right Nabla) source
      QuarkDisplacement<LatticePropagator>* rightNablaDisplace(XMLReader& xml_in,
							       const std::string& path)
      {
	return new RightNablaDisplace<LatticePropagator>(ParamsDir(xml_in, path));
      }

      //! Construct (right D) source
      QuarkDisplacement<LatticePropagator>* rightDDisplace(XMLReader& xml_in,
							   const std::string& path)
      {
	return new RightDDisplace<LatticePropagator>(ParamsDir(xml_in, path));
      }

      //! Construct (right B) source
      QuarkDisplacement<LatticePropagator>* rightBDisplace(XMLReader& xml_in,
							   const std::string& path)
      {
	return new RightBDisplace<LatticePropagator>(ParamsDir(xml_in, path));
      }

      //! Construct (right E) source
      QuarkDisplacement<LatticePropagator>* rightEDisplace(XMLReader& xml_in,
							   const std::string& path)
      {
	return new RightEDisplace<LatticePropagator>(ParamsDir(xml_in, path));
      }

      //! Construct (right Laplacian) source
      QuarkDisplacement<LatticePropagator>* rightLapDisplace(XMLReader& xml_in,
							     const std::string& path)
      {
	return new RightLapDisplace<LatticePropagator>(Params(xml_in, path));
      }



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


    // Construct (right Nabla) source
    // See corresponding .h file for doxygen comments
    template<>
    void
    RightNablaDisplace<LatticePropagator>::operator()(LatticePropagator& tmp,
						      const multi1d<LatticeColorMatrix>& u,
						      enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin;
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \nabla_i\f$
      fin = rightNabla(tmp,u,params.deriv_dir,length);
      tmp = fin;

      END_CODE();
    }


    // Construct (right D) source
    // See corresponding .h file for doxygen comments
    template<>
    void
    RightDDisplace<LatticePropagator>::operator()(LatticePropagator& tmp,
						  const multi1d<LatticeColorMatrix>& u,
						  enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin;
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv D_i\f$
      fin = rightD(tmp,u,params.deriv_dir,length);
      tmp = fin;

      END_CODE();
    }


    // Construct (right B) source
    // See corresponding .h file for doxygen comments
    template<>
    void
    RightBDisplace<LatticePropagator>::operator()(LatticePropagator& tmp,
						  const multi1d<LatticeColorMatrix>& u,
						  enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin;
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv B_i\f$
      fin = rightB(tmp,u,params.deriv_dir,length);
      tmp = fin;

      END_CODE();
    }



    // Construct (right E) source
    // See corresponding .h file for doxygen comments
    template<>
    void
    RightEDisplace<LatticePropagator>::operator()(LatticePropagator& tmp,
						  const multi1d<LatticeColorMatrix>& u,
						  enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin;
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv E_alpha\f$
      fin = rightE(tmp,u,params.deriv_dir,length);
      tmp = fin;

      END_CODE();
    }



    // Construct (right Laplacian) source
    // See corresponding .h file for doxygen comments
    template<>
    void
    RightLapDisplace<LatticePropagator>::operator()(LatticePropagator& tmp,
						    const multi1d<LatticeColorMatrix>& u,
						    enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin;
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv Laplacian\f$
      fin = rightLap(tmp,u,length);
      tmp = fin;

      END_CODE();
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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \gamma_5\nabla_i\f$
      fin = Gamma(G5) * rightNabla(tmp,u,params.deriv_dir,length);
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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \nabla_i\f$
      fin = rightNabla(tmp,u,params.deriv_dir,length);
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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \gamma_4 \nabla_i\f$
      fin = Gamma(1 << 3) * rightNabla(tmp,u,params.deriv_dir,length);
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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \gamma_i\nabla_i\f$
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * rightNabla(tmp,u,k,length);
      
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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (antiSymTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(antiSymTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightNabla(tmp,u,k,length));
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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv s_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (symTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(symTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightNabla(tmp,u,k,length));
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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \gamma_5\gamma_i \nabla_i\f$  
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * rightNabla(tmp,u,k,length);
      
      tmp = Gamma(G5) * fin;

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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (symTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(symTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightNabla(tmp,u,k,length));
	}
      
      tmp = Gamma(G5) * fin;

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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \gamma_5 S_{\alpha jk}\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  Real e = ETensor3d(params.deriv_dir,j,k);
	  if (toBool(e != 0.0))
	    fin += e * (Gamma(1 << j) * rightNabla(tmp,u,k,length));
	}
      
      tmp = Gamma(G5) * fin;

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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \gamma_4\gamma_5\epsilon_{ijk}\gamma_j \nabla_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (antiSymTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(antiSymTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightD(tmp,u,k,length));
	}
      
      tmp = Gamma(1 << 3) * (Gamma(G5) * fin);

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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \gamma_4 D_i\f$  
      fin = Gamma(1 << 3) * rightD(tmp,u,params.deriv_dir,length);
      
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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \gamma_5\gamma_i D_i\f$  
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * rightD(tmp,u,k,length);
      
      tmp = Gamma(G5) * fin;

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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \gamma_5 S_{\alpha jk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  Real e = ETensor3d(params.deriv_dir,j,k);
	  if (toBool(e != 0.0))
	    fin += e * (Gamma(1 << j) * rightD(tmp,u,k,length));
	}
      
      tmp = Gamma(G5) * fin;

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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (symTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(symTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightD(tmp,u,k,length));
	}
      
      tmp = Gamma(G5) * fin;

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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \gamma_5\epsilon_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (antiSymTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(antiSymTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightD(tmp,u,k,length));
	}
      
      tmp = Gamma(G5) * fin;

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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 \gamma_i D_i\f$  
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * rightD(tmp,u,k,length);

      tmp = Gamma(1 << 3) * (Gamma(G5) * fin);

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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 S_{\alpha jk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  Real e = ETensor3d(params.deriv_dir,j,k);
	  if (toBool(e != 0.0))
	    fin += e * (Gamma(1 << j) * rightD(tmp,u,k,length));
	}

      tmp = Gamma(1 << 3) * (Gamma(G5) * fin);

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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 s_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (symTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(symTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightD(tmp,u,k,length));
	}

      tmp = Gamma(1 << 3) * (Gamma(G5) * fin);

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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 \epsilon_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (antiSymTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(antiSymTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightD(tmp,u,k,length));
	}

      tmp = Gamma(1 << 3) * (Gamma(G5) * fin);

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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \gamma_i D_i\f$  
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * rightD(tmp,u,k,length);
      
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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv s_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (symTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(symTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightD(tmp,u,k,length));
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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j D_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (antiSymTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(antiSymTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightD(tmp,u,k,length));
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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \gamma_4\gamma_5 D_i\f$  
      fin = Gamma(1 << 3) * (Gamma(G5) * rightD(tmp,u,params.deriv_dir,length));
      
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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \gamma_5 B_i\f$  
      fin = Gamma(G5) * rightB(tmp,u,params.deriv_dir,length);
      
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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \epsilon_{ijk}\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (antiSymTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(antiSymTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightB(tmp,u,k,length));
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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv s_{ijk}\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (symTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(symTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightB(tmp,u,k,length));
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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \gamma_5 \gamma_i B_i\f$  
      for(int k=0; k < 3; ++k)
	fin += Gamma(1 << k) * rightB(tmp,u,k,length);
      
      tmp = Gamma(G5) * fin;

      END_CODE();
    }


    //! Construct (A1xB_T1) source
    template<>
    void
    MesA1xBT1Displace<LatticePropagator>::operator()(LatticePropagator& tmp,
						     const multi1d<LatticeColorMatrix>& u,
						     enum PlusMinus isign) const
    {
      START_CODE();

      LatticePropagator fin = zero;
      int G5 = Ns*Ns-1;
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \gamma_5 \epsilon_{ijk}\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (antiSymTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(antiSymTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightB(tmp,u,k,length));
	}
      
      tmp = Gamma(G5) * fin;

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
      int length = plusMinus(isign) * params.deriv_length;

      // \f$\Gamma_f \equiv \gamma_5 s_{ijk}\gamma_j B_k\f$  
      for(int j=0; j < 3; ++j)
	for(int k=0; k < 3; ++k)
	{
	  if (symTensor3d(params.deriv_dir,j,k) != 0)
	    fin += Real(symTensor3d(params.deriv_dir,j,k)) * (Gamma(1 << j) * rightB(tmp,u,k,length));
	}
      
      tmp = Gamma(G5) * fin;

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
	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("NABLA-DERIV"),
									     rightNablaDisplace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("D-DERIV"),
									     rightDDisplace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("B-DERIV"),
									     rightBDisplace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("E-DERIV"),
									     rightEDisplace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("LAP-DERIV"),
									     rightLapDisplace);


	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("PIONxNABLA_T1-DERIV"),
									     mesPionxNablaT1Displace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("A0xNABLA_T1-DERIV"),
									     mesA0xNablaT1Displace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("A0_2xNABLA_T1-DERIV"),
									     mesA02xNablaT1Displace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("RHOxNABLA_A1-DERIV"),
									     mesRhoxNablaA1Displace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("RHOxNABLA_T1-DERIV"),
									     mesRhoxNablaT1Displace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("RHOxNABLA_T2-DERIV"),
									     mesRhoxNablaT2Displace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("A1xNABLA_A1-DERIV"),
									     mesA1xNablaA1Displace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("A1xNABLA_T2-DERIV"),
									     mesA1xNablaT2Displace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("A1xNABLA_E-DERIV"),
									     mesA1xNablaEDisplace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("B1xNABLA_T1-DERIV"),
									     mesB1xNablaT1Displace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("A0_2xD_T2-DERIV"),
									     mesA02xDT2Displace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("A1xD_A2-DERIV"),
									     mesA1xDA2Displace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("A1xD_E-DERIV"),
									     mesA1xDEDisplace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("A1xD_T1-DERIV"),
									     mesA1xDT1Displace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("A1xD_T2-DERIV"),
									     mesA1xDT2Displace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("B1xD_A2-DERIV"),
									     mesB1xDA2Displace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("B1xD_E-DERIV"),
									     mesB1xDEDisplace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("B1xD_T1-DERIV"),
									     mesB1xDT1Displace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("B1xD_T2-DERIV"),
									     mesB1xDT2Displace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("RHOxD_A2-DERIV"),
									     mesRhoxDA2Displace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("RHOxD_T1-DERIV"),
									     mesRhoxDT1Displace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("RHOxD_T2-DERIV"),
									     mesRhoxDT2Displace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("PIONxD_T2-DERIV"),
									     mesPionxDT2Displace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("PIONxB_T1-DERIV"),
									     mesPionxBT1Displace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("RHOxB_T1-DERIV"),
									     mesRhoxBT1Displace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("RHOxB_T2-DERIV"),
									     mesRhoxBT2Displace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("A1xB_A1-DERIV"),
									     mesA1xBA1Displace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("A1xB_T1-DERIV"),
									     mesA1xBT1Displace);

	success &= Chroma::ThePropDisplacementFactory::Instance().registerObject(string("A1xB_T2-DERIV"),
									     mesA1xBT2Displace);

	registered = true;
      }
      return success;
    }
  }  // end namespace DerivQuarkDisplacementEnv

}  // end namespace Chroma



  
