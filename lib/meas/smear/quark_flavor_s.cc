// $Id: quark_flavor_s.cc,v 1.1 2006-11-20 22:13:28 kostas Exp $
/*! \file
 *  \brief Derivative displacements
 */

#include "chromabase.h"

#include "meas/smear/quark_flavor_s.h"
#include "meas/smear/quark_displacement_factory.h"
#include "meas/smear/quark_displacement_aggregate.h"

namespace Chroma 
{

  // Read parameters
  void read(XMLReader& xml, const string& path, StaggeredQuarkFlavorOpEnv::Params& param)
  {
    StaggeredQuarkFlavorOpEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const StaggeredQuarkFlavorOpEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Staggered Flavored sources
  /*! \ingroup sources */
  namespace StaggeredQuarkFlavorOpEnv
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


      //-------------------- callback functions ------------------

      //! Construct (right Nabla) source
      QuarkDisplacement<LatticeStaggeredPropagator>* scalar(XMLReader& xml_in,const std::string& path){
	return new StaggeredScalarOp<LatticeStaggeredPropagator>(Params(xml_in, path));
      }

      QuarkDisplacement<LatticeStaggeredPropagator>* vector(XMLReader& xml_in,const std::string& path){
	return new StaggeredVectorOp<LatticeStaggeredPropagator>(ParamsOneIndex(xml_in, path));
      }

      QuarkDisplacement<LatticeStaggeredPropagator>* axial_vector(XMLReader& xml_in,const std::string& path){
	return new StaggeredAxialVectorOp<LatticeStaggeredPropagator>(ParamsOneIndex(xml_in, path));
      }

      QuarkDisplacement<LatticeStaggeredPropagator>* tensor(XMLReader& xml_in,const std::string& path){
	return new StaggeredTensorOp<LatticeStaggeredPropagator>(ParamsTwoIndex(xml_in, path));
      }

      QuarkDisplacement<LatticeStaggeredPropagator>* pseudo_scalar(XMLReader& xml_in,const std::string& path){
	return new StaggeredPseudoScalarOp<LatticeStaggeredPropagator>(Params(xml_in, path));
      }


    } // end anonymous namespace


    //! Initialize
    Params::Params()
    {    }

    //! Read parameters
    Params::Params(XMLReader& xml, const string& path)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
	read(paramtop, "FlavorOp",  FlavorOp);
	break;

      default:
	QDPIO::cerr << __func__ << ": parameter version " << version 
		    << " unsupported." << endl;
	QDP_abort(1);
      }

    }

    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      QDP::write(xml, "version", version);

      write(xml,  "FlavorOp",  FlavorOp);
     
      pop(xml);
    }


//! Read parameters
    ParamsOneIndex::ParamsOneIndex(XMLReader& xml, const string& path)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
	read(paramtop, "FlavorOp",  FlavorOp);
	read(paramtop, "mu",  mu);

	break;

      default:
	QDPIO::cerr << __func__ << ": parameter version " << version 
		    << " unsupported." << endl;
	QDP_abort(1);
      }

    }

    // Writer
    void ParamsOneIndex::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      QDP::write(xml, "version", version);

      write(xml,  "FlavorOp",  FlavorOp);
      write(xml,  "mu",  mu);

      pop(xml);
    }


//! Read parameters
    ParamsTwoIndex::ParamsTwoIndex(XMLReader& xml, const string& path)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
	read(paramtop, "FlavorOp",  FlavorOp);
	read(paramtop, "mu",  mu);
	read(paramtop, "nu",  mu);

	break;

      default:
	QDPIO::cerr << __func__ << ": parameter version " << version 
		    << " unsupported." << endl;
	QDP_abort(1);
      }

    }

    // Writer
    void ParamsTwoIndex::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      QDP::write(xml, "version", version);

      write(xml,  "FlavorOp",  FlavorOp);
      write(xml,  "mu",  mu);
      write(xml,  "nu",  nu);

      pop(xml);
    }

    // Construct scalar source
    // See corresponding .h file for doxygen comments
    template<>
    void
    StaggeredScalarOp<LatticeStaggeredPropagator>::operator()(LatticeStaggeredPropagator& tmp,
								   const multi1d<LatticeColorMatrix>& u,
								   enum PlusMinus isign) const
    {
      START_CODE();

      // do nothing
      // good job! there are no bugs in this code!!

      END_CODE();
    }


    // Construct pseudoscalar source
    // See corresponding .h file for doxygen comments
    template<>
    void
    StaggeredPseudoScalarOp<LatticeStaggeredPropagator>::operator()(LatticeStaggeredPropagator& tmp,
								   const multi1d<LatticeColorMatrix>& u,
								   enum PlusMinus isign) const
    {
      START_CODE();

      FlavorPseudoScalar(tmp, tmp, u) ;

      END_CODE();
    }

    // Construct vector source
    // See corresponding .h file for doxygen comments
    template<>
    void
    StaggeredVectorOp<LatticeStaggeredPropagator>::operator()(LatticeStaggeredPropagator& tmp,
							      const multi1d<LatticeColorMatrix>& u,
							      enum PlusMinus isign) const
    {
      START_CODE();

      FlavorVector(tmp, tmp, u, params.mu) ;

      END_CODE();
    }

    // Construct axial vector source
    // See corresponding .h file for doxygen comments
    template<>
    void
    StaggeredAxialVectorOp<LatticeStaggeredPropagator>::operator()(LatticeStaggeredPropagator& tmp,
							      const multi1d<LatticeColorMatrix>& u,
							      enum PlusMinus isign) const
    {
      START_CODE();

      FlavorAxialVector(tmp, tmp, u, params.mu) ;
      
      END_CODE();
    }


    // Construct vector source
    // See corresponding .h file for doxygen comments
    template<>
    void
    StaggeredTensorOp<LatticeStaggeredPropagator>::operator()(LatticeStaggeredPropagator& tmp,
							      const multi1d<LatticeColorMatrix>& u,
							      enum PlusMinus isign) const
    {
      START_CODE();

      FlavorTensor(tmp, tmp, u, params.mu, params.nu) ;

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
	success &= Chroma::TheStagPropDisplacementFactory::Instance().registerObject(string("SCALAR"), scalar);
	success &= Chroma::TheStagPropDisplacementFactory::Instance().registerObject(string("VECTOR"), vector);
	success &= Chroma::TheStagPropDisplacementFactory::Instance().registerObject(string("AXIAL_VECTOR"), axial_vector);
	success &= Chroma::TheStagPropDisplacementFactory::Instance().registerObject(string("TENSOR"), tensor);
	success &= Chroma::TheStagPropDisplacementFactory::Instance().registerObject(string("PSEUDO_SCALAR"), pseudo_scalar);

	registered = true;
      }
      return success;
    }
  }  // end namespace StaggeredQuarkFlavorOpEnv

}  // end namespace Chroma



  
