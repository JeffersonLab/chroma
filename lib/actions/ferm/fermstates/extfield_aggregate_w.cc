// $Id: extfield_aggregate_w.cc,v 1.1 2006-09-19 17:53:36 edwards Exp $
/*! \file
 *  \brief External field aggregate
 */

#include "chromabase.h"

#include "actions/ferm/fermstates/extfield_factory_w.h"
#include "actions/ferm/fermstates/extfield_aggregate_w.h"

namespace Chroma 
{

#if 0
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
#endif


  //! External fields
  /*! \ingroup fermstates */
  namespace ExternalFieldEnv
  { 
    //! Anonymous namespace
    namespace
    {

      //-------------------- callback functions ---------------------------------------

      //! Construct zero
      ExternalField* zeroFunc(XMLReader& xml_in,
			      const std::string& path)
      {
	return new ZeroExternalField();
      }

#if 0
      //! Construct linear term
      ExternalField* linearFunc(XMLReader& xml_in,
				const std::string& path)
      {
	return new LinearExternalField(LinearParams(xml_in, path));
      }
#endif
      
    } // end anonymous namespace


#if 0
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
#endif


    // Construct linear function
    LatticeComplex
    ZeroExternalField::operator()(int dir) const
    {
      START_CODE();

      LatticeComplex d = 1;

      END_CODE();
      return d;
    }



    // Register all the possible external field funcs
    bool registerAll(void) 
    {
      bool foo = true;

      //! Register all the factories
      foo &= Chroma::TheExternalFieldFactory::Instance().registerObject(string("ZERO"),
									zeroFunc);
      fprintf(stderr,"registered ZERO\n");
      
//      //! Register all the factories
//      foo &= Chroma::TheExternalFieldFactory::Instance().registerObject(string("LINEAR"),
//									linearFunc);
      
      return foo;
    }

    const bool registered = registerAll();

  }  // end namespace ExternalFieldEnv

}  // end namespace Chroma



  
