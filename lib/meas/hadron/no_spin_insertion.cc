// $Id: no_spin_insertion.cc,v 1.2 2006-09-20 20:28:01 edwards Exp $
/*! \file
 *  \brief No spin insertion
 */

#include "chromabase.h"

#include "meas/hadron/spin_insertion_factory.h"
#include "meas/hadron/no_spin_insertion.h"

namespace Chroma 
{

  // Read parameters
  void read(XMLReader& xml, const string& path, NoSpinInsertionEnv::Params& param)
  {
    NoSpinInsertionEnv::Params tmp(xml, path);
    param = tmp;
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const string& path, const NoSpinInsertionEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Hooks to register the class
  namespace NoSpinInsertionEnv
  {
    namespace
    {
      //! Callback function
      SpinInsertion<LatticePropagator>* createProp(XMLReader& xml_in,
						   const std::string& path)
      {
	return new SpinInsert<LatticePropagator>(Params(xml_in, path));
      }

      //! Callback function
      SpinInsertion<LatticeFermion>* createFerm(XMLReader& xml_in,
						const std::string& path)
      {
	return new SpinInsert<LatticeFermion>(Params(xml_in, path));
      }
    
      //! Local registration flag
      bool registered = false;
    }

    //! Name to be used
    const std::string name = "NONE";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::ThePropSpinInsertionFactory::Instance().registerObject(name, createProp);
	success &= Chroma::TheFermSpinInsertionFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }


    //! Parameters for running code
    Params::Params(XMLReader& xml, const string& path)
    {
    }


    //! Parameters for running code
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);
    
      write(xml, "SpinInsertionType", name);

      pop(xml);
    }

  }  // end namespace
}  // end namespace Chroma

