// $Id: simple_spin_insertion_w.cc,v 1.1 2006-05-24 21:09:41 edwards Exp $
/*! \file
 *  \brief Gamma insertion
 */

#include "chromabase.h"

#include "meas/hadron/spin_insertion_factory.h"
#include "meas/hadron/simple_spin_insertion_w.h"

namespace Chroma 
{

  // Read parameters
  void read(XMLReader& xml, const string& path, SimpleSpinInsertionEnv::Params& param)
  {
    SimpleSpinInsertionEnv::Params tmp(xml, path);
    param = tmp;
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const string& path, const SimpleSpinInsertionEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Hooks to register the class
  namespace SimpleSpinInsertionEnv
  {
    //! Callback function
    SpinInsertion<LatticePropagator>* leftSpinProp(XMLReader& xml_in,
						   const std::string& path)
    {
      return new LeftSpinInsert<LatticePropagator>(Params(xml_in, path));
    }

    //! Callback function
    SpinInsertion<LatticeFermion>* leftSpinFerm(XMLReader& xml_in,
						const std::string& path)
    {
      return new LeftSpinInsert<LatticeFermion>(Params(xml_in, path));
    }
    
    //! Callback function
    SpinInsertion<LatticePropagator>* rightSpinProp(XMLReader& xml_in,
						    const std::string& path)
    {
      return new RightSpinInsert<LatticePropagator>(Params(xml_in, path));
    }

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= Chroma::ThePropSpinInsertionFactory::Instance().registerObject("LEFT_GAMMA_INSERTION", 
									    leftSpinProp);
      foo &= Chroma::ThePropSpinInsertionFactory::Instance().registerObject("RIGHT_GAMMA_INSERTION", 
									    rightSpinProp);
      foo &= Chroma::TheFermSpinInsertionFactory::Instance().registerObject("GAMMA_INSERTION", 
									    leftSpinFerm); 
      return foo;
    }

    //! Register the source construction
    const bool registered = registerAll();


    //! Parameters for running code
    Params::Params(XMLReader& xml, const string& path)
    {
      XMLReader paramtop(xml, path);

      read(paramtop, "gamma", gamma);
    }


    //! Parameters for running code
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);
    
      write(xml, "gamma", gamma);

      pop(xml);
    }

  }  // end namespace
}  // end namespace Chroma

