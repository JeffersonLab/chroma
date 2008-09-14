// $Id: inline_rng.cc,v 3.1 2008-09-14 02:26:11 edwards Exp $
/*! \file
 * \brief Inline task to read and write RNG seed
 */

#include "chromabase.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/inline_rng.h"

namespace Chroma 
{ 
  namespace InlineSetRNGEnv 
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineMeas(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "SET_RNG";

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
    // Param stuff
    Params::Params() { frequency = 0; }

    Params::Params(XMLReader& xml_in, const std::string& path) 
    {
      try 
      {
	XMLReader paramtop(xml_in, path);

	if (paramtop.count("Frequency") == 1)
	  read(paramtop, "Frequency", frequency);
	else
	  frequency = 1;

	// Parameters for source construction
	read(paramtop, "RNG", ran_seed);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << __func__ << ": caught Exception reading XML: " << e << endl;
	QDP_abort(1);
      }
    }


    void
    Params::writeXML(XMLWriter& xml_out, const std::string& path) 
    {
      push(xml_out, path);
    
      // Parameters for source construction
      write(xml_out, "RNG", ran_seed);

      pop(xml_out);
    }


    void 
    InlineMeas::operator()(unsigned long update_no, XMLWriter& xml_out) 
    {
      START_CODE();

      push(xml_out, "SetRNG");
      write(xml_out, "update_no", update_no);

      // Set the RNG seed
      QDP::RNG::setrn(params.ran_seed);
      write(xml_out,"RNG", params.ran_seed);
    
      QDPIO::cout << name << ": ran successfully" << endl;

      pop(xml_out);  // read_named_obj

      END_CODE();
    } 

  }
}
