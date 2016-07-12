/*! \file
 * \brief Inline measurement of qpropadd
 *
 * Addition of props
 */

#include "meas/inline/hadron/inline_qpropdiff_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"

#include "meas/inline/io/named_objmap.h"
#include "util/ferm/transf.h"

namespace Chroma 
{ 
  //! Propagator parameters
  void read(XMLReader& xml, const std::string& path, InlineQpropDiffEnv::Params::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);


    read(inputtop, "propA", input.propA);
    read(inputtop, "propB", input.propB);

  }

  //! Propagator parameters
  void write(XMLWriter& xml, const std::string& path, const InlineQpropDiffEnv::Params::NamedObject_t& input)
  {
    push(xml, path);


    write(xml, "propA", input.propA);
    write(xml, "propB", input.propB);

    pop(xml);
  }


  namespace InlineQpropDiffEnv 
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

    const std::string name = "QPROP_DIFF";

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
	// Read in the output propagator/source configuration info
	read(paramtop, "NamedObject", named_obj);
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
	QDP_abort(1);
      }
    }


    void
    Params::writeXML(XMLWriter& xml_out, const std::string& path) 
    {
      push(xml_out, path);
    
      // Write out the output propagator/source configuration info
      write(xml_out, "NamedObject", named_obj);

      pop(xml_out);
    }

    //--------------------------------------------------------------


    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      START_CODE();

      push(xml_out, "qprop_diff");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << "QPROP_DIFF: Check differnce of two props component by component" << std::endl;

      // Write out the input
      params.writeXML(xml_out, "Input");

      //
      // Read in the source along with relevant information.
      // 
      XMLReader propA_file_xml, propA_record_xml;
    
      LatticePropagator propA ;
      LatticePropagator propB ;

      QDPIO::cout << "Snarf the props from a named buffer" << std::endl;
      try
      {
	// Try the cast to see if this is a valid source
	propA = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.propA);

	TheNamedObjMap::Instance().get(params.named_obj.propA).getFileXML(propA_file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.propA).getRecordXML(propA_record_xml);

	propB = TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.propB);
      }    
      catch (std::bad_cast)
      {
	QDPIO::cerr << name << ": caught dynamic cast error" 
		    << std::endl;
	QDP_abort(1);
      }
      catch (const std::string& e) 
      {
	QDPIO::cerr << name << ": error extracting source_header: " << e << std::endl;
	QDP_abort(1);
      }

      for(int spin=0; spin < Ns; ++spin){
    	  for(int color=0; color < Nc; ++color ) {
    			LatticeFermion fermA=zero;
    			LatticeFermion fermB=zero;
    			PropToFerm(propA, fermA, color, spin);
    			PropToFerm(propB, fermB, color, spin);
    			LatticeFermion diff = fermA - fermB;
    			Double diff_L2=sqrt(norm2(diff));
    		//	Double diff_Linf=toDouble(globalMax(diff));
    			Double A_L2 = sqrt(norm2(fermA));
    		//	Double A_Linf=toDouble(globalMax(fermA));
    			Double B_L2 = sqrt(norm2(fermB));
    		//	Double B_Linf=toDouble(globalMax(fermB));
    			QDPIO::cout << "QPROP_DIFF: spin="<<spin<<" col="<<color
    						<< " ||diff||="<<diff_L2
							<< " ||diff||/||A||=" << diff_L2/A_L2
							<< " ||diff||/||B||="<<diff_L2/B_L2 << std::endl;

#if 0
    	      	QDPIO::cout << "QPROP_DIFF: spin="<<spin<<" col="<<color
			<< " ||diff||_inf="<<diff_Linf
			<< " ||diff||_inf/||A||_inf=" << diff_Linf/A_Linf
			<< " ||diff||_inf/||B||_inf="<<diff_Linf/B_Linf << std::endl;

#endif
    	      	QDPIO::cout << std::endl;
    	  }
      }


      pop(xml_out);   // qpropdiff
        
      QDPIO::cout << "QPROP_DIFF: ran successfully" << std::endl;

      END_CODE();
    }

  }

}  
