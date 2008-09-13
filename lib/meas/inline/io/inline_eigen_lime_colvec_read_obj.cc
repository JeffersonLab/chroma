// $Id: inline_eigen_lime_colvec_read_obj.cc,v 3.2 2008-09-13 19:56:40 edwards Exp $
/*! \file
 * \brief Inline task to read an object from a named buffer
 *
 * Named object reading
 */

#include "chromabase.h"
#include "qdp_iogauge.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/inline_eigen_lime_colvec_read_obj.h"
#include "meas/inline/io/named_objmap.h"

#include "util/ferm/subset_vectors.h"

namespace Chroma 
{ 


  //! Object buffer
  void write(XMLWriter& xml, const string& path, const InlineEigenLimeColVecReadNamedObjEnv::Params::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "object_id", input.object_id);

    pop(xml);
  }

  //! File output
  void write(XMLWriter& xml, const string& path, const InlineEigenLimeColVecReadNamedObjEnv::Params::File_t& input)
  {
    push(xml, path);

    write(xml, "file_names", input.file_names);

    pop(xml);
  }


  //! Object buffer
  void read(XMLReader& xml, const string& path, InlineEigenLimeColVecReadNamedObjEnv::Params::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "object_id", input.object_id);
  }

  //! File output
  void read(XMLReader& xml, const string& path, InlineEigenLimeColVecReadNamedObjEnv::Params::File_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "file_names", input.file_names);
  }


  namespace InlineEigenLimeColVecReadNamedObjEnv 
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

    const std::string name = "EIGENINFO_LIME_COLORVEC_READ_NAMED_OBJECT";

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
	read(paramtop, "NamedObject", named_obj);

	// Read in the destination
	read(paramtop, "File", file);
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
      write(xml_out, "NamedObject", named_obj);

      // Write out the destination
      write(xml_out, "File", file);

      pop(xml_out);
    }


    void 
    InlineMeas::operator()(unsigned long update_no, XMLWriter& xml_out) 
    {
      START_CODE();

      push(xml_out, "eigeninfo_lime_colorvec_read_named_obj");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": object reader" << endl;
      StopWatch swatch;

      // Read the object
      QDPIO::cout << "Attempt to read object name = " << params.named_obj.object_id << endl;
      write(xml_out, "object_id", params.named_obj.object_id);
      try
      {
	swatch.reset();
	
	typedef LatticeColorVector T;

	TheNamedObjMap::Instance().create< SubsetVectors<T> >(params.named_obj.object_id);
	SubsetVectors<T>& eigen = TheNamedObjMap::Instance().getData< SubsetVectors<T> >(params.named_obj.object_id);

	eigen.getEvalues().resize(params.file.file_names.size());
	eigen.getEvectors().resize(params.file.file_names.size());
	eigen.getDecayDir() = Nd-1;

	// Read the object
	swatch.start();

	//Make the final file xml which will only contain
	//gauge info, and eigen calculation info
	XMLBufferWriter final_file_xml;

	XMLBufferWriter final_record_xml;
	push(final_record_xml, "SubsetVectors");
	push(final_record_xml, "InfoArray");

	for(int i=0; i < params.file.file_names.size(); ++i)
	{
	  XMLReader curr_file_xml;
	  XMLReader curr_record_xml;

	  std::string filename = params.file.file_names[i];
	  QDPFileReader rdr(curr_file_xml, filename, QDPIO_SERIAL);

	  read(rdr, curr_record_xml, eigen.getEvectors()[i]);
		
	  multi1d<Real> evals;
	  read(curr_record_xml, "/LaplaceEigInfo/EigenValues", evals);

	  eigen.getEvalues()[i].weights = evals; // Copies all the weights

	  write(final_record_xml, "elem", curr_record_xml);

	  //Plop some info into the file xml only once 
	  if (i == 0 )
	  {
	    write(final_file_xml, "Input", curr_file_xml);
	  }
	
	}

	pop(final_record_xml);
	pop(final_record_xml);
	
	swatch.stop();

	TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(final_file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(final_record_xml);

	QDPIO::cout << "Object successfully read: time= " 
		    << swatch.getTimeInSeconds() 
		    << " secs" << endl;
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << name << ": cast error" << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << name << ": error message: " << e << endl;
	QDP_abort(1);
      }
    
      QDPIO::cout << name << ": ran successfully" << endl;

      pop(xml_out);  // read_named_obj

      END_CODE();
    } 

  }
}
