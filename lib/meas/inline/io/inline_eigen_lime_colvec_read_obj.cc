// $Id: inline_eigen_lime_colvec_read_obj.cc,v 3.5 2008-10-30 14:47:19 edwards Exp $
/*! \file
 * \brief Inline task to read an object from a named buffer
 *
 * Named object reading
 */

#include "chromabase.h"
//#include "qdp_iogauge.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/inline_eigen_lime_colvec_read_obj.h"
#include "meas/inline/io/named_objmap.h"

#include "util/ferm/subset_vectors.h"

#include "util/ferm/map_obj/map_obj_aggregate_w.h"
#include "util/ferm/map_obj/map_obj_factory_w.h"

namespace Chroma 
{ 
  namespace InlineEigenLimeColVecReadNamedObjEnv 
  { 
    //! Object buffer
    void write(XMLWriter& xml, const string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "object_id", input.object_id);
      xml << input.object_map.xml;

      pop(xml);
    }

    //! File output
    void write(XMLWriter& xml, const string& path, const Params::File_t& input)
    {
      push(xml, path);
      
      write(xml, "file_names", input.file_names);

      pop(xml);
    }


    //! Object buffer
    void read(XMLReader& xml, const string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "object_id", input.object_id);

      // User Specified MapObject tags
      input.object_map = readXMLGroup(inputtop, "ColorVecMapObject", "MapObjType");
    }

    //! File output
    void read(XMLReader& xml, const string& path, Params::File_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "file_names", input.file_names);
    }


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
	success &= MapObjectWilson4DEnv::registerAll();
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

	std::istringstream  xml_s(params.named_obj.object_map.xml);
	XMLReader MapObjReader(xml_s);
	
	// Create the entry
	Handle< QDP::MapObject<int,EVPair<LatticeColorVector> > > eigen(
	  TheMapObjIntKeyColorEigenVecFactory::Instance().createObject(params.named_obj.object_map.id,
								       MapObjReader,
								       params.named_obj.object_map.path) );

	TheNamedObjMap::Instance().create< Handle< QDP::MapObject<int,EVPair<LatticeColorVector> > > >(params.named_obj.object_id);
	TheNamedObjMap::Instance().getData< Handle< QDP::MapObject<int,EVPair<LatticeColorVector> > > >(params.named_obj.object_id) = eigen;

	// Argh, burying this in here
	const int decay_dir = Nd-1;
	const int Lt = QDP::Layout::lattSize()[decay_dir];

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
	  EVPair<T> readpair;

	  read(rdr, curr_record_xml, readpair.eigenVector);
	  if (curr_record_xml.count("/LaplaceEigInfo/EigenValues") != 0)
	    read(curr_record_xml, "/LaplaceEigInfo/EigenValues", readpair.eigenValue.weights);
	  else if (curr_record_xml.count("/LaplaceEigInfo/EigParams/EigenValues") != 0)
	    read(curr_record_xml, "/LaplaceEigInfo/EigParams/EigenValues", readpair.eigenValue.weights);
	  else
	  {
	    QDPIO::cerr << __func__ << ": LaplaceEigInfo tag for EigenValues not found\n" << std::endl;
	    QDP_abort(1);
	  }

	  eigen->insert(i,readpair);

	  write(final_record_xml, "elem", curr_record_xml);

	  //Plop some info into the file xml only once 
	  if (i == 0 )
	  {
	    write(final_file_xml, "Input", curr_file_xml);
	  }
	
	}

	pop(final_record_xml);
	pop(final_record_xml);
	eigen->flush();

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
