// $Id: inline_eigen_bin_colvec_read_obj.cc,v 3.4 2008-09-13 19:56:40 edwards Exp $
/*! \file
 * \brief Inline task to read an object from a named buffer
 *
 * Named object reading
 */

#include "chromabase.h"
#include "qdp_iogauge.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/inline_eigen_bin_colvec_read_obj.h"
#include "meas/inline/io/named_objmap.h"

#include "util/ferm/subset_vectors.h"

#include "util/ferm/map_obj/map_obj_aggregate_w.h"
#include "util/ferm/map_obj/map_obj_factory_w.h"

namespace Chroma 
{ 


  //! Object buffer
  void write(XMLWriter& xml, const string& path, const InlineEigenBinColVecReadNamedObjEnv::Params::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "object_id", input.object_id);
    xml << input.object_map.xml;

    pop(xml);
  }

  //! File output
  void write(XMLWriter& xml, const string& path, const InlineEigenBinColVecReadNamedObjEnv::Params::File_t& input)
  {
    push(xml, path);

    write(xml, "file_names", input.file_names);

    pop(xml);
  }


  //! Object buffer
  void read(XMLReader& xml, const string& path, InlineEigenBinColVecReadNamedObjEnv::Params::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "object_id", input.object_id);

    // User Specified MapObject tags
    input.object_map = readXMLGroup(inputtop, "ColorVecMapObject", "MapObjType");
  }

  //! File output
  void read(XMLReader& xml, const string& path, InlineEigenBinColVecReadNamedObjEnv::Params::File_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "file_names", input.file_names);
  }


  namespace InlineEigenBinColVecReadNamedObjEnv 
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

    const std::string name = "EIGENINFO_BIN_COLORVEC_READ_NAMED_OBJECT";

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

      push(xml_out, "eigeninfo_bin_colorvec_read_named_obj");
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
	for(int i=0; i < params.file.file_names.size(); ++i)
	{
	  BinaryFileReader bin(params.file.file_names[i]);

	  EVPair<T> read_pair;
	  // Read the vector
	  read(bin, read_pair.eigenVector);

	  // Resize the eigenvalue.weights array
	  read_pair.eigenValue.weights.resize(Lt);

	  // Read the weights
	  for(int t=0; t < Lt; ++t)
	  {
	    read(bin, read_pair.eigenValue.weights[t]);
	  }

	  // Insert into new map
	  eigen->insert(i,read_pair);

	}
	eigen->flush();

	swatch.stop();

	XMLBufferWriter file_xml;
	push(file_xml, "SubsetVectors");
	pop(file_xml);

	XMLBufferWriter record_xml;
	push(record_xml, "SubsetVectors");
	pop(record_xml);

	TheNamedObjMap::Instance().get(params.named_obj.object_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.object_id).setRecordXML(record_xml);

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
