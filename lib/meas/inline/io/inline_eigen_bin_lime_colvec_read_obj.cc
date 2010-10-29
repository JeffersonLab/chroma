// $Id: inline_eigen_bin_lime_colvec_read_obj.cc,v 3.3 2009-07-17 14:53:13 edwards Exp $
/*! \file
 * \brief Inline task to read an object from a named buffer
 *
 * Named object reading
 */

#include "chromabase.h"
//#include "qdp_iogauge.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/inline/io/inline_eigen_bin_lime_colvec_read_obj.h"
#include "meas/inline/io/named_objmap.h"

#include "util/ferm/subset_vectors.h"

#include "util/ferm/map_obj/map_obj_aggregate_w.h"
#include "util/ferm/map_obj/map_obj_factory_w.h"

namespace Chroma 
{ 


  //! Object buffer
  void read(XMLReader& xml, const string& path, InlineEigenBinLimeColVecReadNamedObjEnv::Params::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "object_id", input.object_id);

    // User Specified MapObject tags
    input.object_map = readXMLGroup(inputtop, "ColorVecMapObject", "MapObjType");
  }

  //! File output
  void read(XMLReader& xml, const string& path, InlineEigenBinLimeColVecReadNamedObjEnv::Params::File_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "file_name", input.file_name);
  }


  namespace InlineEigenBinLimeColVecReadNamedObjEnv 
  { 
		
		
//convience routine to return an array given space and time coords
    multi1d<int> coords(const int x, const int y, const int z, const int t)
    {
      multi1d<int> ret(4);

      ret[0] = x;
      ret[1] = y;
      ret[2] = z;
      ret[3] = t;

      return ret;
    }

    //Turn an array into a timeslice of a LatticeColorVector
    void unserialize(LatticeColorVector &cvec, const multi1d< Complex > &vec,
		     const int& t)	
    {
      int vsize = QDP::Layout::lattSize()[0] * QDP::Layout::lattSize()[0] * QDP::Layout::lattSize()[0] * Nc;

      if (vec.size() != vsize) 
      {
	cerr << "in unserialize: invalid size of serialized vector" 
	     << endl;

	exit(0);
      }

      //Loop over space coords 
      for (int x = 0 ; x < QDP::Layout::lattSize()[0] ; ++x)
	for (int y = 0 ; y < QDP::Layout::lattSize()[0] ; ++y)
	  for (int z = 0 ; z < QDP::Layout::lattSize()[0] ; ++z)
	  {
	    ColorVector sitevec = zero;

	    for (int c = 0 ; c < Nc ; ++c)
	    {
	      Complex temp = vec[ c + Nc*(z + QDP::Layout::lattSize()[0]*(y + QDP::Layout::lattSize()[0]*x)) ];

	      pokeColor(sitevec, temp, c);
	    }

	    pokeSite(cvec, sitevec, coords(x, y, z, t) );
	  }

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

      const std::string name = "EIGENINFO_BIN_LIME_COLORVEC_READ_NAMED_OBJECT";
    }

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
    InlineMeas::operator()(unsigned long update_no, XMLWriter& xml_out) 
    {
      START_CODE();

      push(xml_out, "eigeninfo_bin_lime_colorvec_read_named_obj");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": object reader" << endl;
      StopWatch swatch;

      // Read the object
      QDPIO::cout << "Attempt to read object name = " << params.named_obj.object_id << endl;
      write(xml_out, "object_id", params.named_obj.object_id);
      try
      {
	swatch.reset();

	std::istringstream  xml_s(params.named_obj.object_map.xml);
	XMLReader MapObjReader(xml_s);
	
	// Create the entry
	Handle< QDP::MapObject<int,EVPair<LatticeColorVector> > > eigen(
	  TheMapObjIntKeyColorEigenVecFactory::Instance().createObject(params.named_obj.object_map.id,
								       MapObjReader,
								       params.named_obj.object_map.path) );

	TheNamedObjMap::Instance().create< Handle< QDP::MapObject<int,EVPair<LatticeColorVector> > > >(params.named_obj.object_id);
	TheNamedObjMap::Instance().getData< Handle< QDP::MapObject<int,EVPair<LatticeColorVector> > > >(params.named_obj.object_id) = eigen;

	int nt = QDP::Layout::lattSize()[Nd-1];
	int ns = QDP::Layout::lattSize()[0];

	int ndim = ns * ns * ns * Nc;

	// Read the object
	swatch.start();

	//Make the final file xml which will only contain
	//gauge info, and eigen calculation info
	XMLBufferWriter final_file_xml;

	XMLBufferWriter final_record_xml;
	push(final_record_xml, "SubsetVectors");
	push(final_record_xml, "InfoArray");

	XMLReader file_xml;

	std::string filename = params.file.file_name;
	QDPFileReader rdr(file_xml, filename, QDPIO_SERIAL);

	//Plop some info into the file xml only once 
	write(final_file_xml, "Input", file_xml);

#if 1
	multi1d<LatticeColorVector> evecs;
	multi1d<SubsetVectorWeight_t> evals;
#endif
	int nev;
			
	// Weird storage (time slower than vector number)
	for (int t = 0 ; t < nt ; ++t)
	{

	  XMLReader curr_record_xml;
	  BinaryBufferReader bin_rdr;
	  read(rdr, curr_record_xml, bin_rdr);

						
	  // Weights for all vectors on a given timeslice t
	  multi1d<Real> evals_t;

	  if (curr_record_xml.count("/LaplaceEigInfo/EigenValues") != 0)
	    read(curr_record_xml, "/LaplaceEigInfo/EigenValues", evals_t);
	  else if (curr_record_xml.count("/LaplaceEigInfo/EigParams/EigenValues") != 0)
	    read(curr_record_xml, "/LaplaceEigInfo/EigParams/EigenValues", evals_t);
	  else
	  {
	    QDPIO::cerr << __func__ << ": LaplaceEigInfo tag for EigenValues not found\n" << std::endl;
	    QDP_abort(1);
	  }



	  if (t == 0)
	  {
	    nev = evals_t.size();
	    evals.resize(nev);
	    evecs.resize(nev);
	    QDPIO::cout << "Initializing eigenpairs" << endl;
	    for (int v = 0 ; v < nev ; ++v)
	    {
	      evecs[v] = zero;
	      evals[v].weights.resize(nt);
	    }		
	  }
	  QDPIO::cout << "Unserealizing evecs for timeslice " << t << endl;
	  for (int n = 0 ; n < nev ; n++)
	  {

	    evals[n].weights[t] = evals_t[n]; // Copies all the weights

	    multi1d<Complex> temp;
	    read(bin_rdr, temp);

	    if (temp.size() != ndim )
	    {
	      QDPIO::cerr << "Invalid array size" << endl;
	      exit(1);
	    }
	    unserialize(evecs[n], temp, t);
	  }
	  write(final_record_xml, "elem", curr_record_xml);

	}//t

	for(int n=0; n < nev; n++) { 
	  QDPIO::cout << "Inserting eval/evec pair " << n << endl;
	  EVPair<LatticeColorVector> pair;
	  pair.eigenValue = evals[n];
	  pair.eigenVector = evecs[n];
	  eigen->insert(n,pair);
	}
	eigen->flush();

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
