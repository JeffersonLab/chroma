// $Id: dilution_quark_source_const_w.cc,v 1.8 2008-01-18 18:50:05 jbulava Exp $
/*! \file
 * \brief Dilution scheme specified by MAKE_SOURCE and PROPAGATOR calls  
 *
 */

#include "handle.h"
#include "meas/hadron/dilution_quark_source_const_w.h"
#include "meas/hadron/dilution_scheme_factory.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/sources/dilutezN_source_const.h"

namespace Chroma 
{ 

  // Read parameters
  void read(XMLReader& xml, const string& path, DilutionQuarkSourceConstEnv::Params& param)
  {
    DilutionQuarkSourceConstEnv::Params tmp(xml, path);
    param = tmp;
  }

  // Writer
  void write(XMLWriter& xml, const string& path, const DilutionQuarkSourceConstEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  /*!
   * \ingroup hadron
   *
   * 
   */
  namespace DilutionQuarkSourceConstEnv
  { 
		//Read Quark dilution files per timeslice
    void read(XMLReader& xml, const string& path, DilutionQuarkSourceConstEnv::Params::QuarkFiles_t::TimeSliceFiles_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "dilution_files", input.dilution_files);
    }

    //Read Quark timeslice files 
    void read(XMLReader& xml, const string& path, DilutionQuarkSourceConstEnv::Params::QuarkFiles_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "timeslice_files", input.timeslice_files);
    }


 //Write Quark dilution files 
    void write(XMLWriter& xml, const string& path, const DilutionQuarkSourceConstEnv::Params::QuarkFiles_t::TimeSliceFiles_t& input)
    {
      push(xml, path);
      write(xml, "dilution_files", input.dilution_files);
      pop(xml);
    }

    //Write Quark time slice files 
    void write(XMLWriter& xml, const string& path, const DilutionQuarkSourceConstEnv::Params::QuarkFiles_t& input)
    {
      push(xml, path);
      write(xml, "timeslice_files", input.timeslice_files);
      pop(xml);
    }

    //! Initialize
    Params::Params()
    {
    }


    //! Read parameters
    Params::Params(XMLReader& xml, const string& path)
    {
      XMLReader paramtop(xml, path);

      int version;
      read(paramtop, "version", version);

      switch (version) 
      {
      case 1:
				/**************************************************************************/
				break;

			default :
				/**************************************************************************/

				QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
				QDP_abort(1);
			}

      read(paramtop, "QuarkFiles", quark_files);
    }


    // Writer
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);

      int version = 1;
      write(xml, "version", version);
      write(xml, "QuarkFiles", quark_files);
      pop(xml);
    }


  //}


  //namespace DilutionQuarkSourceConstEnv 
  //{ 
    // Anonymous namespace for registration
    namespace
    {
      DilutionScheme<LatticeFermion>* createScheme(XMLReader& xml_in, 
						       const std::string& path) 
      {
				return new ConstDilutionScheme(Params(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "DILUTION_QUARK_SOURCE_CONST_FERM";

		//! Register all the factories
		bool registerAll() 
		{
			bool success = true; 
			if (! registered)
			{
				success &= TheFermDilutionSchemeFactory::Instance().registerObject(name, createScheme);
				registered = true;
			}
			return success;
		}



		bool operator!= (const QuarkSourceSolutions_t::TimeSlices_t::Dilutions_t & dilA,
			const QuarkSourceSolutions_t::TimeSlices_t::Dilutions_t & dilB)
		{
			bool val = false;

			//Check if the seeds are the same 
			Seed seedA, seedB;
			
			read(dilA.source_header, "ran_seed" , seedA);
			read(dilB.source_header, "ran_seed" , seedB);

			val |= (seedA != seedB);

			
			//Check Spatial mask and spatial mask size
			
			multi1d<int> mask_sizeA, mask_sizeB, maskA, maskB;

			read(dilA.source_header, "/spatial_mask_size" , mask_sizeA);
			read(dilB.source_header, "/spatial_mask_size" , mask_sizeB);
			val |= (mask_sizeA != mask_sizeB);

			read(dilA.source_header, "/spatial_mask" , maskA);
			read(dilB.source_header, "/spatial_mask" , maskB);
			val |= (maskA != maskB);

			read(dilA.source_header, "/color_mask" , maskA);
			read(dilB.source_header, "/color_mask" , maskB);
			val |= (maskA != maskB);

			read(dilA.source_header, "/spin_mask" , maskA);
			read(dilB.source_header, "/spin_mask" , maskB);
			val |= (maskA != maskB);

			return val;

		}
    //-------------------------------------------------------------------------------
    // Function call
    void ConstDilutionScheme::init()
    {
      START_CODE();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      //
      // Read the soultion headers
      //
      StopWatch swatch;
      swatch.reset();
      swatch.start();

      try
      {

				bool initq = false;

				QDPIO::cout << "Attempt to read solutions " << endl;
		
				quark.time_slices.resize( params.quark_files.timeslice_files.size() );

				QDPIO::cout<< "time_slices.size = " << quark.time_slices.size() << endl;

				
				for (int t0 = 0 ; t0 < quark.time_slices.size() ; ++t0 )
				{
					quark.time_slices[t0].dilutions.resize(
							params.quark_files.timeslice_files[t0].dilution_files.size() );
				
					QDPIO::cout << "dilutions.size = " << 
						quark.time_slices[t0].dilutions.size() << endl;
					
					int time = 0;

					for(int dil = 0; dil < quark.time_slices[t0].dilutions.size(); ++dil)
					{				
						quark.time_slices[t0].dilutions[dil].soln_file =
							params.quark.timeslice_files[t0].dilution_files[dil];

	    			XMLReader file_xml, record_xml, record_xml_source;

	    			QDPIO::cout << "reading file = " << 
							quark.time_slices[t0].dilutions[dil].soln_file << endl;

	    			QDPFileReader from(file_xml, 
								quark.time_slices[t0].dilutions[dil].soln_file, QDPIO_SERIAL);

						//Use the record xml only, throw away the lattice fermion
						LatticeFermion junk;
						read(from, record_xml, junk);
	    			close(from);

	    			read(record_xml, "/Propagator/PropSource", 
								quark.dilutions[dil].source_header);
	    			read(record_xml, "/Propagator/ForwardProp", 
								quark.dilutions[dil].prop_header);


	    			if (!initq)
	    			{
	      			read(record_xml, "/Propagator/PropSource/Source/ran_seed",
		   				quark.seed);
					
	    				quark.decay_dir = 
								quark.time_slices[t0].dilutions[dil].source_header.j_decay;

	      			initq = true;
	    			}

	   
						//Test that t0 is the same for all dilutions on this timeslice
						//grab the first t0 to prime the process
					
						if (dil == 0)
						{
							time = quark.time_slices[t0].dilutions[dil].source_header.t0;
						}

						if (time != quark.time_slices[t0].dilutions[dil].source_header.t0)
						{
							QDPIO::cerr << "t0's DO NOT MATCH FOR ALL DILUTIONS ON TIME SLICE "
							<< t0 << endl;

							QDP_abort(1);
						}

					
						//Test that each t0 has the same dilutions per timeslice 
						if ( quark.time_slices[t0].dilutions[dil] != 
								quark.time_slices[0].dilutions[dil] )
						{
							QDPIO::cerr << "Dilutions do not match on time slice " << 
								t0 << endl;

							QDP_abort(1);
						}

	  			}//dil
					
					quark.time_slices[t0].t0 = time;

					

				}//t0

				

      } //try
      catch (const string& e) 
      {
	QDPIO::cerr << "Error extracting headers: " << e << endl;
	QDP_abort(1);
      }
      swatch.stop();

      QDPIO::cout << "Source and solution headers successfully read: time = "
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;


      END_CODE();
    } // init


		
		//Create and return the diluted source
		//For gauge invariance testing reasons, this routine could be changed
		//to get the source from the named object map
		LatticeFermion ConstDilutionScheme::dilutedSource(int t0, int dil ) const 
		{
			QuarkSourceSolutions_t::TimeSlices_t::Dilutions_t &qq = 
				quark.time_slices[t0].dilutions[dil];

			//Dummy gauge field to send to the source construction routine;
			multi1d<LatticeColorMatrix> dummy = zero;
			

			// Build source construction
	    QDPIO::cout << "Source_id = " << qq.source_header.source.id << endl;
	    QDPIO::cout << "Source = XX" << qq.source_header.source.xml << "XX" << endl;

	    std::istringstream  xml_s(qq.source_header.source.xml);
	    XMLReader  sourcetop(xml_s);

	    if (qq.source_header.source.id != DiluteZNQuarkSourceConstEnv::name)
	    {
				QDPIO::cerr << "Expected source_type = " << DiluteZNQuarkSourceConstEnv::name << endl;
				QDP_abort(1);
	    }

	    // Manually create the params so I can peek into them and use the source constructor
	    DiluteZNQuarkSourceConstEnv::Params  srcParams(sourcetop, 
							     qq.source_header.source.path);
	    DiluteZNQuarkSourceConstEnv::SourceConst<LatticeFermion>  srcConst(srcParams);
      
			QDP::RNG::setrn( quark.seed );

			 
			return srcConst(dummy);
	      
		} //dilutedSource		


		LatticeFermion ConstDilutionScheme::dilutedSolution(int t0, int dil) const 
		{
			
			const std::string & filename = quark.time_slices[t0].dilutions[dil].soln_file;

			LatticeFermion soln; 

		  XMLReader file_xml, record_xml, record_xml_source;

	    QDPIO::cout << "reading file = " << filename << endl;
	    QDPFileReader from(file_xml, filename, QDPIO_SERIAL);

			read(from, record_xml, soln);

			return soln;
		}


	
  } // namespace DilutionQuarkSourceConstEnv

}// namespace Chroma 
