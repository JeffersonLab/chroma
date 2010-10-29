// $Id: inline_create_colorvecs.cc,v 3.11 2009-07-17 14:52:36 edwards Exp $
/*! \file
 * \brief Construct colorvectors via power iteration of the laplacian
 */

#include "fermact.h"
#include "actions/boson/operator/klein_gord.h"
#include "meas/inline/hadron/inline_create_colorvecs.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/gaus_smear.h"
#include "meas/glue/mesplq.h"
#include "util/ferm/map_obj/map_obj_aggregate_w.h"
#include "util/ferm/map_obj/map_obj_factory_w.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlineCreateColorVecsEnv 
  {
    //! Propagator input
    void read(XMLReader& xml, const string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "colorvec_id", input.colorvec_id);

      // User Specified MapObject tags
      input.colorvec_obj = readXMLGroup(inputtop, "ColorVecMapObject", "MapObjType");
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "colorvec_id", input.colorvec_id);
      xml << input.colorvec_obj.xml;

      pop(xml);
    }

    //! Propagator input
    void read(XMLReader& xml, const string& path, Params::Param_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "num_vecs", input.num_vecs);
      read(inputtop, "decay_dir", input.decay_dir);
      read(inputtop, "num_iter", input.num_iter);
      read(inputtop, "num_orthog", input.num_orthog);
      read(inputtop, "width", input.width);
      input.link_smear = readXMLGroup(inputtop, "LinkSmearing", "LinkSmearingType");
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const Params::Param_t& out)
    {
      push(xml, path);

      write(xml, "num_vecs", out.num_vecs);
      write(xml, "decay_dir", out.decay_dir);
      write(xml, "num_iter", out.num_iter);
      write(xml, "num_orthog", out.num_orthog);
      write(xml, "width", out.width);
      xml << out.link_smear.xml;

      pop(xml);
    }


    //! Propagator input
    void read(XMLReader& xml, const string& path, Params& input)
    {
      Params tmp(xml, path);
      input = tmp;
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const Params& input)
    {
      push(xml, path);
    
      write(xml, "Param", input.param);
      write(xml, "NamedObject", input.named_obj);

      pop(xml);
    }
  } // namespace InlineCreateColorVecsEnv 


  namespace InlineCreateColorVecsEnv 
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
      
    const std::string name = "CREATE_COLORVECS";

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


    //-------------------------------------------------------------------------
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
	read(paramtop, "Param", param);

	// Read in the output propagator/source configuration info
	read(paramtop, "NamedObject", named_obj);

	// Possible alternate XML file pattern
	if (paramtop.count("xml_file") != 0) 
	{
	  read(paramtop, "xml_file", xml_file);
	}
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << endl;
	QDP_abort(1);
      }
    }



    // Function call
    void 
    InlineMeas::operator()(unsigned long update_no,
			   XMLWriter& xml_out) 
    {
      // If xml file not empty, then use alternate
      if (params.xml_file != "")
      {
	string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "CreateColorVecs");
	write(xml_out, "update_no", update_no);
	write(xml_out, "xml_file", xml_file);
	pop(xml_out);

	XMLFileWriter xml(xml_file);
	func(update_no, xml);
      }
      else
      {
	func(update_no, xml_out);
      }
    }


    // Real work done here
    void 
    InlineMeas::func(unsigned long update_no,
		     XMLWriter& xml_out) 
    {
      START_CODE();

      StopWatch snoop;
      snoop.reset();
      snoop.start();

      // Test and grab a reference to the gauge field
      multi1d<LatticeColorMatrix> u;
      XMLBufferWriter gauge_xml;
      try
      {
	u = TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
	TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
      }
      catch( std::bad_cast ) 
      {
	QDPIO::cerr << name << ": caught dynamic cast error" << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << name << ": map call failed: " << e << endl;
	QDP_abort(1);
      }

      push(xml_out, "CreateColorVecs");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": Create color vectors" << endl;

      proginfo(xml_out);    // Print out basic program info

      // Write out the input
      write(xml_out, "Input", params);

      // Write out the config header
      write(xml_out, "Config_info", gauge_xml);

      push(xml_out, "Output_version");
      write(xml_out, "out_version", 1);
      pop(xml_out);

      // Calculate some gauge invariant observables just for info.
      MesPlq(xml_out, "Observables", u);

      //
      // Smear the gauge field if needed
      //
      multi1d<LatticeColorMatrix> u_smr = u;
      
      try
      {
	std::istringstream  xml_l(params.param.link_smear.xml);
	XMLReader  linktop(xml_l);
	QDPIO::cout << "Link smearing type = " 
		    << params.param.link_smear.id
		    << endl;
		
	Handle< LinkSmearing >
	  linkSmearing(TheLinkSmearingFactory::Instance().createObject(params.param.link_smear.id, 
								       linktop,params.param.link_smear.path));
	(*linkSmearing)(u_smr);
      }
      catch(const std::string& e){
	QDPIO::cerr << name << ": Caught Exception link smearing: "<<e<< endl;
	QDP_abort(1);
      }

      // Record the smeared observables
      MesPlq(xml_out, "Smeared_Observables", u_smr);


      //
      // Create the output files
      //
      try
      {
	std::istringstream  xml_s(params.named_obj.colorvec_obj.xml);
	XMLReader MapObjReader(xml_s);
	
	// Create the entry
	TheNamedObjMap::Instance().create< Handle< QDP::MapObject<int,EVPair<LatticeColorVector> > > >(params.named_obj.colorvec_id);
	TheNamedObjMap::Instance().getData< Handle< QDP::MapObject<int,EVPair<LatticeColorVector> > > >(params.named_obj.colorvec_id) =
	  TheMapObjIntKeyColorEigenVecFactory::Instance().createObject(params.named_obj.colorvec_obj.id,
								       MapObjReader,
								       params.named_obj.colorvec_obj.path);
      }
      catch (std::bad_cast)
      {
	QDPIO::cerr << name << ": caught dynamic cast error" << endl;
	QDP_abort(1);
      }
      catch (const string& e) 
      {
	QDPIO::cerr << name << ": error creating prop: " << e << endl;
	QDP_abort(1);
      }

      // Cast should be valid now
      MapObject<int,EVPair<LatticeColorVector> >& color_vecs = 
	*(TheNamedObjMap::Instance().getData< Handle< MapObject<int,EVPair<LatticeColorVector> > > >(params.named_obj.colorvec_id));

      // The code goes here
      StopWatch swatch;
      swatch.reset();
      swatch.start();

      /** 
       *Loop for num_iter
       *    apply smearing
       *    if orthonormalize is true to a GramSchmit on them
       *    Compute the v'Smearing v matrix elements for all color vectors
       *    call them evals and monitor convergence on one time slice
       *EndLoop

       **/

      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, params.param.decay_dir);
      
      const int num_vecs = params.param.num_vecs;


      // color_vecs.resizeEvectors(num_vecs);
      multi1d<SubsetVectorWeight_t> evals(num_vecs);
      // Temporary
      multi1d<LatticeColorVector> evecs(num_vecs);


      QDPIO::cout << "Got here" << endl;
      for(int i(0);i<num_vecs;i++){
	evals[i].weights.resize(phases.numSubsets());
      }

      //
      // Initialize the color vectors with gaussian numbers
      // Gaussian smear in an orthogonalization loop
      //
      for(int hit=0; hit <= params.param.num_orthog; ++hit)
      {
	for(int i=0; i < num_vecs; ++i)
	{
	  QDPIO::cout << name << ": Doing colorvec: "<<i << " hit no: "<<hit<<endl;
	  if (hit == 0) {
	    gaussian(evecs[i]);
	  }
	  else {
	    gausSmear(u_smr, 
		      evecs[i],
		      params.param.width, params.param.num_iter, params.param.decay_dir);
	  }

	  for(int k=0; k < i; ++k) {
	    multi1d<DComplex> cc = 
	      sumMulti(localInnerProduct(evecs[k],
					 evecs[i]),
		       phases.getSet());

	    for(int t(0);t<phases.numSubsets();t++) {
	      evecs[i][phases.getSet()[t]] -= cc[t]*evecs[k];
	    }

	  }

	  multi1d<Double> norm2 = sumMulti(localNorm2(evecs[i]),phases.getSet());

	  for(int t=0; t < phases.numSubsets(); ++t) {
	    evecs[i][phases.getSet()[t]] /= sqrt(norm2[t]);
	  }

	}
      }


      //
      // Compute the version of eigenvalues (which they are not)
      // 
      {
	multi1d< multi1d<Double> > source_corrs(evecs.size());
	for(int m=0; m < source_corrs.size(); ++m)
	  source_corrs[m] = sumMulti(localNorm2(evecs[m]), phases.getSet());
	push(xml_out, "Source_correlators");
	write(xml_out, "source_corrs", source_corrs);
	pop(xml_out);
	//compute the "eigenvalues"
	push(xml_out,"SmearingEvals");


	for(int i=0; i < num_vecs; ++i) {
	  LatticeColorVector Svec;
	  klein_gord(u_smr, evecs[i], Svec, Real(0), params.param.decay_dir);

	  multi1d<DComplex> cc = 
	    sumMulti(localInnerProduct(evecs[i], 
				       Svec),  
		     phases.getSet());
	  
	  for(int t=0; t < phases.numSubsets(); ++t) {
	    evals[i].weights[t] = real(cc[t]);
	  }

	  push(xml_out,"Vector");
	  write(xml_out, "VecNo",i);
	  write(xml_out, "Evals", evals[i].weights);
	  pop(xml_out);
	}
	pop(xml_out);
      }
      

      for(int i=0; i < num_vecs; i++) { 
	EVPair<LatticeColorVector> pair;
	pair.eigenValue=evals[i];
	pair.eigenVector=evecs[i];
	color_vecs.insert(i, pair);
      }
      color_vecs.flush();
      
      swatch.stop();
      QDPIO::cout << name << ": time for colorvec construction = "
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;      

      pop(xml_out);  // CreateColorVecs

      /**/
      // Write the meta-data for this operator
      {
	XMLBufferWriter file_xml;

	push(file_xml, "LaplaceEigInfo");
	write(file_xml, "num_vecs", num_vecs); 
	write(file_xml, "Params", params.param);
	write(file_xml, "Config_info", gauge_xml);
	pop(file_xml);

	XMLBufferWriter record_xml;
	push(record_xml, "SubsetVectors");
	for(int i(0);i<num_vecs;i++){
	  push(record_xml, "EigenPair");
	  write(record_xml, "EigenPairNumber", i); 
	  write(record_xml, "EigenValues", evals[i].weights); 
	  pop(record_xml);
	}
	pop(record_xml);
	
	// Write the propagator xml info
	TheNamedObjMap::Instance().get(params.named_obj.colorvec_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.colorvec_id).setRecordXML(record_xml);
      }
      /**/

      snoop.stop();
      QDPIO::cout << name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;

      QDPIO::cout << name << ": ran successfully" << endl;

      END_CODE();
    } 

  }

} // namespace Chroma
