// $Id: inline_laplace_eigs.cc,v 1.2 2009-06-15 18:04:29 jbulava Exp $
/*! \file
 * \brief Use the IRL method to solve for eigenvalues and eigenvectors 
 * of the gauge-covariant laplacian.  
 */

#include "fermact.h"
#include "meas/inline/hadron/inline_laplace_eigs.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/smear/link_smearing_factory.h"
#include "meas/smear/link_smearing_aggregate.h"
#include "meas/smear/gaus_smear.h"
#include "meas/glue/mesplq.h"
#include "util/ferm/subset_vectors.h"
#include "util/ferm/map_obj.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlineLaplaceEigsEnv 
  {
    //! Propagator input
    void read(XMLReader& xml, const string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
      read(inputtop, "colorvec_id", input.colorvec_id);
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);
      write(xml, "colorvec_id", input.colorvec_id);

      pop(xml);
    }

    //! Propagator input
    void read(XMLReader& xml, const string& path, Params::Param_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "num_vecs", input.num_vecs);
      read(inputtop, "decay_dir", input.decay_dir);
      read(inputtop, "max_iter", input.max_iter);
      read(inputtop, "tol", input.tol);
      input.link_smear = readXMLGroup(inputtop, "LinkSmearing", "LinkSmearingType");
    }

    //! Propagator output
    void write(XMLWriter& xml, const string& path, const Params::Param_t& out)
    {
      push(xml, path);

      write(xml, "num_vecs", out.num_vecs);
      write(xml, "decay_dir", out.decay_dir);
      write(xml, "max_iter", out.max_iter);
      write(xml, "tol", out.tol);
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


  namespace InlineLaplaceEigsEnv 
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
      
    const std::string name = "LAPLACE_EIGS";

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

	push(xml_out, "LaplaceEigs");
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

void gram(const multi1d<LatticeColorVector>& init, multi1d<LatticeColorVector>& ortho)
	{
	for(int i = 0; i < init.size() ; i++){
	  ortho[i] = init[i];
	  for(int j = 0; j < i; j++){
	    ortho[i] -= innerProduct(ortho[j], init[i])/innerProduct(ortho[j], ortho[j]) * ortho[j];
	  }
	}
	//still need to normalize (take another parameter of LCV's?)
	
	
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

      push(xml_out, "LaplaceEigs");
      write(xml_out, "update_no", update_no);

      QDPIO::cout << name << ": Use the IRL method to solve for laplace eigenpairs" << endl;

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
	TheNamedObjMap::Instance().create< SubsetVectors<LatticeColorVector> >(params.named_obj.colorvec_id);
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
      SubsetVectors<LatticeColorVector>& color_vecs =
	TheNamedObjMap::Instance().getData< SubsetVectors<LatticeColorVector> >(params.named_obj.colorvec_id);

      // The code goes here
      StopWatch swatch;
      swatch.reset();
      swatch.start();

      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, params.param.decay_dir);
      
      const int num_vecs = params.param.num_vecs;

      color_vecs.getEvectors().resize(num_vecs);
      color_vecs.getEvalues().resize(num_vecs);
      for(int i(0);i<num_vecs;i++){
	color_vecs.getEvalues()[i].weights.resize(phases.numSubsets());
      }

      color_vecs.getDecayDir() = params.param.decay_dir;


	  //Choose starting vector 
	  LatticeColorVector start;

	
	  //Build Krlov subspace
	  int kdim = 3 * params.param.num_vecs;
	 
	  multi1d<LatticeColorVector> krylov_vecs(kdim);

	  krylov_vecs[0] = start;

	  for (int k = 1 ; k < kdim  ; ++k)
	  {	
		 //Perform orthogonalization wrt previous vector
		 orthog(u, krylov_vecs[k-1], krylov_vecs[k]);
		  
	  	//Reorthogonalize?
		//gram(krylov,krylov,k);

		 if (k > 1)
		 {
			LatticeColorVector beta_km1; 
			getBeta(beta_km1, krylov_vecs[k-1]);

			 krylov_vecs[k] -= krylov_vecs[k-1] * beta_km1; 
		 }

		 LatticeColorVector beta_k;
		 getBeta(beta_k, krylov_vecs[k]);

		 if (norm2(beta_k) < params.param.tol)
		 {
			 break;
		 }
		 
		 krylov_vecs /= beta_k; 
	  }

	  	





	  //BRANDON: CODE GOES HERE
      /*
	v1 = v/||v|| - choose starting vector in direction of interest
	Lanczos factorization - AVm = VmTm + rmem

	Loop to convergence - (Tk = Dk)
	compute o(Tm) and p shifts (u1, u2, up)
	initialize Q = Im
	
	for j = 1 to p{
	QR factorize: QjRj = Tm - ujI
	update Tm = Qj*TmQj, Q=QQj
	}

	rk = vk+1Bk + rmok, where Bk = Tm(k+1,k) and ok = Q(m,k)
	Vk = VmQ(:, 1 : k); Tk = Tm(1 : k, 1 : k)
	k step Lanczos factorization
	AVk = VkTk + rkek*
	apply p additional steps of Lanczos to obtain a new mstep Lanczos
	AVm = VmTm + rmem*
	end

	Aside: m step Lanczos
	start with r = rk, prev residual or start vec, Bk = ||rk||
	for j = k+1, m {
	vj = r/Bj-1
	r = Avj - vj-1/Bj-1
	aj = bj*r
	r = r - vjaj
	if(||r|| < p*sqrt(aj^2 + Bj-1^2){
	s = Vj*r
	r = r - Vjs
	aj = aj + sj, Bj = Bj + sj-1
}
      */

		  pop(xml_out);
      
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
	  write(record_xml, "EigenValues", color_vecs.getEvalues()[i].weights); 
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
