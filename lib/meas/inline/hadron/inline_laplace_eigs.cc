// $Id: inline_laplace_eigs.cc,v 1.4 2009-06-23 15:12:42 jbulava Exp $
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
#include "actions/boson/operator/klein_gord.h"

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

    /*
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
    */  
   
     template<typename T>                                                                                                                                           
     void partitionedInnerProduct(const T& phi, const T& chi, multi1d<DComplex>& inner_prod, const Set& product_set){                       
                                                                                                                                                                    
          inner_prod = sumMulti(localInnerProduct(phi,chi),product_set);                                                                                            
     }
 
    
    /*
     void partitionedInnerProduct(const LatticeColorVector& phi, const LatticeColorVector& chi, multi1d<DComplex>& inner_prod, Set& product_set)
     {
     
       inner_prod = sumMulti(localInnerProduct(phi,chi),product_set);
     
     }
*/

     template<typename T>                                                                                                                                           
     void laplacian(const multi1d<LatticeColorMatrix>& u, const T& psi, T& chi, int j_decay)
    {                                                                                                                                                               
          T temp;                                                                                                                                                   
                                                                                                                                                                    
          Real minus_one = -1.;                                                                                                                                     
                                                                                                                                                                    
          temp = chi * minus_one;                                                                                                                                   
                                                                                                                                                                    
          klein_gord(u, temp, chi, 0.0, j_decay);                                                                                                                   
    }  


    void q(const multi1d<LatticeColorMatrix>& u,
	   const LatticeColorVector& psi,
	   LatticeColorVector& chi,
	   int j_decay)
    {
      laplacian(u, psi, chi, j_decay);
      chi *= Real(-2.0 / 14.0);
      chi += psi * Real(-1.0 - 2.0 * .35 / 14.0);
    }

    //12th order chebyshev
    void chebyshev(const multi1d<LatticeColorMatrix>& u,
		   const LatticeColorVector& psi,
		   LatticeColorVector& chi,
		   int j_decay)
    {
      int n = 12;                                                                                                                                               
      double chebCo[6] = {-72.0, 840.0, -3584.0, 6912.0, -6144.0, 2048.0};
      LatticeColorVector tmp = psi;
      LatticeColorVector prev;
      LatticeColorVector final = zero;

      for(int i = 2; i <= n; i += 2){
	if(i > 2)
	  tmp = prev;

	q(u, tmp, chi, j_decay);
	tmp = chi;
	q(u, tmp, chi, j_decay);
	prev = chi;

	final += chebCo[i/2-1]*chi;
      }

      final += psi;
      chi = final;
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
      int nt = phases.numSubsets();

      color_vecs.getEvectors().resize(num_vecs);
      color_vecs.getEvalues().resize(num_vecs);
      for(int i(0);i<num_vecs;i++){
	color_vecs.getEvalues()[i].weights.resize(phases.numSubsets());
      }

      color_vecs.getDecayDir() = params.param.decay_dir;

		// Choose the starting eigenvectors to have identical 
		// components and unit norm. 
		// The norm is evaluated time slice by time slice
		LatticeColorVector starting_vectors;
		ColorVector ones;
		pokeColor(ones, Complex(1.0), 0);
		pokeColor(ones, Complex(1.0), 1);
		pokeColor(ones, Complex(1.0), 2);
	
		starting_vectors = ones;

		// Norm of the eigenvectors on each time slice
		// vector_norms is declared complex, this allows 
		// us to use partitionedInnerProduct but it may be a
		// design flaw
		multi1d<DComplex> vector_norms;

		// This function gives the norms squared
		partitionedInnerProduct(starting_vectors,starting_vectors,vector_norms,phases.getSet());
		// Apply the square root to get the true norm
		// and normalise the starting vectors
		for(int t=0; t<nt; ++t)
		{
		  vector_norms[t]  = Complex(sqrt(Real(real(vector_norms[t]))));
			starting_vectors[phases.getSet()[t]] /= vector_norms[t];
		}


	  //Build Krlov subspace
	  int kdim = 3 * params.param.num_vecs;
	  int j_decay = params.param.decay_dir;

		// beta should really be an array of Reals	
		multi1d< multi1d<DComplex> > beta(kdim-1);	
		multi1d< multi1d<DComplex> > alpha(kdim);

		for (int k = 0 ; k < kdim ; ++k)
		{			
			alpha[k].resize(nt);

			if (k < kdim-1)
			{
				beta[k].resize(nt);
			}
		}


		multi1d<LatticeColorVector> lanczos_vectors(kdim);
		lanczos_vectors[0] = starting_vectors;

		// Yields alpha[0] ... alpha[kdim-2]
		// 				beta[0] ... beta[kdim-2] 
		// 				lanczos_vector[0] ... lanczos_vector[kdim-1]
		// After the last iteration compute alpha[kdim-1]
		for(int k=0; k<kdim-1; ++k)
		{
		        //temporary seems to be defined as a single element but is used as both an array and a single element?
			LatticeColorVector temporary;
			// Apply the spatial Laplace operator; j_decay denotes the temporal direction		
			laplacian(u,lanczos_vectors[k],temporary,j_decay); 
			
			if(k > 0){	
			 	for(int t=0; t<nt; ++t){
				  temporary[phases.getSet()[t]] -= beta[k-1][t]*lanczos_vectors[k-1];
				}
			}

			partitionedInnerProduct(lanczos_vectors[k],temporary,alpha[k],phases.getSet());
			
			for(int t=0; t<nt; ++t){
			  temporary[phases.getSet()[t]] -= alpha[k][t]*lanczos_vectors[k];
			}


			// Reorthogonalise - this may be unnecessary
			if(k>0){
			 for(int t=0; t<nt; ++t){
			   temporary[phases.getSet()[t]] -= alpha[k-1][t]*lanczos_vectors[k-1];
			 }
			}

			for(int t=0; t<nt; ++t){
			  temporary[phases.getSet()[t]] -= alpha[k][t]*lanczos_vectors[k];
			} //
			
		
			// Global reorthogonalisation to go here?	
			// .......
			// .....	

			partitionedInnerProduct(temporary,temporary,beta[k],phases.getSet());

			for(int t=0; t<nt; ++t)
			{
			  beta[k][t] 			= Complex(sqrt(Real(real(beta[k][t]))));
				lanczos_vectors[k+1][phases.getSet()[t]] = temporary/beta[k][t];
			}
		}
		// Loop over k is complete, now compute alpha[kdim-1]

		// Finally compute eigenvectors and eigenvalues




		
		/*
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
		*/
		
		// On a gauge field configuration {U}, the action of the
		// three dimensional laplacian on a colour field (eigenvector) Psi
		// is 
		// -2*(Nd-1)*Psi +  
		// sum{mu ne j_decay} [ U_mu(x) * Psi(x+mu)  +
		// 										  U^{dagger}_mu(x-mu) * Psi(x-mu) ]
		//
		// \param j_decay  'temporal direction'
		
		// We will apply the Laplacian using the 
		// "klein_gord" function :
		// void klein_gord(const multi1d<LatticeColorMatrix>& u,
		// 								 const T& psi,
		// 								 T& chi,
		// 								 const Real& mass_sq, int j_decay)
		//
		// where 
		// Chi := (mass_sq + 2*(Nd-1))*Psi -
		// 				sum_{mu ne j_decay} [U_mu(x) * Psi(x+mu) +
		// 														 U^dagger_mu(x-mu) * Psi(x-mu) ]
		//
		// \param  mass_sq 'klein-gordan mass*mass' 
		//
		// To get the action of the laplacian on  
		// the field chi, we write the (templated) function:
		//
		// template<typename T>
		// void laplacian(const multi1d<LatticeColorMatrix>& u,
		// 								T& psi, T& chi, int j_decay)
		//{
		//	T temp;
		//
		// 	Real minus_one = -1.;
		//
		// 	temp = chi * minus_one;
		//	
		//	klein_gord(u, temp, chi, 0.0, j_decay);
		//}
		//
		//
		//
		// It is also useful to define a function which evaluates the inner 
		// product of two fields time slice by time slice. 
		// The function should return an array of complex numbers. 
		// Each element of the array will hold the inner 
		// product on a single timeslice. 
		// 
		// We will use Sftmom to partition the lattice into 
		// a set of timeslices. 
		// Let phases be an instance of Sftmom, and let inner_prod
		// hold the inner products of the fields phi and chi 
		// evaluated on a set of subsets defined by phases.
		// 
		//
		// multi1d<DComplex> inner_prod = 
		// sumMulti(localInnerProduct,phi,chi),
		// 					phases.getSet());
		//
		// Then our function could be 
		// template<typename T>
		// void partitionedInnerProduct(const T& phi, const T& chi, 
		// 															multi1d<DComplex>& inner_prod,
		// 															Set& product_set){
		//
		//	inner_prod = sumMulti(localInnerProduct(phi,chi),product_set);
		// }
		//
		// Note: I'm not sure whether we have to pass inner_prod by reference 
		// 	
		// 
		//
		//
		//
		//
		//
			





	  //BRANDON: CODE GOES HERE

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
