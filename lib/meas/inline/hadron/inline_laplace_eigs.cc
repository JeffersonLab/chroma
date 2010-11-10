// $Id: inline_laplace_eigs.cc,v 1.15 2009-07-08 21:44:23 jbulava Exp $
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
#include "util/ferm/map_obj/map_obj_aggregate_w.h"
#include "util/ferm/map_obj/map_obj_factory_w.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"
#include "actions/boson/operator/klein_gord.h"
#include <qdp-lapack.h>

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
    
    template<typename T> 
    void partitionedInnerProduct(const T& phi, const T& chi, multi1d<DComplex>& inner_prod, const Set& product_set)
    {
      inner_prod = sumMulti(localInnerProduct(phi,chi),product_set);
    }
    
    template<typename T>
    void laplacian(const multi1d<LatticeColorMatrix>& u, const T& psi, T& chi, int j_decay)
    {
      T temp;
      Real minus_one = -1.;
      temp = psi * minus_one;
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
      catch( std::bad_cast )  {
	QDPIO::cerr << name << ": caught dynamic cast error" << endl;
	QDP_abort(1);
      }
      catch (const string& e) {
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
      
      try  { 
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
      try {
	std::istringstream  xml_s(params.named_obj.colorvec_obj.xml);
	XMLReader MapObjReader(xml_s);
	
	// Create the entry
	TheNamedObjMap::Instance().create< Handle< QDP::MapObject<int,EVPair<LatticeColorVector> > > >(params.named_obj.colorvec_id);
	TheNamedObjMap::Instance().getData< Handle< QDP::MapObject<int,EVPair<LatticeColorVector> > > >(params.named_obj.colorvec_id) =
	  TheMapObjIntKeyColorEigenVecFactory::Instance().createObject(params.named_obj.colorvec_obj.id,
								       MapObjReader,
								       params.named_obj.colorvec_obj.path);
      }
      catch (std::bad_cast) {
	QDPIO::cerr << name << ": caught dynamic cast error" << endl;
	QDP_abort(1);
      }
      catch (const string& e) {
	
	QDPIO::cerr << name << ": error creating prop: " << e << endl;
	QDP_abort(1);
      }
      
      // Cast should be valid now
      // Cast should be valid now
      QDP::MapObject<int,EVPair<LatticeColorVector> >& color_vecs = 
	*(TheNamedObjMap::Instance().getData< Handle< QDP::MapObject<int,EVPair<LatticeColorVector> > > >(params.named_obj.colorvec_id));

      
      // The code goes here
      StopWatch swatch;
      StopWatch fossil;
      fossil.reset();
      swatch.reset();
      swatch.start();
	  
      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, params.param.decay_dir);
      
      int num_vecs = params.param.num_vecs;
      int nt = phases.numSubsets();
      multi1d<EVPair<LatticeColorVector> >  ev_pairs(num_vecs);
      for(int n=0; n < num_vecs; ++n) { 
	ev_pairs[n].eigenValue.weights.resize(nt);
      }
      
	  
      // Choose the starting eigenvectors to have identical 
      // components and unit norm. 
      // The norm is evaluated time slice by time slice
      LatticeColorVector starting_vectors;
	  
      /*
	ColorVector ones;
	pokeColor(ones, Complex(1.0), 0);
	pokeColor(ones, Complex(1.0), 1);
	pokeColor(ones, Complex(1.0), 2);
	
	starting_vectors = ones;
      */			
      
      gaussian(starting_vectors);
      
      // Norm of the eigenvectors on each time slice
      // vector_norms is declared complex, this allows 
      // us to use partitionedInnerProduct but it may be a
      // design flaw
      multi1d<DComplex> vector_norms;
      
      
      QDPIO::cout << "Normalizing starting vector" << endl;
      
      // This function gives the norms squared
      partitionedInnerProduct(starting_vectors,starting_vectors,vector_norms,phases.getSet());
      // Apply the square root to get the true norm
      // and normalise the starting vectors
      
      QDPIO::cout << "Nt = " << nt << endl;
      
      for(int t=0; t<nt; ++t) {
	//QDPIO::cout << "vector_norms[" << t << "] = " << vector_norms[t] << endl; 
	
	vector_norms[t]  = Complex(sqrt(Real(real(vector_norms[t]))));
	starting_vectors[phases.getSet()[t]] /= vector_norms[t];
      }
	  
	  
      //Build Krlov subspace
      int kdim = 3 * params.param.num_vecs;
      int j_decay = params.param.decay_dir;
      
      QDPIO::cout << "Krylov Dim = " << kdim << endl; 
      
      // beta should really be an array of Reals	
      multi1d< multi1d<DComplex> > beta(kdim-1);	
      multi1d< multi1d<DComplex> > alpha(kdim);
      
      multi1d<double*> d(nt);
      multi1d<double*> e(nt);
      multi1d<double*> z(nt);
      
      for (int t = 0 ; t < nt ; ++t) {
	d[t] = new double[kdim];
	e[t] = new double[kdim - 1];
	z[t] = new double[(kdim) * (kdim)];
      }
      
      for (int k = 0 ; k < kdim ; ++k) {
	
	alpha[k].resize(nt);
	
	if (k < kdim - 1) {
	  
	  beta[k].resize(nt);
	}
      }
	  

      multi1d<LatticeColorVector> lanczos_vectors(kdim);
      lanczos_vectors[0] = starting_vectors;
      
      // Yields alpha[0] ... alpha[kdim-2]
      // 				beta[0] ... beta[kdim-2] 
      // 				lanczos_vector[0] ... lanczos_vector[kdim-1]
      // After the last iteration compute alpha[kdim-1]
      for(int k=0; k<kdim-1; ++k) {
	
	//QDPIO::cout << "k = " << k << endl; 
	
	
	//temporary seems to be defined as a single element but is used as both an array and a single element?
	LatticeColorVector temporary;
	// Apply the spatial Laplace operator; j_decay denotes the temporal direction		
	//laplacian(u_smr,lanczos_vectors[k],temporary,j_decay); 
	chebyshev(u_smr,lanczos_vectors[k],temporary,j_decay); 
	
	if(k > 0){	
	  for(int t=0; t<nt; ++t){
	    temporary[phases.getSet()[t]] -= beta[k-1][t]*lanczos_vectors[k-1];
	  }
	}
	    
	partitionedInnerProduct(lanczos_vectors[k],temporary,alpha[k],phases.getSet());
	
	for(int t=0; t<nt; ++t){
	  //QDPIO::cout << "alpha[k][" << t << "] = " << alpha[k][t] << endl;
	  temporary[phases.getSet()[t]] -= alpha[k][t]*lanczos_vectors[k];
	}
	    
	    
	QDPIO::cout << "Reorthogonalizing" << endl;	
	multi1d<DComplex> alpha_temp(nt);
	// Reorthogonalise - this may be unnecessary
	if(k>0){
	  
	  partitionedInnerProduct(lanczos_vectors[k-1],temporary,alpha_temp,phases.getSet());
	  for(int t=0; t<nt; ++t){
	    temporary[phases.getSet()[t]] -= alpha_temp[t]*lanczos_vectors[k-1];
	  }
	}
	    
	partitionedInnerProduct(lanczos_vectors[k],temporary,alpha_temp,phases.getSet());
	
	for(int t=0; t<nt; ++t){
	  temporary[phases.getSet()[t]] -= alpha_temp[t]*lanczos_vectors[k];
	} //
	
	    
	// Global reorthogonalisation to go here?	
	// .......
	// .....	
	
	partitionedInnerProduct(temporary,temporary,beta[k],phases.getSet());
	    
	for(int t=0; t<nt; ++t) {
	  //QDPIO::cout << "beta[k][" << t << "] = " << beta[k][t] << endl;
	  beta[k][t] = Complex(sqrt(Real(real(beta[k][t]))));
	  lanczos_vectors[k+1][phases.getSet()[t]] = temporary/beta[k][t];
	  //if (k < kdim - 1)
	  d[t][k] = toDouble(Real(real(alpha[k][t])));
	  
	  //if (k < kdim - 2)
	  e[t][k] = toDouble(Real(real(beta[k][t])));
	}
	/*
	  QDPIO::cout << "Checking orthogonality of vector " << k+1 << endl;
	  for(int m = 0; m <= k; m++){
	  multi1d<DComplex> tmp(nt);
	  partitionedInnerProduct(lanczos_vectors[k+1],lanczos_vectors[m],tmp,phases.getSet());
	  
	  
	  for(int t = 0; t < nt; t++){
	  QDPIO::cout << "   t = " << t << ": " << tmp[t] << endl;
	  }
	  }
	*/
	
      }
      // Loop over k is complete, now compute alpha[kdim-1]
      
      LatticeColorVector tmp;
      
      chebyshev(u_smr, lanczos_vectors[kdim-1], tmp, j_decay);
      
      for(int t = 0; t < nt; t++){
	tmp[phases.getSet()[t]] -= beta[kdim-2][t]*lanczos_vectors[kdim-2];
      }
      
      partitionedInnerProduct(lanczos_vectors[kdim-1],tmp,alpha[kdim-1],phases.getSet());
      
      // Finally compute eigenvectors and eigenvalues
      
      //Is AL = LT, up to small corrections? 
      
      
      /*
	QDPIO::cout << "Testing AL = LT" << endl;
	
	
	for (int k = 0 ; k < kdim -1 ; ++k)
	{
	QDPIO::cout << "Row " << k << endl;
	
	LatticeColorVector al = zero;
	
	//laplacian(u_smr, lanczos_vectors[k], al, j_decay);
	chebyshev(u_smr, lanczos_vectors[k], al, j_decay);
	
	LatticeColorVector lt = zero;
	
	for (int t = 0 ; t < nt ; ++t)
	{	
	if (k != 0)
	{
	lt[ phases.getSet()[t] ] += toDouble(Real(real(beta[k-1][t])))
	* lanczos_vectors[k-1];
	}
	
	if (k != (kdim - 2))
	{
	
	lt[phases.getSet()[t]] += toDouble(Real(real(beta[k][t]))) * 
	lanczos_vectors[k+1];
	}
	
	lt[ phases.getSet()[t] ] += toDouble(Real(real(alpha[k][t]))) * 
	lanczos_vectors[k];
	} //t
	
	LatticeColorVector ldiff = al - lt; 
	
	multi1d<DComplex> ldcnt(nt); 
	partitionedInnerProduct(ldiff, ldiff, ldcnt, phases.getSet());
	
	for (int t = 0 ; t < nt ; t++)
	
	if (toDouble(Real(real(ldcnt[t]))) > 1e-5)
	QDPIO::cout << "   dcnt[" << t << "] = " << ldcnt[t] << endl; 
	
	} //k
	    
      */
	  
      //parameters for dsteqr
      char compz = 'I';
      
      double* work  = new double[2*(kdim) - 2];
      
      int info = 0;
      int ldz = kdim;
      
      
      multi1d< multi1d< multi1d<double> > > evecs(nt);
      multi1d< multi1d<double> > evals(nt);
      
      for(int t = 0; t < nt; t++) {
	
	//QDPIO::cout << "Starting QR factorization t = " << t << endl;
	fossil.reset();
	fossil.start();
	
	QDPLapack::dsteqr(&compz, &ldz, d[t], e[t], z[t], &ldz, work, &info);
	
	fossil.stop();
	
	QDPIO::cout << "LAPACK routine completed: " << fossil.getTimeInSeconds() << " sec" << endl;
	
	QDPIO::cout << "info = " << info << endl;
	
	
	evecs[t].resize(kdim);
	evals[t].resize(kdim);
	
	for (int v = 0 ; v < kdim ; ++v) {
	  evals[t][v] = d[t][v]; 
	  
	  //QDPIO::cout << "Eval[ " << v << "] = " << evals[t][v] << endl;
	  
	  evecs[t][v].resize(kdim);
	  
	  for (int n = 0 ; n < kdim ; ++n) {
	    evecs[t][v][n] = z[t][v * (kdim ) + n ];
	  }
	  
	  //Apply matrix to vector
	  multi1d<double> Av(kdim );
	  
	  double dcnt = 0;
	  
	  for (int n = 0 ; n < kdim  ; ++n) {
	    Av[n] = 0;
	    if (n != 0) {
	      Av[n] += toDouble(Real(real(beta[n-1][t]))) * 
		evecs[t][v][n-1];
	    }
	    
	    if (n != (kdim - 1)) {
	      
	      Av[n] += toDouble(Real(real(beta[n][t]))) * 
		evecs[t][v][n+1];
	    }
	    
	    Av[n] += toDouble(Real(real(alpha[n][t]))) * 
	      evecs[t][v][n];
	    
	    
	    dcnt += (Av[n] - evals[t][v]*evecs[t][v][n]) *
	      (Av[n] - evals[t][v]*evecs[t][v][n]);
	  }//n
	  
	  //QDPIO::cout << "Vector " << v << " : dcnt = " << dcnt << endl;
	  
	}//v
	    
      }//t
      
      
      //Get Eigenvectors

      QDPIO::cout << "Obtaining eigenvectors of the laplacian" << endl;
      for (int k = 0 ; k < params.param.num_vecs ; ++k) {
	LatticeColorVector vec_k = zero;
	
	//LatticeColorVector lambda_v = zero;
	    
	for (int t = 0 ; t < nt ; ++t) {
	  //QDPIO::cout << "t = " << t << endl;
	  
	  for (int n = 0 ; n < kdim  ; ++n) {
	    vec_k[phases.getSet()[t] ] += 
	      Real(evecs[t][kdim - 1 - k][n]) * lanczos_vectors[n];
	  }
	  
	  //QDPIO::cout << "Made evector" << endl;
	  
	  //lambda_v[phases.getSet()[t] ] += 
	  //	Real(evals[t][kdim - 3 - k]) * vec_k; 
	  
	  
	  //QDPIO::cout << "Eval[" << k << "] = " <<
	  
	}
	    
	ev_pairs[k].eigenVector = vec_k;
	    
	//Test if this is an eigenvector
	//	LatticeColorVector avec = zero;
	
	/*
	//laplacian(u_smr, vec_k, avec, j_decay);
	chebyshev(u_smr, vec_k, avec, j_decay);
	
	multi1d< DComplex > dcnt_arr(nt);
	
	LatticeColorVector diffs = avec - lambda_v;
	partitionedInnerProduct( diffs, diffs, dcnt_arr, phases.getSet());  
	*/
	
	//QDPIO::cout << "Testing Lap. eigvec " << k << endl;
	/*
	  for (int t = 0 ; t < nt ; ++t)
	  {
	  if (toDouble(Real(real(dcnt_arr[t]))) > 1e-5) 
	  QDPIO::cout << "dcnt[" << t << "] = " << dcnt_arr[t] 
	  << endl;
	  }
	*/
	    
	multi1d< DComplex > temp(nt);
	multi1d< DComplex > temp2(nt);
	
	LatticeColorVector bvec = zero;
	laplacian(u_smr, vec_k, bvec, j_decay);
	partitionedInnerProduct(vec_k, bvec, temp, phases.getSet());
	partitionedInnerProduct(vec_k, vec_k, temp2, phases.getSet());
	
	//Obtaining eigenvalues of the laplacian
	for(int t = 0; t < nt; t++){
	  Complex temp3 = temp[t] / temp2[t];
	  
	  evals[t][k] = -1.0 * toDouble(Real(real(temp3)));
	  
	  ev_pairs[k].eigenValue.weights[t] = 
	    Real(evals[t][k]);
	  
	  QDPIO::cout << "t = " << t << endl;
	  QDPIO::cout << "lap_evals[" << k << "] = " << evals[t][k] << endl;
	}
	
	LatticeColorVector lambda_v2 = zero;
	
	for(int t = 0; t < nt; t++){
	  
	  lambda_v2[phases.getSet()[t]] = Real(evals[t][k]) * vec_k;
	  
	}
	
	multi1d< DComplex > dcnt_arr2(nt);
	
	LatticeColorVector diffs2 = bvec + lambda_v2;
	
	partitionedInnerProduct(diffs2, diffs2, dcnt_arr2, phases.getSet());
	
	QDPIO::cout << "Testing Laplace Eigenvalue " << k << endl;
	for(int t = 0; t < nt; t++){
	  if(toDouble(Real(real(dcnt_arr2[t]))) > 1e-5)
	    QDPIO::cout << "dcnt[" << k << "] = " << dcnt_arr2[t] << endl;
	}
	
      }//k

      for(int n=0; n < num_vecs; n++) { 
	color_vecs.insert(n, ev_pairs[n]);
      }
      color_vecs.flush();
      
      //pop(xml_out);
      
      swatch.stop();
      QDPIO::cout << name << ": time for colorvec construction = "
		  << swatch.getTimeInSeconds() 
		  << " secs" << endl;      
      
      pop(xml_out);  // CreateColorVecs
      
      
      
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
	  write(record_xml, "EigenValues", ev_pairs[i].eigenValue.weights); 
	  pop(record_xml);
	}
	pop(record_xml);
	
	// Write the propagator xml info
	TheNamedObjMap::Instance().get(params.named_obj.colorvec_id).setFileXML(file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.colorvec_id).setRecordXML(record_xml);
      }
      
      
      snoop.stop();
      QDPIO::cout << name << ": total time = "
		  << snoop.getTimeInSeconds() 
		  << " secs" << endl;
      
      QDPIO::cout << name << ": ran successfully" << endl;
      
      END_CODE();
    } 
    
  }
  
} // namespace Chroma
