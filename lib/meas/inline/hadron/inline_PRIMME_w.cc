// Use PRIMME from Andreas to find the low mode eigenspace of MdagM. - Arjun

#include "fermact.h"
#include "meas/inline/hadron/inline_PRIMME_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/smear/gaus_smear.h"
#include "meas/glue/mesplq.h"
#include "util/ferm/subset_vectors.h"
#include "util/ferm/map_obj/map_obj_aggregate_w.h"
#include "util/ferm/map_obj/map_obj_factory_w.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"
#include "actions/boson/operator/klein_gord.h"
#include "meas/inline/io/named_objmap.h"
//Added new includes below.
#include <complex> 
#include <primme.h>
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/linop/unprec_clover_linop_w.h"
#include "util/info/unique_id.h"
//#include "io/enum_io/enum_qdpvolfmt_io.h"
#ifdef BUILD_QOP_MG
#include "actions/ferm/invert/qop_mg/syssolver_mdagm_qop_mg_w.h"
#endif


namespace Chroma 
{ 
  namespace InlinePrimmeEnv 
  {

    //! input
    void read(XMLReader& xml, const std::string& path, Params::NamedObject_t& input)
    {
      XMLReader inputtop(xml, path);

      read(inputtop, "gauge_id", input.gauge_id);
    }

    //! output
    void write(XMLWriter& xml, const std::string& path, const Params::NamedObject_t& input)
    {
      push(xml, path);

      write(xml, "gauge_id", input.gauge_id);

      pop(xml);
    }

    //! input
    void read(XMLReader& xml, const std::string& path, Params::Param_t& input)
    {
      XMLReader inputtop(xml, path);

      //Fixed code to include all primme params.

      //multi1d<int> boundary;

      //read(inputtop, "boundary", boundary);

      XMLReader xml_tmp(inputtop, "FermionAction");

      XMLBufferWriter xml_out;
      push(xml_out,"FermionAction");
      xml_out << xml_tmp;
      //write(xml_out, "boundary", boundary);
      pop(xml_out);

      XMLReader xml_inn(xml_out);
      input.fermact = readXMLGroup(xml_inn, "/FermionAction", "FermAct");

      read(inputtop, "Nev", input.nev);
      read(inputtop, "Tol", input.tol);
      read(inputtop, "PrintLevel", input.printLevel);
      read(inputtop, "LambdaC", input.lambda_C);
      read(inputtop, "LambdaMax", input.lambda_max);
      read(inputtop, "NCheb", input.n_cheb);
      read(inputtop, "file_name", input.file_name);
      read(inputtop, "maxBasisSize", input.maxBasisSize);
      read(inputtop, "primme_method", input.primme_method);
#ifdef BUILD_QOP_MG
      read(inputtop, "MG_Params", input.invParam);
#endif
    }

    //! output
    void write(XMLWriter& xml, const std::string& path, const Params::Param_t& out)
    {
      push(xml, path);

      //Fixed code to include all primme params.
      write(xml, "Nev" , out.nev);
      write(xml, "Tol", out.tol);
      write(xml, "LambdaC", out.lambda_C);
      write(xml, "LambdaMax", out.lambda_max);
      write(xml, "NCheb", out.n_cheb);
      write(xml, "file_name", out.file_name);
      write(xml, "maxBasisSize", out.maxBasisSize);
      write(xml, "primme_method", out.primme_method);
#ifdef BUILD_QOP_MG
      write(xml, "MG_Params", out.invParam);
#endif

      pop(xml);
    }

    //! input
    void read(XMLReader& xml, const std::string& path, Params& input)
    {
      Params tmp(xml, path);
      input = tmp;
    }

    //! output
    void write(XMLWriter& xml, const std::string& path, const Params& input)
    {
      push(xml, path);
    
      write(xml, "Param", input.param);
      write(xml, "NamedObject", input.named_obj);

      pop(xml);
    }
  } // namespace InlineCreateColorVecsEnv 


  namespace InlinePrimmeEnv 
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
      
    const std::string name = "PRIMME";

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

	// Read in the output source configuration info
	read(paramtop, "NamedObject", named_obj);

	// Possible alternate XML file pattern
	if (paramtop.count("xml_file") != 0) 
	{
	  read(paramtop, "xml_file", xml_file);
	}
      }
      catch(const std::string& e) 
      {
	QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
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
	std::string xml_file = makeXMLFileName(params.xml_file, update_no);

	push(xml_out, "PRIMME");
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

    //Holds the input parameters 
    struct EigenInput_t
    {
  	int    nev;          // Number of eigenvalues to be calculated
  	double tol;          // Relatice accuracy of eigenvalues and vectors 
  	int    printLevel;   // Debug

  	double lambda_C;     // cutoff of D^t*D
  	double lambda_max;   // max eig of D^t*D
  	int    n_cheb;       // Rank of T_n(x)
    };

    //A full Handle to an unpreconditioned Dirac operator is used. For that, these types are needed.
    typedef LatticeFermion               T;    typedef multi1d<LatticeColorMatrix>  P;    typedef multi1d<LatticeColorMatrix>  Q;

    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    //------------------------------------------------------------------------------------
    // Shifted q(x) = 1 + 2*(x + lambda_C)/(lambda_max - lambda_C)
    void applyRemap(/*Handle<UnprecLinearOperator<T,P,Q> >*/Handle<LinearOperator<T> > linop,
		    const LatticeFermion& psi,
		    LatticeFermion& chi,
		    double lambda_C,
		    double lambda_max)
    {
	LatticeFermion vec_tmp;
	LatticeFermion Dv;
	(*linop)(Dv, psi, PLUS);
  	(*linop)(vec_tmp, Dv, MINUS);
  
  	Real ftmp1 = 2.0 / (lambda_max - lambda_C);
  	Real ftmp2 = 1.0 + ftmp1*lambda_C;

  	// Shifted q(x) = 1 + 2*(x + lambda_C)/(lambda_max + lambda_C)
  	chi = ftmp2*psi + ftmp1*vec_tmp;
    }

    //------------------------------------------------------------------------------------
    // Routine for T[n, q]
    void applyChebyshev(/*Handle<UnprecLinearOperator<T,P,Q> >*/Handle<LinearOperator<T> > linop,
		        const LatticeFermion& psi,
		        LatticeFermion& chi,
		        double lambda_C,
		        double lambda_max,
		        int n)
    {
  	if (n < -1)
  	{
    		QDP_error_exit("buggered T[n,x]");
  	}
  	else if (n == -1)
  	{
		LatticeFermion Dv;
		(*linop)(Dv, psi, PLUS);
    		(*linop)(chi, Dv, MINUS);
  	}
  	else if (n == 0)
  	{
    		chi = psi;
  	}
  	else if (n == 1)
  	{
    		applyRemap(linop, psi, chi, lambda_C, lambda_max);
  	}
  	else
  	{
    		LatticeFermion tmp1, tmp2, tmp3;
    		applyChebyshev(linop, psi, tmp2, lambda_C, lambda_max, n-2);
    		applyChebyshev(linop, psi, tmp3, lambda_C, lambda_max, n-1);
    		applyRemap(linop, tmp3, tmp1, lambda_C, lambda_max);

    		chi = Real(2)*tmp1 - tmp2;
  	}
    }


    struct LoadedPrimmeParams_t : public primme_params
    {
  	LoadedPrimmeParams_t(const primme_params& p, 
		             /*Handle<UnprecLinearOperator<T,P,Q> >*/Handle<LinearOperator<T> > linop_,
		       	     double C_, double max_, int n_ 
#ifdef BUILD_QOP_MG 
,Handle< MdagMSysSolverQOPMG> MdagMinv
#endif
) : 
    	primme_params(p), linop(linop_), lambda_C(C_), lambda_max(max_), n_cheb(n_) 
#ifdef BUILD_QOP_MG 
,MdagMinv(MdagMinv)
#endif
{}

  	// Hold the other useful stuff
  	/*Handle<UnprecLinearOperator<T,P,Q> >*/Handle<LinearOperator<T> >  linop;  // The operator to apply
#ifdef BUILD_QOP_MG
	Handle< MdagMSysSolverQOPMG> MdagMinv; //Multigrid inverter to be used as a preconditioner.
#endif

  	double lambda_C;     // cutoff of -lap
  	double lambda_max;   // max eig of -lap
  	int    n_cheb;       // Rank of T_n(x)
    };

    /******************************************************************************
     * Applies the matrix vector multiplication on a block of vectors.
     ******************************************************************************/
    void MatrixMatvec(void *x, void *y, int *blockSize, primme_params *primme) 
    {
    //  QDPIO::cout << __PRETTY_FUNCTION__ << ": entering" << std::endl;

    	LoadedPrimmeParams_t* params = static_cast<LoadedPrimmeParams_t*>(primme);

      	typedef LatticeFermion::Subtype_t Type_t;
	
	//QDPIO::cout<<"The block size here is: "<<*blockSize<<std::endl;
  	for(int i = 0; i < *blockSize; i++) 
  	{
    		Complex_Z* xvec = reinterpret_cast<Complex_Z*>(x);
    		Complex_Z* yvec = reinterpret_cast<Complex_Z*>(y);

    		Type_t* x_p = reinterpret_cast<Type_t*>(xvec + primme->nLocal*i);
    		Type_t* y_p = reinterpret_cast<Type_t*>(yvec + primme->nLocal*i);

		//QDPIO::cout<<"i is: "<<i<<std::endl;		
		//QDPIO::cout<<"nLocal is: "<<primme->nLocal<<std::endl;
    		LatticeFermion cvec_x(x_p, 0.0);
    		LatticeFermion cvec_y(y_p, 0.0);

    		applyChebyshev(params->linop, cvec_x, cvec_y, params->lambda_C, params->lambda_max, params->n_cheb);
  	}

	//  QDPIO::cout << __PRETTY_FUNCTION__ << ": exiting" << std::endl;
    }

#ifdef BUILD_QOP_MG
    void DoMG( Handle< MdagMSysSolverQOPMG> MdagMinv,
		        const LatticeFermion& psi,
		        LatticeFermion& chi)
    {
	(*MdagMinv)(chi, psi);
    }
    void MultiGridPreconditioner(void *x, void *y, int *blockSize, primme_params *primme) 
    {
    //  QDPIO::cout << __PRETTY_FUNCTION__ << ": entering" << std::endl;

    	LoadedPrimmeParams_t* params = static_cast<LoadedPrimmeParams_t*>(primme);

      	typedef LatticeFermion::Subtype_t Type_t;
	
  	for(int i = 0; i < *blockSize; i++) 
  	{
    		Complex_Z* xvec = reinterpret_cast<Complex_Z*>(x);
    		Complex_Z* yvec = reinterpret_cast<Complex_Z*>(y);

    		Type_t* x_p = reinterpret_cast<Type_t*>(xvec + primme->nLocal*i);
    		Type_t* y_p = reinterpret_cast<Type_t*>(yvec + primme->nLocal*i);

    		LatticeFermion cvec_x(x_p, 0.0);
    		LatticeFermion cvec_y(y_p, 0.0);

    		DoMG(params->MdagMinv, cvec_x, cvec_y);
  	}
    }
#endif

    // This is a HACK
    #define BLAS_DCOPY  dcopy_
    extern "C" {
    	void BLAS_DCOPY(int *n, double *x, int *incx, double *y, int *incy);
    }

    #ifdef USE_QMP
    void globalSumDouble(void *sendBuf, void *recvBuf, int *count, primme_params *params) {
   	int ONE = 1;
    	QMP_status_t f;

    	BLAS_DCOPY(count, (double *) sendBuf, &ONE, (double *) recvBuf, &ONE);
    	f = QMP_sum_double_array((double* ) recvBuf, *count);
    }
    /*#else
    void globalSumDouble(void *sendBuf, void *recvBuf, int *count, primme_params *params) 
    {
   	int ONE = 1;

   	BLAS_DCOPY(count, (double *) sendBuf, &ONE, (double *) recvBuf, &ONE);

    }*/
    #endif
    
    //*****************************************************************************
    // Routine used to perform the inversion
    // evecs : list of nev eigenvectors 
    // evals : list of eigenvalues
    // nev   : number of eigenvalues/vectors desired
    void getEigenPairs(/*Handle<UnprecLinearOperator<T,P,Q> >*/Handle<LinearOperator<T> > linop,
		       multi1d< multi1d< std::complex<double> > >& evecs, 
		       multi1d<double> & evals,
		       const EigenInput_t &param,
		       const int maxBasisSize,
           	       const std::string primme_method,
		       multi1d<LatticeFermion>& final_evecs
#ifdef BUILD_QOP_MG 
,Handle< MdagMSysSolverQOPMG> MdagMinv
#endif
)
    {
  	int nev = param.nev;  // target number of desired evs

  	// Size of evec
	int n = QDP::Layout::vol()*Nc*Ns;
	int nLocal = QDP::Layout::sitesOnNode()*Nc*Ns;
	QDPIO::cout<<"The evecs are of size "<<n<<std::endl;
	QDPIO::cout<<"The evecs are of local size "<<nLocal<<std::endl;

  	QDPIO::cout << "Nev = " << nev << "  Tol=" << param.tol << "  Ncheb= " << param.n_cheb << "  lambda_C= " << param.lambda_C << "  lambda_max= " << param.lambda_max << std::endl;

  	if ( nev > n) 
  	{
    		std::cerr << "ERROR: nev must be <= the dimension of the matrix" << std::endl;
    		exit(1);

  	}

  	// Primme params
  	primme_params primme;

 	// Primme solver setup
  	// Initialize defaults in primme
  	primme_initialize(&primme);
	
  	// Preset method
  	primme_preset_method method;
	
	//Method is an xml input.
  	//method = DYNAMIC;
  	primme.numEvals = nev;
  	primme.printLevel = param.printLevel;
  	primme.n = n;
	primme.nLocal = nLocal;
  	primme.aNorm = 1.0;
  	primme.eps = param.tol;
  	primme.target = primme_smallest;
	primme.procID = Layout::nodeNumber();
	//Needed for parallelism.

  	// Correction parameters
  	primme.correctionParams.precondition = 0;
	
  	//Preconditioner either 0 or multgrid.
  	primme.matrixMatvec = MatrixMatvec;
#ifdef BUILD_QOP_MG
	primme.applyPreconditioner = MultiGridPreconditioner;
	primme.correctionParams.precondition = 1;
	QDPIO::cout<<"MultiGrid Preconditioner set!"<<std::endl;
#else
  	primme.applyPreconditioner = 0;
#endif
	
	//Some parscalar calls are set here.
	#ifdef USE_QMP
	QDPIO::cout<<"Parallel version initializing."<<std::endl;
	primme.globalSumDouble = globalSumDouble;
	primme.numProcs = n/nLocal;
	#endif

	//Section below based on user defined inputs.
	if(maxBasisSize != 0)
		primme.maxBasisSize = maxBasisSize;
	if(primme_method == "DYNAMIC")	
		method = DYNAMIC;
	else if(primme_method == "DEFAULT_MIN_TIME")
		method = DEFAULT_MIN_TIME;
	else if(primme_method == "DEFAULT_MIN_MATVECS")
		method = DEFAULT_MIN_MATVECS;
	else if(primme_method == "Arnoldi")
		method = Arnoldi;
	else if(primme_method == "GD")
		method = GD;
	else if(primme_method == "GD_plusK")
		method = GD_plusK;
	else if(primme_method == "GD_Olsen_plusK")
		method = GD_Olsen_plusK;
	else if(primme_method == "JD_Olsen_plusK")
		method = JD_Olsen_plusK;
	else if(primme_method == "RQI")
		method = RQI; 
	else if(primme_method == "JDQR")
		method = JDQR;
	else if(primme_method == "JDQMR")
		method = JDQMR;
	else if(primme_method == "JDQMR_ETol")
		method = JDQMR_ETol;
	else if(primme_method == "SUBSPACE_ITERATION")
		method = SUBSPACE_ITERATION;
	else if(primme_method == "LOBPCG_OrthoBasis")
		method = LOBPCG_OrthoBasis;
	else if(primme_method == "LOBPCG_OrthoBasis_Window")
		method = LOBPCG_OrthoBasis_Window;
	else
	{
		QDPIO::cout<<"The primme method you are trying to call doesn't actually exist..."<<std::endl;
		QDPIO::cout<<"Have a nice day."<<std::endl;
		QDPIO::cout<<"PEACING OUT!"<<std::endl;
		QDP_abort(1);
	} 

  	// Should set lots of defaults
  	if (primme_set_method(method, &primme) < 0 ) 
  	{
		QDPIO::cout<<"Error number is "<<primme_set_method(method, &primme)<<std::endl;
    		QDPIO::cerr << __func__ << ": invalid preset method\n";
    		//QDP_abort(1);
  	}

  	// Do not use
  	primme.matrix         = 0;
  	primme.preconditioner = 0;

  	// Optional: report memory requirements
  	int ret = zprimme(NULL,NULL,NULL,&primme);
  	QDPIO::cout << "return code = " << ret << "\n";

  	QDPIO::cout << "\nPRIMME will allocate the following memory:\n";
  	QDPIO::cout << "real workspace, " << primme.realWorkSize << " bytes\n";
  	QDPIO::cout << "int  workspace, " << primme.intWorkSize << " bytes\n";

  	// Print params that PRIMME defaulted
	if (Layout::primaryNode())
  		primme_display_params(primme);
	//Only head node should be printing params if multiple processors running.

  	// Allocate space for converged Ritz values and residual norms
  	double* primme_evals = (double *)primme_calloc(primme.numEvals, sizeof(double), "evals");
	//nLocal used instead of n below.
  	Complex_Z* primme_evecs = (Complex_Z*)primme_valloc(primme.nLocal*
					       (primme.numEvals+primme.maxBlockSize)*sizeof(Complex_Z), 
					       "evecs");
  	double* primme_rnorms = (double *)primme_calloc(primme.numEvals, sizeof(double), "rnorms");
  

  	// Initial guess (optional)
	//nLocal replaced n in loop below, this matches Andreas's code.
  	for(int i = 0; i < primme.nLocal; i++) 
  	{
    		//primme_evecs[i].r = 1.0/sqrt(primme.n);
		//Create a random initial guess vector using a QDP function call.
		Double x;
		random(x);
		double d = reinterpret_cast<double &>(x);
		primme_evecs[i].r = d;
    		primme_evecs[i].i = 0.0L;
  	}

  	// Make a bigger structure holding 
  	LoadedPrimmeParams_t big_primme_params(primme, 
					       linop, 
					       param.lambda_C, param.lambda_max, param.n_cheb
#ifdef BUILD_QOP_MG 
,MdagMinv
#endif
); 


  	// Call primme
  	QDPIO::cout << __func__ << ": Calling zprimme" << std::endl;

  	StopWatch swatch;
  	swatch.reset(); swatch.start();

  	//  double wt1 = primme_get_wtime(); 
	//Only give the smaller, proper primme struct.
  	ret = zprimme(primme_evals, primme_evecs, primme_rnorms, &big_primme_params);
  	//  double wt2 = primme_get_wtime(); 

  	swatch.stop();
  	QDPIO::cout << "zprimme routine completed: " << swatch.getTimeInSeconds() << " sec" << std::endl;

  	// Reporting
  	primme_PrintStackTrace(primme);

  	if (primme.procID == 0) 
  	{
    		for (int i=0; i < primme.numEvals; i++) 
    		{
      			fprintf(stdout, "Eval[%d]: %-22.15g rnorm: %-22.15g\n", i+1, abs(primme_evals[i]), primme_rnorms[i]); 
    		}
		//Since the bigger struct is fed into primme, it is also used for printing.
    		fprintf(stdout, " %d eigenpairs converged\n", /*primme*/big_primme_params.initSize);

    		fprintf(stdout, "Tolerance : %-22.15E\n", /*primme*/big_primme_params.aNorm*primme.eps);
    		fprintf(stdout, "Iterations: %-d\n", /*primme*/big_primme_params.stats.numOuterIterations); 
    		fprintf(stdout, "Restarts  : %-d\n", /*primme*/big_primme_params.stats.numRestarts);
    		fprintf(stdout, "Matvecs   : %-d\n", /*primme*/big_primme_params.stats.numMatvecs);
    		fprintf(stdout, "Preconds  : %-d\n", /*primme*/big_primme_params.stats.numPreconds);
    
    		//    fprintf(stdout, "\n\n#,%d,%.1f\n\n", primme.stats.numMatvecs, wt2-wt1); 
    		fprintf(stdout, "\n\n#,%d\n\n", /*primme*/big_primme_params.stats.numMatvecs); 

		QDPIO::cout<<"DynamicMethodSwitch is: "<</*primme*/big_primme_params.dynamicMethodSwitch<<std::endl; 

    		switch (/*primme*/big_primme_params.dynamicMethodSwitch) 
    		{
    		case -1: fprintf(stdout,
			"Recommended method for next run: DEFAULT_MIN_MATVECS\n"); break;
    		case -2: fprintf(stdout,
		     	"Recommended method for next run: DEFAULT_MIN_TIME\n"); break;
    		case -3: fprintf(stdout,
			"Recommended method for next run: DYNAMIC (close call)\n"); break;
    		}
  	}

	if (ret != 0) 
      	{
    		QDPIO::cerr << "Error: zprimme returned with nonzero exit status\n"<<std::endl;
		QDPIO::cout<<"ret is equal to "<<ret<<std::endl;		
    		QDP_abort(1);
  	}

  	QDPIO::cout << __func__ << ": clean up" << std::endl;

  	// Cleanup
  	primme_Free(&primme);
	//Bigger struct also freed.
	primme_Free(&big_primme_params);

  	QDPIO::cout << __func__ << ": reorganize eigs\n" << std::endl;

  	evals.resize(nev);
  	evecs.resize(nev);

  	for(int val = 0 ; val < nev ; ++val)
  	{
    		LatticeFermion temp;
    		std::complex<double>* temp_p = reinterpret_cast<std::complex<double>*>(temp.getF());

		//nLocal used instead of n below for the same reasons as loop.
    		evecs[val].resize(nLocal);

		//nLocal instead of n used in loop below, this only matters for parscalar.
    		for (int m = 0 ; m < nLocal ; ++m)
    		{
      			Complex_Z zz = primme_evecs[m + primme.nLocal*val];
      			evecs[val][m] = std::complex<double>(zz.r, zz.i);
      			temp_p[m]     = evecs[val][m];
    		}

		//Back door way of storing the eigenvectors when running with many processors.
		final_evecs[val] = temp;

    		LatticeFermion temp2;
		
		LatticeFermion Dv;
		(*linop)(Dv, temp, PLUS);    		
		(*linop)(temp2, Dv, MINUS);

    		Double   val1 = innerProductReal(temp, temp2);
    		Double   val2 = innerProductReal(temp, temp);

    		Double   val3 = val1 / val2;

    		evals[val] = toDouble(real(val3));
	
    		QDPIO::cout << "mdagm_evals[" << val << "] = " << evals[val] << std::endl;
  		}
  	QDPIO::cout << "\n" << std::endl;
    }


    //*****************************************************************************
    //Test if the calculated values/vectors actually 
    //satisfy Ax = lambda * x
    void test( /*Handle<UnprecLinearOperator<T,P,Q> >*/Handle<LinearOperator<T> > linop,
	       const multi1d< multi1d< std::complex<double> > >& eigen_vecs,
	       const multi1d<double>& eigen_vals, double tol)
    {
  	QDPIO::cout << "Testing converged eigenvalues and eigenvectors" << std::endl;

  	int nev = eigen_vecs.size();

  	for (int i = 0 ; i < nev ; ++i)
  	{
    		LatticeFermion x = zero;
    		std::complex<double>* x_p = reinterpret_cast< std::complex<double>* >(x.getF());
    		for (int d = 0 ; d < eigen_vecs[i].size() ; ++d)
    		{
      			x_p[d] = eigen_vecs[i][d];
    		}

    		Double lambda = eigen_vals[i];

		//This function is not called for now.
    		LatticeFermion lambda_x = lambda * x;
    		LatticeFermion Ax = zero;		

		LatticeFermion Dv;
		(*linop)(Dv, x, PLUS);
    		(*linop)(Ax, Dv, MINUS);

    		Real dcnt = norm2(Ax - lambda_x);
    		QDPIO::cout << "evec_num = " << i << "  |diff|^2 = " << dcnt << std::endl;

    		if ( toBool( (Real(tol) < dcnt) ) )
    		{
      			QDPIO::cerr << "ERROR: diff > 10 * tol, INVALID EIGEN PAIR" << std::endl;
      			//exit(0);
    		}
  	}
    }

    //Test the eigenvalue equation again, after everything has been converted to the proper QDP types.
    void test2( /*Handle<UnprecLinearOperator<T,P,Q> >*/Handle<LinearOperator<T> > linop,
	       multi1d<LatticeFermion>& eigen_vecs,
	       multi1d<Real>& eigen_vals, double tol)
    {
  	QDPIO::cout << "Testing converged eigenvalues and eigenvectors" << std::endl;

  	int nev = eigen_vecs.size();

  	for (int i = 0 ; i < nev ; ++i)
  	{
    		LatticeFermion x = eigen_vecs[i];

    		Double lambda = eigen_vals[i];
		QDPIO::cout<<"lambda "<<i<<": "<<lambda<<std::endl;

    		LatticeFermion lambda_x = lambda * x;
    		LatticeFermion Ax;		

		LatticeFermion Dv;
		(*linop)(Dv, x, PLUS);
    		(*linop)(Ax, Dv, MINUS);

    		Real dcnt = norm2(Ax - lambda_x);
    		QDPIO::cout << "evec_num = " << i << "  |diff|^2 = " << dcnt << std::endl;

    		if ( toBool( (Real(tol) < dcnt) ) )
    		{
      			QDPIO::cerr << "ERROR: diff > 10 * tol, INVALID EIGEN PAIR" << std::endl;
      			//exit(0);
    		}
  	}
    }

    //*****************************************************************************
    void getEigenV(multi1d<LatticeFermion>& final_evecs,
	           multi1d<Real>& final_evals,
	           /*Handle<UnprecLinearOperator<T,P,Q> >*/Handle<LinearOperator<T> > linop,
	           const EigenInput_t& eig_param,
		   const int maxBasisSize,
		   const std::string primme_method
#ifdef BUILD_QOP_MG 
,Handle< MdagMSysSolverQOPMG> MdagMinv
#endif
)
    {
  	int nev = eig_param.nev;

  	multi1d< multi1d< std::complex<double> > > eigen_vecs;
  	multi1d<double> eigen_vals;

  	StopWatch swatch;
  	swatch.reset();

  	QDPIO::cout << "\n\nCalling EigenPairs" << std::endl; 
  	swatch.start();

	//Store eigenvalues and eigenvectors for writing
  	final_evals.resize(nev);
  	final_evecs.resize(nev);

  	// Find eigenvectors
	getEigenPairs(linop, eigen_vecs, eigen_vals, eig_param, maxBasisSize, primme_method, final_evecs
#ifdef BUILD_QOP_MG 
,MdagMinv
#endif
);

	//getEigenPairs directly takes in final_evecs so resizing must be done first.

  	swatch.stop();
  	QDPIO::cout << "EigenPairs completed: time = " << swatch.getTimeInSeconds() << " sec" << std::endl;

	//Only the test with LatticeFermions is needed, in parallel reinterpret casting again may cause issues.	  	
	//QDPIO::cout << "\n\nTesting Results..." << std::endl;
  	//swatch.reset();
  	//swatch.start();

  	//Test results
  	//test(linop, eigen_vecs, eigen_vals, eig_param.tol);
	//Ax and lambda*x need to have x be the eigenvector.

  	//swatch.stop();
  	//QDPIO::cout << "Test completed: time = " << swatch.getTimeInSeconds() << " sec" << std::endl;

  	//QDPIO::cout << "evals:" << std::endl;

  	for (int i = 0 ; i < nev ; ++i)
  	{
    		final_evals[i] = Real( std::real(eigen_vals[i]) );
  	}

	//This only works in scalar.
  	/*for (int i = 0 ; i < nev ; ++i)
  	{
    		//Put eigenvector into final LatticeFermion
    		std::complex<double>* final_p = reinterpret_cast< std::complex<double>* >(final_evecs[i].getF());
    		for (int d = 0 ; d < eigen_vecs[i].size() ; ++d)
    		{
      			final_p[d] = eigen_vecs[i][d];
    		}

  	}*/

	QDPIO::cout<<"\n\nVectors successfully copied and converted to LatticeFermions!"<<std::endl;	
	QDPIO::cout << "\n\nTesting Results for LatticeFermions..." << std::endl;
  	swatch.reset();
  	swatch.start();

  	//Test results	
  	test2(linop, final_evecs, final_evals, eig_param.tol);
	//Ax and lambda*x need to have x be the eigenvector.

  	swatch.stop();
  	QDPIO::cout << "Test 2 completed: time = " << swatch.getTimeInSeconds() << " sec" << std::endl;
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
  
	//Only a Dirac operator with Clover fermions is needed for now.
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
		QDPIO::cerr << name << ": caught dynamic cast error" << std::endl;
		QDP_abort(1);
        }
        catch (const std::string& e) 
        {
		QDPIO::cerr << name << ": map call failed: " << e << std::endl;
		QDP_abort(1);
        }
	//
	// Initialize fermion action
	//
	std::istringstream  xml_s(params.param.fermact.xml);
	XMLReader  fermacttop(xml_s);
	QDPIO::cout << "FermAct = " << params.param.fermact.id << std::endl;

	// Generic Wilson-Type stuff
	//Handle< FermionAction<T,P,Q> >
	  //S_f(TheFermionActionFactory::Instance().createObject(params.param.fermact.id,
							       //fermacttop,
							       //params.param.fermact.path));
	//If a more general Dirac Operator is needed, this is how to do it.
	Handle< FermAct4D<T,P,Q> >
	  S_f(TheWilsonTypeFermActFactory::Instance().createObject(params.param.fermact.id,
							       fermacttop,
							       params.param.fermact.path));
	Handle< FermState<T,P,Q> > state(S_f->createState(u));
	/*Handle<UnprecLinearOperator<T,P,Q> > linop;
	if(params.param.fermact.id == "CLOVER"){
		CloverFermActParams cp(fermacttop,params.param.fermact.path);
		linop=new UnprecCloverLinOp(state,cp) ;
		QDPIO::cout<<"CLOVER Dirac Operator has been created successfully."<<std::endl;
	}
	else
	{
		QDPIO::cout<<"This action is not currently supported, sorry..."<<std::endl;
		QDPIO::cout<<"PEACING OUT!"<<std::endl;
		QDP_abort(1);
	}*/

	Handle<LinearOperator<T> > linop=S_f->linOp(state) ;
#ifdef BUILD_QOP_MG
	//SysSolverQOPMGParams invParam;
	Handle< MdagMSysSolverQOPMG> MdagMinv;
	//MdagMSysSolverQOPMG MdagMinv(linop, state, invParam);
	MdagMinv = new MdagMSysSolverQOPMG(linop,state,params.param.invParam);
#endif

  	// Print the operator coefficients
  	std::cout << "\n\n";
  	//linop->printOp(std::cout);
	//The unpreconditioned operator used here does not have a printOp function.

  	//Final eigenvectors 
  	multi1d<LatticeFermion> final_evecs;

  	//Final eigenvalues
  	multi1d<Real> final_evals;

  	QDPIO::cout << "Calculating eigenpairs" << std::endl;
	
	EigenInput_t local_eig_param;
	local_eig_param.nev = params.param.nev;
	local_eig_param.tol = params.param.tol;
	local_eig_param.printLevel = params.param.printLevel;
	local_eig_param.lambda_C = params.param.lambda_C;
	local_eig_param.lambda_max = params.param.lambda_max;
	local_eig_param.n_cheb = params.param.n_cheb;
  	getEigenV(final_evecs, final_evals, linop, local_eig_param, params.param.maxBasisSize, params.param.primme_method
#ifdef BUILD_QOP_MG 
,MdagMinv
#endif
);
	
  	QDPIO::cout << "Writing Results" << std::endl;

  	//
  	// Write out the eigenvectors
  	//
  	QDPIO::cout << "Write out eigenvectors" << std::endl;

  	XMLBufferWriter file_xml;

  	push(file_xml, "MODMetaData");

  	// All the eigenvalues, but pack in an array
  	if (1)
  	{
    	multi1d< multi1d<Real> > temp_evals(final_evals.size());

    	// Resize first
    	for (int i = 0 ; i < final_evals.size(); i++)
    	{
      		temp_evals[i].resize(1);
      		temp_evals[i][0] = final_evals[i];
    	}

    	write(file_xml, "Weights", temp_evals);
  	}
  	pop(file_xml);

	//Fake OPTeigCG way of writing to disk.
	//This requires creation of specific data types.
	multi1d<Complex> H(final_evals.size()*final_evals.size());
      	multi1d<Complex> HU(final_evals.size()*final_evals.size());
	for(int index = 0; index < H.size(); index ++)
	{
		H(index) = 0;
		HU(index) = 0;
	}
	for(int index = 0; index < final_evals.size(); index ++)
	{
		H(index*(final_evals.size() + 1)) = final_evals(index);
		HU(index*(final_evals.size() + 1)) = sqrt(final_evals[index]);
	}
	// File XML                                            
        XMLBufferWriter OPTeigCG_file_xml;
      	push(OPTeigCG_file_xml, "OptEigInfo");
      	write(OPTeigCG_file_xml, "id", uniqueId());
      	write(OPTeigCG_file_xml, "N", final_evecs.size());
      	write(OPTeigCG_file_xml, "ncurEvals", final_evecs.size());
      	write(OPTeigCG_file_xml, "restartTol", -1);
      	write(OPTeigCG_file_xml, "lde", final_evals.size());
      	write(OPTeigCG_file_xml, "ldh", final_evals.size());
      	pop(OPTeigCG_file_xml);

	// Open file                                
        QDPFileWriter to(OPTeigCG_file_xml,
		         params.param.file_name,
		         QDPIO_SINGLEFILE,
		         QDPIO_PARALLEL/*QDPIO_SERIAL*/,QDPIO_OPEN);
	//Writing in parallel!
	
	for(int v(0);v<final_evecs.size();v++){
		XMLBufferWriter record_xml;
		push(record_xml, "EigenVector");
		write(record_xml,"no",v);
		pop(record_xml);
		write(to, record_xml, final_evecs[v]);
      	}
      
      	{
		XMLBufferWriter record_xml;
		push(record_xml, "EigenValues");
		pop(record_xml);
		write(to, record_xml, final_evals);
      	}
      	{
		XMLBufferWriter record_xml;
		push(record_xml, "H");
		pop(record_xml);
		write(to, record_xml, H);
      	}
      	{
		XMLBufferWriter record_xml;
		push(record_xml, "HU");
		pop(record_xml);
		write(to, record_xml, HU);
      	}

  	QDPIO::cout << "File Written " << std::endl;

  	snoop.stop();
  	QDPIO::cout << "Inline PRIMME Completed: time = " << snoop.getTimeInSeconds() << " sec " << std::endl;


  	pop(xml_out);     
      
      	QDPIO::cout << name << ": ran successfully" << std::endl;
      
      	END_CODE();
    } 
    
  }
  
} // namespace Chroma
