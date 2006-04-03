// $Id: inline_ritz_H_w.cc,v 3.0 2006-04-03 04:59:01 edwards Exp $
/*! \file
 * \brief Inline construction of eigenvalues (Ritz)
 *
 * Eigenvalue calculations
 */

#include "fermact.h"
#include "meas/inline/eig/inline_ritz_H_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermacts/fermacts_aggregate_w.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "util/ferm/eigeninfo.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/eig/eig_spec.h"
#include "actions/ferm/linop/lopscl.h"
#include "io/eigen_io.h"

#include "meas/inline/io/default_gauge_field.h"

namespace Chroma 
{ 
  namespace InlineRitzEnv 
  { 
    AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					    const std::string& path) 
    {
      return new InlineRitz(InlineRitzParams(xml_in, path));
    }

    bool registerAll()
    {
      bool foo = true;
      foo &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
      foo &= WilsonTypeFermActsEnv::registered;
      return foo;
    }

    const std::string name = "RITZ_KS_HERM_WILSON";
    const bool registered = registerAll();
  };


  //! Eigeninfo input
  void read(XMLReader& xml, const string& path, InlineRitzParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    input.gauge_id = InlineDefaultGaugeField::readGaugeId(inputtop, "gauge_id");
    read(inputtop, "eigen_id", input.eigen_id);
  }

  //! Eigeninfo output
  void write(XMLWriter& xml, const string& path, const InlineRitzParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "eigen_id", input.eigen_id);

    pop(xml);
  }


  // Param stuff
  InlineRitzParams::InlineRitzParams()
  { 
    frequency = 0; 
    named_obj.gauge_id = InlineDefaultGaugeField::getId();
  }

  InlineRitzParams::InlineRitzParams(XMLReader& xml_in, const std::string& path) 
  {
    try 
    {
      XMLReader paramtop(xml_in, path);



      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;

//       // Parameters for source construction
//       read(paramtop, "Param", param);

      read(paramtop, "Param/version", version);

      XMLReader ferm_act_reader(paramtop, "Param/FermionAction");
      std::ostringstream fermact_os;
      ferm_act_reader.print(fermact_os);
      fermact = fermact_os.str();
      QDPIO::cout << "Read Ferm Act" << fermact << endl;


      read(paramtop, "Param/RitzParams", ritz_params);

      
      // Read any auxiliary state information
      if( paramtop.count("StateInfo") == 1 ) {
	XMLReader xml_state_info(paramtop, "StateInfo");
	std::ostringstream os;
	xml_state_info.print(os);
	stateInfo = os.str();
      }
      else { 
	XMLBufferWriter s_i_xml;
	push(s_i_xml, "StateInfo");
	pop(s_i_xml);
	stateInfo = s_i_xml.printCurrentContext();
      }

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


  void
  InlineRitzParams::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    QDP::write(xml_out, "Frequency", frequency);
    
    push(xml_out, "Param");
    QDP::write(xml_out, "version", version);

    std::istringstream  fa_is(fermact);
    XMLReader fa_reader(fa_is);
    xml_out << fa_reader;

    Chroma::write(xml_out, "RitzParams", ritz_params);
    pop(xml_out); // Param
      

    istringstream header_is(stateInfo);
    XMLReader xml_header(header_is);
    xml_out << xml_header;
    
    Chroma::write(xml_out, "NamedObject", named_obj);
   
    pop(xml_out); //  Path
  }


  void RitzCode4DHw(Handle< LinearOperator<LatticeFermion> >& MM,
		    Handle< LinearOperator<LatticeFermion> >& H,
		    const RitzParams_t& params,
		    XMLWriter& xml_out,
		    EigenInfo& eigenvec_val);


  // Function call
  void 
  InlineRitz::operator()(unsigned long update_no,
			 XMLWriter& xml_out) 
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
      string xml_file = makeXMLFileName(params.xml_file, update_no);


      push(xml_out, "RitzEigenHw");
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
  InlineRitz::func(unsigned long update_no,
		   XMLWriter& xml_out) 
  {
    START_CODE();

    QDP::StopWatch snoop;
    snoop.reset();
    snoop.start();

    // Grab the gauge field
    XMLBufferWriter gauge_xml;
    multi1d<LatticeColorMatrix> u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
    TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);

    push(xml_out, "RitzEigenHw");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << InlineRitzEnv::name << ": RitzEigenHw" << endl;

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.write(xml_out, "Input");

    XMLBufferWriter record_xml; 
    params.write(record_xml, "RecordXML");

    // Write out the config header
    write(xml_out, "Config_info", gauge_xml);

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 1);
    pop(xml_out);

    // Calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);


    TheNamedObjMap::Instance().create<EigenInfo>(params.named_obj.eigen_id);
    EigenInfo& eigenvec_val = 
      TheNamedObjMap::Instance().getData<EigenInfo>(params.named_obj.eigen_id);


    // File XML - the name of the measurement and and ID
    XMLBufferWriter file_xml;
    push(file_xml, "RitzEigenHw");
    int id = 0;    // NEED TO FIX THIS - SOMETHING NON-TRIVIAL NEEDED
    write(file_xml, "id", id);
    pop(file_xml);


    TheNamedObjMap::Instance().get(params.named_obj.eigen_id).setFileXML(file_xml);
    TheNamedObjMap::Instance().get(params.named_obj.eigen_id).setRecordXML(record_xml);
    //
    // Initialize fermion action
    //
    std::istringstream  xml_s(params.fermact);
    XMLReader  fermacttop(xml_s);
    const string fermact_path = "/FermionAction";
    string fermact;

    try
    {
      read(fermacttop, fermact_path + "/FermAct", fermact);
    }
    catch (const std::string& e) 
    {
      QDPIO::cerr << "Error reading fermact: " << e << endl;
      QDP_abort(1);
    }

    QDPIO::cout << "FermAct = " << fermact << endl;

    // Make a reader for the stateInfo
    std::istringstream state_info_is(params.stateInfo);
    XMLReader state_info_xml(state_info_is);
    string state_info_path="/StateInfo";

    bool success = false;

    if (! success)
    {

      try { 
	StopWatch swatch;
	swatch.reset();
	QDPIO::cout << "Try the various factories" << endl;

	// Typedefs to save typing
	typedef LatticeFermion               T;
	typedef multi1d<LatticeColorMatrix>  P;
	typedef multi1d<LatticeColorMatrix>  Q;

	Handle< WilsonTypeFermAct<T,P,Q> >
	  S_f(TheWilsonTypeFermActFactory::Instance().createObject(fermact,
								   fermacttop,
								   fermact_path));
	
	Handle< FermState<T,P,Q> > state(S_f->createState(u,
							  state_info_xml,
							  state_info_path));

	Handle< LinearOperator<LatticeFermion> > MM(S_f->lMdagM(state));

	Handle< LinearOperator<LatticeFermion> > H(S_f->hermitianLinOp(state));
	swatch.start();
	RitzCode4DHw(MM, H, params.ritz_params, xml_out, eigenvec_val);
	swatch.stop();
	QDPIO::cout << "Eigenvalues/-vectors computed: time= " 
		    << swatch.getTimeInSeconds() 
		    << " secs" << endl;
      
	success = true;
      }

      catch (const std::string& e) 
      {
	QDPIO::cout << InlineRitzEnv::name << ": caught exception around quarkprop: " << e << endl;
      }
    }

    if (! success)
    {
      QDPIO::cerr << "Error: no fermact found" << endl;
      QDP_abort(1);
    }
 
    snoop.stop();
    QDPIO::cout << InlineRitzEnv::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << InlineRitzEnv::name << ": ran successfully" << endl;

    pop(xml_out);

    END_CODE();
  }





  void RitzCode4DHw(Handle< LinearOperator<LatticeFermion> >& MM,
		    Handle< LinearOperator<LatticeFermion> >& H,
		    const RitzParams_t& params,
		    XMLWriter& xml_out,
		    EigenInfo& eigenvec_val)
  {

    // Try and get lowest eigenvalue of MM
    const OrderedSubset& s = MM->subset();
  
    multi1d<Real> lambda(params.Neig+params.Ndummy);
    multi1d<Real> check_norm(params.Neig);
    multi1d<LatticeFermion> psi(params.Neig
				+params.Ndummy);
      
  
    for(int i =0; i < params.Neig + params.Ndummy; i++){
      psi[i] = zero;
      gaussian(psi[i],s);
      lambda[i] = Real(1);
    }

    
    int n_CG_count;
    Real delta_cycle = Real(1);
    XMLBufferWriter eig_spec_xml;
    int n_KS_count = 0;
    int n_jacob_count = 0;
    EigSpecRitzKS(*MM, 
		  lambda, 
		  psi, 
		  params.Neig,
		  params.Ndummy,                // No of dummies
		  params.Nrenorm, 
		  params.MinKSIter, 
		  params.MaxKSIter,             // Max iters / KS cycle
		  params.MaxKS,            // Max no of KS cycles
		  params.GammaFactor,       // Gamma factor
		  params.MaxCG,
		  params.RsdR,
		  params.RsdA,  
		  params.RsdZero,
		  params.ProjApsiP,
		  n_CG_count,
		  n_KS_count,
		  n_jacob_count,
		  eig_spec_xml);
  
    // Dump output
    xml_out << eig_spec_xml;
    write(xml_out, "lambda_Msq", lambda); 
  
    // Check norms
    for(int i=0; i < params.Neig; i++) { 
      LatticeFermion Me;
      LatticeFermion lambda_e;
      (*MM)(Me, psi[i], PLUS);
      lambda_e[s] = lambda[i]*psi[i];
    
    
      LatticeFermion r_norm;
      r_norm[s] = Me - lambda_e;
    
      check_norm[i] = norm2(r_norm,s);
      check_norm[i] = sqrt(check_norm[i]);
    }
    write(xml_out, "check_norm", check_norm);
  
    for(int i=0; i < params.Neig; i++) {
      check_norm[i] /= fabs(lambda[i]);
    }
    write(xml_out, "check_norm_rel", check_norm);
  
    // Fix to ev-s of gamma_5 wilson...
    // Try to get one:
    multi1d<bool> valid_eig(params.Neig);
    int n_valid;
    int n_jacob;
  
    fixMMev2Mev(*H, 
		lambda, 
		psi, 
		params.Neig, 
		params.RsdR,
		params.RsdA, 
		params.RsdZero, 
		valid_eig, 
		n_valid, 
		n_jacob);
  
    push(xml_out, "eigFix");
    write(xml_out, "lambda_Hw", lambda);
    write(xml_out, "n_valid", n_valid);
    write(xml_out, "valid_eig", valid_eig);
  
    for(int i=0; i < params.Neig; i++) { 
      LatticeFermion Me;
      (*H)(Me, psi[i], PLUS);
    
      bool zeroP = toBool( fabs(lambda[i]) <  params.RsdZero );
      if( zeroP ) {
	check_norm[i] = norm2(Me,s);
	check_norm[i] = sqrt(check_norm[i]);
      }
      else {
	LatticeFermion lambda_e;
	LatticeFermion r_norm;
      
	lambda_e[s] = lambda[i]*psi[i];
	r_norm[s] = Me - lambda_e;
      
      
	check_norm[i] = norm2(r_norm,s);
	check_norm[i] = sqrt(check_norm[i]);
      }
   
      QDPIO::cout << "lambda_lo[" << i << "] = " << lambda[i] << "  "; 
      QDPIO::cout << "check_norm["<<i<<"] = " << check_norm[i] << endl;
    }
    write(xml_out, "check_norm", check_norm);
  
    for(int i=0; i < params.Neig; i++) { 
      check_norm[i] /= fabs(lambda[i]);
      QDPIO::cout << "check_norm_rel["<< i <<"] = " << check_norm[i] << endl;
    }
    QDPIO::cout << flush ;
    write(xml_out, "check_norm_rel", check_norm);
    pop(xml_out);
  
  
    // Now get the absolute value of the  highest e-value
    // Work with H^{dag}H = M^{dag}M
    Real hi_RsdR = 1.0e-4;
    Real hi_RsdA = 1.0e-4;
  
    multi1d<Real> lambda_high_aux(1);
    multi1d<LatticeFermion> lambda_high_vec(1);

    gaussian(lambda_high_vec[0],s);

    lambda_high_vec[0][s] /= sqrt(norm2(lambda_high_vec[0],s));

    int n_cg_high;
    XMLBufferWriter high_xml;
  
    Handle< LinearOperator<LatticeFermion> > MinusMM = new lopscl<LatticeFermion, Real>(MM, Real(-1.0));
    // Initial guess -- upper bound on spectrum
    lambda_high_aux[0] = Real(8);
  
  
    push(high_xml, "LambdaHighRitz");
      
    // Minus MM ought to produce a negative e-value
    // since MM is herm_pos_def
    // ie minus MM is hermitian -ve definite

    EigSpecRitzCG( *MinusMM,
		   lambda_high_aux,
		   lambda_high_vec,
		   1,
		   params.Nrenorm,
		   params.MinKSIter,
		   params.MaxCG,
		   hi_RsdR,
		   hi_RsdA,
		   params.RsdZero,
		   params.ProjApsiP,
		   n_cg_high,
		   high_xml);

    QDPIO::cout << "Got Here" << endl << flush ;

    lambda_high_aux[0] = sqrt(fabs(lambda_high_aux[0]));
    QDPIO::cout << "|| lambda_hi || = " << lambda_high_aux[0]  << " hi_Rsd_r = " << hi_RsdR << endl;
  
    pop(high_xml);
    xml_out << high_xml;
  
    push(xml_out, "Highest");
    write(xml_out, "lambda_hi", lambda_high_aux[0]);
    pop(xml_out);

    eigenvec_val.getEvalues().resize(params.Neig);
    eigenvec_val.getEvectors().resize(params.Neig);
    for (int i=0; i<params.Neig; i++)
      eigenvec_val.getEvalues()[i]=lambda[i];
    eigenvec_val.getLargest()=lambda_high_aux[0];
    for (int i=0; i<params.Neig; i++)
      eigenvec_val.getEvectors()[i]=psi[i];


//   QDPIO::cout << "Writing low eigenvalues/vectors" << endl;
//   writeEigen(input, lambda, psi, lambda_high_aux[0], QDPIO_SERIAL);

  }

}
