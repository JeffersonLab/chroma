// $Id: t_ovlap_bj.cc,v 1.8 2004-01-02 03:19:41 edwards Exp $

#include <iostream>
#include <sstream>
#include <string>

#include <cstdio>

#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include "chroma.h"
#include "state.h"
#include "actions/ferm/linop/lovlapms_w.h"
#include "actions/ferm/fermacts/zolotarev_state.h"
#include "actions/ferm/fermacts/zolotarev4d_fermact_bj_w.h"
#include "actions/ferm/linop/lovlapms_w.h"

using namespace QDP;
using namespace std;

enum BC { ANTIPERIODIC=-1,  PERIODIC =1 };
enum GaugeStartType { HOT_START = 0, COLD_START = 1, FILE_START = 2 };
enum GaugeFormat { SZIN_GAUGE_FORMAT = 0, NERSC_GAUGE_FORMAT = 1 };



// Struct for test parameters
//
typedef struct {
  multi1d<int> nrow;
  multi1d<int> boundary;
  multi1d<int> rng_seed;
  multi1d<Real> lambda;
  Real lambda_max;
  int gauge_start_type;
  int gauge_file_format;
  string gauge_filename;
  
  Real  wilson_mass;
  Real  quark_mass;
  int  approx_order;
  double approx_min;
  double approx_max;
  Real rsd_cg;
  int max_cg;
} Param_t;

// Declare routine to read the parameters
void readParams(const string& filename, Param_t& params)
{
  XMLReader reader(filename);

  try {
    // Read Params
    read(reader, "/params/lattice/nrow", params.nrow);
    read(reader, "/params/lattice/boundary", params.boundary);
    read(reader, "/params/RNG/seed", params.rng_seed);
    if ( reader.count("/params/eValues") > 0 ) { 
      read(reader, "/params/eValues/lambda", params.lambda);
      read(reader, "/params/eValues/lambdaMax", params.lambda_max);
    }
    else { 
      params.lambda.resize(0);
    }
   read(reader, "/params/Cfg/startType", params.gauge_start_type);
   if( params.gauge_start_type == FILE_START ) { 
      read(reader, "/params/Cfg/gaugeFilename", params.gauge_filename);
      read(reader, "/params/Cfg/gaugeFileFormat", params.gauge_file_format);
    }
   

   read(reader, "/params/zolotarev/approxOrder", params.approx_order);

   if( reader.count("/params/zolotarev/approxMin") == 1 ) {
     read(reader, "/params/zolotarev/approxMin", params.approx_min);
   }
   else { 
     params.approx_min  = -1;
   }

   if( reader.count("/params/zolotarev/approxMax") == 1 ) { 
     read(reader, "/params/zolotarev/approxMax", params.approx_max);
   } else {
     params.approx_max = -1;
   }

   read(reader, "/params/zolotarev/wilsonMass", params.wilson_mass);
   read(reader, "/params/zolotarev/quarkMass",  params.quark_mass);
   read(reader, "/params/zolotarev/RsdCG", params.rsd_cg);
   read(reader, "/params/zolotarev/MaxCG", params.max_cg);



  }
  catch(const string& e) { 
    throw e;
  }
}

void dumpParams(XMLWriter& writer, Param_t& params)
{
  push(writer, "params");
  push(writer, "lattice");
  write(writer, "nrow", params.nrow);
  write(writer, "boundary", params.boundary);
  pop(writer); // lattice
  push(writer, "RNG");
  write(writer, "seed", params.rng_seed);
  pop(writer); // RNG

  if( params.lambda.size() >  0 ) { 
    push(writer, "eValues");
    write(writer, "lambda", params.lambda);
    write(writer, "lambdaMax", params.lambda_max);
    pop(writer); // eValues
  }
  push(writer, "Cfg");
  write(writer, "startType", params.gauge_start_type);
  if( params.gauge_start_type == FILE_START ) { 
    write(writer, "gaugeFileFormat", params.gauge_file_format);
    write(writer, "gaugeFilename", params.gauge_filename);
  }
  pop(writer); // Cfg

  push(writer, "zolotarev");
  write(writer, "approxOrder", params.approx_order);
  if( params.approx_min > 0 ) {
    write(writer, "approxMin", params.approx_min);
  }

  if( params.approx_max > 0 ) {
    write(writer, "approxMax", params.approx_max);
  }

  write(writer, "wilsonMass", params.wilson_mass);
  write(writer, "quarkMass",  params.quark_mass);
  write(writer, "RsdCG",      params.rsd_cg);
  write(writer, "MaxCG",      params.max_cg);
  pop(writer); // zolotarev

  pop(writer); // params
}


void readEigenVecs(const multi1d<LatticeColorMatrix>& u,
		   const UnprecWilsonFermAct& S_aux,
		   const multi1d<Real>& lambda_lo, 
		   multi1d<LatticeFermion>& eigen_vec,
		   XMLWriter& xml_out)
{


  // Create a connect State
  Handle<const ConnectState>  s(new SimpleConnectState(u));
  Handle<const LinearOperator<LatticeFermion> > D_w( S_aux.linOp(s) );


  // Create Space for the eigen vecs
  eigen_vec.resize(lambda_lo.size());
  char filename[60];
  
  // Create Space for the eigenvector norms
  multi1d<Real> e_norms(lambda_lo.size());
  multi1d<Real> evec_norms(lambda_lo.size());

  for(int i = 0; i < lambda_lo.size(); i++) { 
 
    sprintf(filename, "eigenvector_%03d", i);
    string filename_string(filename);

    cout << "Reading eigenvector: " << filename_string << endl;

    readSzinFerm(eigen_vec[i], filename_string);
    evec_norms[i] = (Real)sqrt(norm2(eigen_vec[i]));
    LatticeFermion D_ev, tmp_ev, lambda_e;

    // D_ew = D ev(i)
    (*D_w)(tmp_ev, eigen_vec[i], PLUS);
    D_ev = Gamma(15)*tmp_ev;

    // Lambda_e 
    lambda_e = lambda_lo[i]*eigen_vec[i];

    D_ev -= lambda_e;
    e_norms[i] = (Real)sqrt(norm2(D_ev));        
  }    
  push(xml_out, "EigenvectorTest");
  push(xml_out, "EigenVecNorms");
  Write(xml_out, evec_norms);
  pop(xml_out);
  push(xml_out, "EigenTestNorms");
  Write(xml_out, e_norms);
  pop(xml_out);
  pop(xml_out);
}
  
int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  Param_t params;

  try { 
    readParams("./DATA", params);
  }
  catch(const string& s) { 
    QDPIO::cerr << "Caught exception " << s << endl;
    exit(1);
  }

  /* Real kappa = Real(1) / ( Real(2)*(Real(Nd) + params.wilson_mass) );
  for(int scaleit=0; scaleit < params.lambda.size(); scaleit++) {
    params.lambda[scaleit] /= ( Real(2)*kappa );
  }
  */

  // Setup the lattice
  Layout::setLattSize(params.nrow);
  Layout::create();

  // Write out the params
  XMLFileWriter xml_out("t_ovlap.xml");
  push(xml_out, "overlapTest");

  dumpParams(xml_out, params);


  // Create a FermBC
  Handle<FermBC<LatticeFermion> >  fbc(new SimpleFermBC<LatticeFermion>(params.boundary));
  
  // The Gauge Field
  multi1d<LatticeColorMatrix> u(Nd);
  
  switch ((GaugeStartType)params.gauge_start_type) { 
  case COLD_START:
    for(int j = 0; j < Nd; j++) { 
      u(j) = Real(1);
    }
    break;
  case HOT_START:
    // Hot (disordered) start
    for(int j=0; j < Nd; j++) { 
      random(u(j));
      reunit(u(j));
    }
    break;
  case FILE_START:

    // Start from File 
    switch( (GaugeFormat)params.gauge_file_format) { 
    case SZIN_GAUGE_FORMAT:
      {
	XMLReader szin_xml;
	readSzin(szin_xml, u, params.gauge_filename);
	try { 
	  push(xml_out, "GaugeInfo");
	  xml_out << szin_xml;
	  pop(xml_out);

	}
	catch(const string& e) {
	  cerr << "Error: " << e << endl;
	}
	
      }
      break;

    case NERSC_GAUGE_FORMAT:
      {
	XMLReader nersc_xml;
	readArchiv(nersc_xml, u, params.gauge_filename);

	try { 
	  push(xml_out, "GaugeInfo");
	  xml_out << nersc_xml;
	  pop(xml_out);

	}
	catch(const string& e) {
	  cerr << "Error: " << e << endl;
	}
      }
      break;

    default:
      ostringstream file_read_error;
      file_read_error << "Unknown gauge file format" << params.gauge_file_format ;
      throw file_read_error.str();
    }
    break;
  default:
    ostringstream startup_error;
    startup_error << "Unknown start type " << params.gauge_start_type <<endl;
    throw startup_error.str();
  }


  // Measure the plaquette on the gauge
  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);
  push(xml_out, "plaquette");
  Write(xml_out, w_plaq);
  Write(xml_out, s_plaq);
  Write(xml_out, t_plaq);
  Write(xml_out, link);
  pop(xml_out);

  //! Wilsoniums;
  // Put this puppy into a handle to allow Zolo to copy it around as a **BASE** class
  // WARNING: the handle now owns the data. The use of a bare S_w below is legal,
  // but just don't delete it.
  Handle<UnprecWilsonTypeFermAct<LatticeFermion> >  S_w(new UnprecWilsonFermAct(fbc, params.wilson_mass));

  Real m_q = 0.0;
  XMLBufferWriter my_writer;

  //! N order Zolo approx, with wilson action.
  Zolotarev4DFermActBj   S(fbc, S_w, 
			   params.quark_mass,
			   params.approx_order, 
			   params.rsd_cg,
			   params.max_cg,
			   my_writer);


  const ConnectState* connect_state_ptr;
  multi1d<LatticeFermion> eigen_vecs;


  // Flick on BC's  - do not do this. Now let it be down in createState
//  phfctr(u);

  if( params.lambda.size() == 0 ) { 

    // Connect State with no eigenvectors
    connect_state_ptr = S.createState(u, 
				      Real(params.approx_min), 
				      Real(params.approx_max));
  }
  else {

    // Connect State with eigenvectors
    readEigenVecs(u, dynamic_cast<UnprecWilsonFermAct&>(*S_w), params.lambda, eigen_vecs, xml_out);

    connect_state_ptr = S.createState(u, 
				      params.lambda,
				      eigen_vecs,
				      params.lambda_max);
  }


  // Stuff the pointer into a handle. Now, the handle owns the data.
  Handle<const ConnectState> connect_state(connect_state_ptr);
						     

  // Make me a linop (this callls the initialise function)
  Handle<const LinearOperator<LatticeFermion> > D_op(S.linOp(connect_state));

  LatticeFermion psi;
  gaussian(psi);
  Double n2 = norm2(psi);
  psi /= n2;

  LatticeFermion s1, s2, s3, tmp2;
  s1 = s2 = s3 = tmp2 = zero;

  (*D_op)(s1,psi,PLUS);
  (*D_op)(s2,psi,MINUS);

  (*D_op)(tmp2, psi, PLUS);
  (*D_op)(s3, tmp2, MINUS);

  s3 *= 2;
  s3 -= s1;
  s3 -= s2;

  Double circle_norm = sqrt(norm2(s3));
  cout << "Circle Norm: " << circle_norm << endl;
  xml_out << S.getWriter();
  write(xml_out, "circleNorm", circle_norm);

  // Time to bolt
  pop(xml_out);
  QDP_finalize();

  exit(0);
}
