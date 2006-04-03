// $Id: t_read_eigen.cc,v 3.0 2006-04-03 04:59:16 edwards Exp $

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

#include <cstdio>

#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include "chroma.h"
#include "io/param_io.h"
#include "io/eigen_io.h"


using namespace Chroma;


struct ReadEigen_t {
  int            version;
  Real           Mass;
  multi1d<int>   boundary;
  multi1d<int>   nrow;
  QDP::Seed      seed; 
  int            Neig;
  Cfg_t          cfg;
  EigenIO_t      eigen_io_params;
};

void read(XMLReader& xml, const string& path, ReadEigen_t& param)
{
  XMLReader paramtop(xml, path);
  read(paramtop, "Param/version", param.version);
  read(paramtop, "Param/Mass",    param.Mass);
  read(paramtop, "Param/boundary", param.boundary);
  read(paramtop, "Param/nrow",     param.nrow);
  read(paramtop, "Param/rng",     param.seed);
  read(paramtop, "Param/Neig",    param.Neig);
  read(paramtop, "Cfg",       param.cfg);
  read(paramtop, "Eigen",      param.eigen_io_params);
}

int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  ReadEigen_t input;
  XMLReader xml_in(Chroma::getXMLInputFileName());

  try { 
    read(xml_in, "/ReadEigen", input);
  }
  catch( const string& e ) { 
    QDPIO::cerr << "Caught Exception: " << e << endl;
    QDP_error_exit("Exiting\n");
  }

  // Setup the lattice
  Layout::setLattSize(input.nrow);
  Layout::create();


  QDP::RNG::setrn(input.seed);

  QDPIO::cout << "ReadEigen" << endl;

  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_file_xml, gauge_xml;

  switch (input.cfg.cfg_type) 
  {
  case CFG_TYPE_SZIN :
    QDPIO::cout << "Reading SZIN Gauge config" << endl;
    readSzin(gauge_xml, u, input.cfg.cfg_file);
    break;

  case CFG_TYPE_SZINQIO:
    QDPIO::cout << "Reading SZIN QIO gauge config" << endl;
    readGauge(gauge_file_xml, gauge_xml, u, input.cfg.cfg_file, QDPIO_SERIAL);
    break;

  case CFG_TYPE_NERSC:
    QDPIO::cout << "Reading NERSC gauge config" << endl;
    readArchiv(gauge_xml, u, input.cfg.cfg_file);
    break;
  default :
    QDP_error_exit("Configuration type is unsupported.");
  }

  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "ReadEigen");

  //  write((XMLWriter &)xml_out, "InputParams", input);

  // Measure the plaquette on the gauge
  MesPlq(xml_out, "Observables", u);
  xml_out.flush();

  //! Wilsoniums;
  // Put this puppy into a handle to allow Zolo to copy it around as a **BASE** class
  // WARNING: the handle now owns the data. The use of a bare S_w below is legal,
  // but just don't delete it.
  // Create a FermBC
  Handle<FermBC<LatticeFermion> >  fbc(new SimpleFermBC<LatticeFermion>(input.boundary));
  UnprecWilsonFermAct  S_w(fbc, input.Mass);

  Handle< const ConnectState > connect_state = S_w.createState(u);

  multi1d<Real> lambda_lo(input.Neig);
  Real lambda_hi;
  multi1d<LatticeFermion> eigv_lo(input.Neig);

  ChromaWilsonRitz_t header;
  readEigen(header, lambda_lo, eigv_lo, lambda_hi, 
	    input.eigen_io_params.eigen_file, 
	    input.Neig,QDPIO_SERIAL);

  write(xml_out, "lambda", lambda_lo);
  write(xml_out, "lambda_hi", lambda_hi);
  
  Handle< const LinearOperator<LatticeFermion> > H = S_w.hermitianLinOp(connect_state);  

  multi1d<Double> check_norm(header.ritz_params.Neig);

  for(int i=0; i < header.ritz_params.Neig; i++) { 
    LatticeFermion Me;
    (*H)(Me, eigv_lo[i], PLUS);

    LatticeFermion lambda_e;
   
    lambda_e = lambda_lo[i]*eigv_lo[i];
    LatticeFermion r_norm = Me - lambda_e;
    check_norm[i] = sqrt(norm2(r_norm))/fabs(lambda_lo[i]);
  }

  write(xml_out, "check_norm", check_norm);
  pop(xml_out);

  Chroma::finalize();
    
  exit(0);
}
