// $Id: seqprop.cc,v 1.2 2003-12-17 05:04:28 edwards Exp $
/*! \file
 *  \brief Main code for sequential propagator generation
 */

#error "CODE NOT READY YET"


#include <iostream>
#include <cstdio>

#define MAIN

#include "chroma.h"

//! Sequential propagator generation
/*
 *  \defgroup propagator Propagator generation
 *  \ingroup main
 *
 *  Read quark propagators, compute the sequential sources needed
 *  for baryon and meson form factors and/or structure function moments.
 *
 *  This routine does not compute the form factors or moments --
 *  that is done in separate routines....
 *
 */

int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  START_CODE("seqprop");

  // Input parameter structure
  Propagator_input_t  input;

  // Instantiate xml reader for DATA
  XMLReader xml_in("DATA");

  // Read data
  read(xml_in, "/propagator", input);

  // Specify lattice size, shape, etc.
  Layout::setLattSize(input.param.nrow);
  Layout::create();

  // Read in the configuration along with relevant information.
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_xml;

  switch (input.param.cfg_type) 
  {
  case CFG_TYPE_SZIN :
    readSzin(gauge_xml, u, input.cfg.cfg_file);
    break;
  default :
    QDP_error_exit("Configuration type is unsupported.");
  }


  // Read in the source along with relevant information.
  LatticePropagator quark_prop_source;
  XMLReader source_xml;

  switch (input.param.prop_type) 
  {
  case PROP_TYPE_SZIN :
//    readSzinQprop(source_xml, quark_prop_source, input.prop.source_file);
    quark_prop_source = 1;
    break;
  default :
    QDP_error_exit("Propagator type is unsupported.");
  }


  // Instantiate XML writer for XMLDAT
  XMLFileWriter xml_out("XMLDAT");
  push(xml_out, "seqprop");

  // Write out the input
  write(xml_out, "Input", xml_in);

  // Write out the config header
  write(xml_out, "Config_info", gauge_xml);

  // Write out the source header
  write(xml_out, "Source_info", source_xml);

  push(xml_out, "Output_version");
  write(xml_out, "out_version", 1);
  pop(xml_out);

  xml_out.flush();


  // Check if the gauge field configuration is unitarized
  unitarityCheck(u);

  // Calculate some gauge invariant observables just for info.
  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq(u, w_plaq, s_plaq, t_plaq, link);

  push(xml_out, "Observables");
  Write(xml_out, w_plaq);
  Write(xml_out, s_plaq);
  Write(xml_out, t_plaq);
  Write(xml_out, link);
  pop(xml_out);

  xml_out.flush();


//---------------------------------------------------------------------------------------
// MOVE THIS CODE AROUND
  LatticePropagator quark_propagator;
  LatticePropagator seq_quark_prop;
  LatticePropagator quark_prop_src;

  LatticeReal psi_sq;
  multi1d<Double> pion(length);
  multi1d<Double> quark_src(length);
  multi1d<Double> seq_quark(length);
  int length;

  LatticeBoolean lbad;	/* not used */

  multi1d<int> boundary(Nd);
  multi1d<int> t_srce(Nd);		/* Coordinates of the quark propagator source. */
  int t_sink;		/* Time coordinate of the sequential quark
				   propagator source (ie. the sink). */
  multi1d<int> sink_mom(Nd-1);	/* Sink baryon momentum */
  multi1d<int> seed(4);		/* Random number seed (see SETRN for meaning) */

  int j_decay;		/* Direction to measure propagators. */
  int cfg_type;		/* Configuration type - szin, Illinois staggered */

  Real Kappa_mes;

  int loop;
  int seq_src_ctr;
  int seq_src_value;

  int SU2AdjP;
  int NcbW;
  int IOper;
  int i;
  int j;
  int mu;
  int cb;
  int numbad;
  String cfg_file;
  String prop_file;

  int version;		/* input parameter version */
  int io_version;		/* input parameter version */
  int out_version;		/* output format version */
  int prop_io_location;	/* FE/BE Location of propagators */

  int Z3_src;		/* Z3 random value source at 2^Z3_src
			   pts/dir */
  int Pt_src;		/* Point source */
  int Pt_snk;		/* Point sink */
  int Sl_src;		/* Shell source */
  int Sl_snk;		/* Shell sink */
  int src_type;		/* Source type */
  int snk_type;		/* Sink type */
  int Wvf_kind;		/* Wave function kind: gauge invariant */
  int wvf_type;		/* Wave function type */
  multi1d<Real> Kappa(numKappa);		/* Array of Kappa's passed to hadron */
  int numKappa;		/* Number of Kappa's passed to hadron */
  multi1d<Real> wvf_param(numKappa);	/* Array of width's or other parameters
					   for "shell" source/sink wave function. */
  multi1d<int> WvfIntPar(numKappa);	/* Array of iter numbers to approx. Gaussian */

  /* or terminate CG inversion for Wuppertal smearing */

  multi1d<Real> RsdCG(numKappa);	        /* CG accuracy */

  int OperProj;            /* Operator used for projection */
  Real Kappa_t;
  int NWilsVec;            /* Number of vectors for projection */
  Real lam_hi;                 /* Largest eigenvalue of OperProj */

  int ncg_had;		/* Number of CG iterations in hadron. */

  LatticeReal lrftmp;		/* Junk slice to fill with random numbers */

  Double w_plaq;		/* Gauge invariant plaquette variables... */
  Double s_plaq;
  Double t_plaq;
  Double link;			/* Gauge-variant link variable */
  multi1d<DComplex> pollp(Nd);		/* Polyakov loop */

  /*
   *  Now we have some information that specifies the particular sequential sources
   *  we are computing
   */

  int numSeq_src;		/* The total number of sequential sources */
  multi1d<int> Seq_src(numSeq_src);	/* An array containing a list of the 
				   sequential source */


  Complex pion_src;		/* Exponentiated propagator
				   evaluated to the source */


  Complex Z3_wght;	/* Not sure this is needed. We keep Z3_src only
			   because of propagator file naming convention! */
  char propagator_file[100];

  /* Initialize dynamic memory allocator. */
  INITIALIZE_ALLOCATOR;

  io_version = 3;		/* Input parameter version */
  out_version = 3;		/* Output format version */

  SU2AdjP = 0;
  OperProj = -1;
  NWilsVec = 0;
  prop_io_location = FE_BINARY_LOCATION;   /* For now, all propagators come 
					      from the "back-end". */
  PolyArgResc = 1;
  AnisoP = NO;
  xi_0 = 1;
  xiF_0 = 1;
  Wilsr_s = 1;

  push(nml_in,"IO_version");
  Read(nml_in, version);
  pop(nml_in);
  push(nml_in,"param");
  Read(nml_in, FermTypeP);
  Read(nml_in, Nd);
  Read(nml_in, Nc);
  Read(nml_in, numKappa);
  pop(nml_in);
  global_Nd = Nd;
  global_Nc = Nc;

    
  switch (FermTypeP)
  {
  case WILSON_FERMIONS:
    FPRINTF(trm_out," seqprop: Sequential Propagators for Wilson fermions\n");
    push(nml_in,"param");
    Read(nml_in, Kappa);
    pop(nml_in);
    for(i = 0; i < numKappa; ++i)
      FPRINTF(trm_out,"Propagator  Kappa: %7.4f\n",Kappa[i]);

    proginfo ();
    push(xml_out,"IO_version");
    Write(xml_out, version);
    pop(xml_out);
    push(xml_out,"Output_version");
    Write(xml_out, out_version);
    pop(xml_out);

    break;
  case STAGGERED_FERMIONS:
    QDP_error_exit("Staggered fermions not supported");
    break;
  default:
    QDP_error_exit("Fermion type incorrect", FermTypeP);
  }

        
  switch (version)
  {
  case 1:

    /*
     *  This case is appropropriate for Wilson-type actions on
     *  an isotropic lattice.
     */

    push(nml_in,"param");
    Read(nml_in, OverMass);
    Read(nml_in, RatPolyDeg);
    Read(nml_in, PolyArgResc);
    Read(nml_in, NWilsVec);
    pop(nml_in);
    push(nml_in,"param");
    Read(nml_in, cfg_type);
    Read(nml_in, j_decay);
    pop(nml_in);
    push(nml_in,"param");
    Read(nml_in, Z3_src);
    pop(nml_in);
    push(nml_in,"param");
    Read(nml_in, Pt_src);
    Read(nml_in, Sl_src);
    pop(nml_in);
    push(nml_in,"param");
    Read(nml_in, Pt_snk);
    Read(nml_in, Sl_snk);
    Read(nml_in, t_sink);
    Read(nml_in, sink_mom);
    pop(nml_in);
    push(nml_in,"param");
    Read(nml_in, InvType);
    Read(nml_in, FermAct);
    Read(nml_in, H_parity);
    Read(nml_in, ClovCoeff);
    Read(nml_in, u0);
    Read(nml_in, MRover);
    pop(nml_in);
    push(nml_in,"param");
    Read(nml_in, MaxCG);
    Read(nml_in, RsdCG);
    pop(nml_in);
    push(nml_in,"param");
    Read(nml_in, Wvf_kind);
    pop(nml_in);
    push(nml_in,"param");
    Read(nml_in, wvf_param);
    pop(nml_in);
    push(nml_in,"param");
    Read(nml_in, WvfIntPar);
    pop(nml_in);

    numSeq_src = 4;		/* Now set up the source */

    
    Seq_src[0] = 0;
    Seq_src[1] = 1;
    Seq_src[2] = 2;
    Seq_src[3] = 3;

    break;

  case 2:
    /*
     *  This case is appropropriate for Wilson-type actions on
     *  an anisotropic lattice.  Thus we need to specify extra parameters.
     */

    push(nml_in,"param");
    Read(nml_in, OverMass);
    Read(nml_in, RatPolyDeg);
    Read(nml_in, PolyArgResc);
    Read(nml_in, NWilsVec);
    pop(nml_in);
    push(nml_in,"param");
    Read(nml_in, cfg_type);
    Read(nml_in, j_decay);
    pop(nml_in);

    push(nml_in,"param");
    Read(nml_in, AnisoP);
    Read(nml_in, t_dir);
    Read(nml_in, xi_0);
    Read(nml_in, xiF_0);
    Read(nml_in, Wilsr_s);
    pop(nml_in);
    push(nml_in,"param");
    Read(nml_in, InvType);
    Read(nml_in, FermAct);
    Read(nml_in, H_parity);
    Read(nml_in, ClovCoeff);
    Read(nml_in, ClovCoeffR);
    Read(nml_in, ClovCoeffT);
    Read(nml_in, u0);
    Read(nml_in, MRover);
    pop(nml_in); /* The aniso params */

    push(nml_in,"param");
    Read(nml_in, Z3_src);
    pop(nml_in);
    push(nml_in,"param");
    Read(nml_in, Pt_src);
    Read(nml_in, Sl_src);
    pop(nml_in);
    push(nml_in,"param");
    Read(nml_in, Pt_snk);
    Read(nml_in, Sl_snk);
    Read(nml_in, t_sink);
    Read(nml_in, sink_mom);
    pop(nml_in);
    push(nml_in,"param");
    Read(nml_in, MaxCG);
    Read(nml_in, RsdCG);
    pop(nml_in);
    push(nml_in,"param");
    Read(nml_in, Wvf_kind);
    pop(nml_in);
    push(nml_in,"param");
    Read(nml_in, wvf_param);
    pop(nml_in);
    push(nml_in,"param");
    Read(nml_in, WvfIntPar);
    pop(nml_in);

    numSeq_src = 4;		/* Now set up the source */

    
    Seq_src[0] = 0;
    Seq_src[1] = 1;
    Seq_src[2] = 2;
    Seq_src[3] = 3;

    break;

  case 3:

    /*
     *  This case is appropropriate for Wilson-type actions on
     *  an anisotropic lattice.  Thus we need to specify extra parameters.
     */

    push(nml_in,"param");
    Read(nml_in, OverMass);
    Read(nml_in, RatPolyDeg);
    Read(nml_in, PolyArgResc);
    Read(nml_in, NWilsVec);
    pop(nml_in);
    push(nml_in,"param");
    Read(nml_in, cfg_type);
    Read(nml_in, j_decay);
    pop(nml_in);

    push(nml_in,"param");
    Read(nml_in, AnisoP);
    Read(nml_in, t_dir);
    Read(nml_in, xi_0);
    Read(nml_in, xiF_0);
    Read(nml_in, Wilsr_s);
    pop(nml_in);
    push(nml_in,"param");
    Read(nml_in, InvType);
    Read(nml_in, FermAct);
    Read(nml_in, H_parity);
    Read(nml_in, ClovCoeff);
    Read(nml_in, ClovCoeffR);
    Read(nml_in, ClovCoeffT);
    Read(nml_in, u0);
    Read(nml_in, MRover);
    pop(nml_in); /* The aniso params */

    push(nml_in,"param");
    Read(nml_in, Z3_src);
    pop(nml_in);
    push(nml_in,"param");
    Read(nml_in, Pt_src);
    Read(nml_in, Sl_src);
    pop(nml_in);
    push(nml_in,"param");
    Read(nml_in, Pt_snk);
    Read(nml_in, Sl_snk);
    Read(nml_in, t_sink);
    Read(nml_in, sink_mom);
    pop(nml_in);
    push(nml_in,"param");
    Read(nml_in, MaxCG);
    Read(nml_in, RsdCG);
    pop(nml_in);
    push(nml_in,"param");
    Read(nml_in, Wvf_kind);
    pop(nml_in);
    push(nml_in,"param");
    Read(nml_in, wvf_param);
    pop(nml_in);
    push(nml_in,"param");
    Read(nml_in, WvfIntPar);
    pop(nml_in);

    /*
     *  Now we read in the information associated with the sequential sources
     */

    push(nml_in,"param");
    Read(nml_in, numSeq_src);
    pop(nml_in);
    
    /*
     *  Now read in the particular Sequential Sources we are evaluating
     */

    push(nml_in,"param");
    Read(nml_in, Seq_src);
    pop(nml_in);

    for(seq_src_ctr = 0; seq_src_ctr < numSeq_src; seq_src_ctr++)
      FPRINTF(trm_out, "Computing sequential source of type %d\n", 
	      Seq_src[seq_src_ctr]);

    break;

  default:
    QDP_error_exit("Unsupported input parameter version", version);
  }

  /* Specify lattice size, shape, etc. */
        
  push(nml_in,"param");
  Read(nml_in, nrow);
  pop(nml_in);
  push(nml_in,"param");
  Read(nml_in, boundary);
  pop(nml_in);
  push(nml_in,"param");
  Read(nml_in, t_srce);
  pop(nml_in);
  push(nml_in,"Cfg");
  Read(nml_in, cfg_file);
  pop(nml_in);

  /* Check if this is the correct binary for these fermions */
  /* Set the following, so that "chkferm" does not crash! */
  KappaMC = 0.5;
  KappaMD = 0.5;
  MassMC = 0.5;
  MassMD = 0.5;
  chkferm ();
  Ns = global_Ns;
  Ns2 = Ns/2;
  KappaMC = 0;
  KappaMD = 0;
  MassMC = 0;
  MassMD = 0;

  /* Check for valid parameter values */
  for ( i=0; i <= Nd-1; i=i+1 )
    if( INTEGER_MOD_FUNCTION(nrow[i],2) != 0 )
      QDP_error_exit("lattice shape is invalid odd linear size not allowed", i, nrow[i]);

#if 0
  for ( i=0; i <= Nd-1; i=i+1 )
    if ( INTEGER_ABS_FUNCTION(boundary[i]) != 1 )
      QDP_error_exit("boundary condition invalid must by +/- 1", i, boundary[i]);
#endif

  for ( i=0; i <= Nd-1; i=i+1 )
    if ( t_srce[i] < 0 || t_srce[i] >= nrow[i]  )
      QDP_error_exit("quark propagator source coordinate incorrect", i, t_srce[i]);

  if ( t_sink < 0 || t_sink >= nrow[j_decay] )
    QDP_error_exit("sink time coordinate incorrect", t_sink);

  FPRINTF(trm_out,"\n     Gauge group: SU(%1d)\n",Nc);

  for ( i=0; i <= numKappa-1; i=i+1 )
    if ( Kappa[i] < 0.0 )
      QDP_error_exit("unreasonable value for Kappa", i, Kappa[i]);

  for ( i=0; i <= numKappa-1; i=i+1 )
    if ( MaxCG < 0 || RsdCG[i] < WORD_VALUE(WORD_RsdCG,SMALL) ||
	 RsdCG[i] > WORD_VALUE(WORD_RsdCG,HALF) )
      QDP_error_exit("unreasonable CG parameters", MaxCG, RsdCG[i]);


  Z3_wght = 0;

  /* Check for unnecessary multiple occurances of kappas and/or wvf_params */
  if ( numKappa > 1 )
  {
    if (Sl_src == YES)
    {
      for( i=1; i <= numKappa-1; i=i+1 )
	for( j=0; j<= i-1; j=j+1 )
	  if( Kappa[j] == Kappa[i] &&
	      wvf_param[j] == wvf_param[i] )
	    QDP_error_exit("Same kappa and wvf_param", i, Kappa[i], wvf_param[i], j, Kappa[j], wvf_param[j]);
    }
    else
    {
      for ( i=1; i <= numKappa-1; i=i+1 )
	for ( j=0; j<= i-1; j=j+1 )
	  if ( Kappa[j] == Kappa[i] )
	    QDP_error_exit("Same kappa without shell source or sink", i, Kappa[i], j, Kappa[j]);
    }
  }

  /* Set the operator used for projection if needed */
  if (NWilsVec != 0)
  {
    operproj (FermAct, OperProj);

    if ( OperProj < 0 )
      QDP_error_exit("Operator used for projection not set. Check fermact", FermAct);
  }

  /* Initialize neighbour communication and boundary conditions */
  INITIALIZE_GEOMETRY;
  INITIALIZE_BOUNDARY;

  FPRINTF(trm_out,"     volume: %2d",nrow[0]);
  for(i=1;i<Nd;i++)
    FPRINTF(trm_out," x %2d",nrow[i]);
  FPRINTF(trm_out,"\n");

  
  /* Read in the configuration along with relevant information. */
  FPRINTF(stderr,"Initialize config\n");

  switch (cfg_type)
  {
  case 1:
    readszin (cfg_file, u, seed);
    break;
  case 2:
    readst (cfg_file, u, seed);
    break;
  case 3:
    readmilc (cfg_file, u);
    push(nml_in,"Seed");
    Read(nml_in, seed);
    pop(nml_in);
    break;
  case 4:
    readfeszin (cfg_file, u, seed);
    break;
  case 5:
    push(xml_out,"Free_Field");
    Write(xml_out, cfg_type);
    pop(xml_out);
    u = 1;
    seed = 0;
    seed[0] = 11;
    break;
  case 6:
    push(xml_out,"Semi-Haar_measure");
    Write(xml_out, cfg_type);
    pop(xml_out);
    seed = 0;
    seed[0] = 11;
    setrn (seed);
    HotSt (u);
    break;
  default:
    QDP_error_exit("configuration type is invalid", cfg_type);
  }
  
  FPRINTF(stderr,"Config successfully initialized\n");


  
  push(xml_out,"param1");
  Write(xml_out, FermTypeP);
  Write(xml_out, Nd);
  Write(xml_out, Nc);
  Write(xml_out, Ns);
  Write(xml_out, numKappa);
  Write(xml_out, Kappa);
  pop(xml_out);
  push(xml_out,"param2");
  Write(xml_out, FermAct);
  Write(xml_out, OverMass);
  Write(xml_out, RatPolyDeg);
  Write(xml_out, PolyArgResc);
  Write(xml_out, NWilsVec);
  Write(xml_out, H_parity);
  Write(xml_out, ClovCoeff);
  Write(xml_out, ClovCoeffR);
  Write(xml_out, ClovCoeffT);
  Write(xml_out, u0);
  pop(xml_out);
  push(xml_out,"param3");
  Write(xml_out, cfg_type);
  Write(xml_out, j_decay);
  pop(xml_out);
  push(xml_out,"param4");
  Write(xml_out, MaxCG);
  Write(xml_out, RsdCG);
  Write(xml_out, InvType);
  Write(xml_out, MRover);
  pop(xml_out);
  push(xml_out,"param5");
  Write(xml_out, Pt_src);
  Write(xml_out, Sl_src);
  Write(xml_out, Z3_src);
  pop(xml_out);
  push(xml_out,"param6");
  Write(xml_out, Pt_snk);
  Write(xml_out, Sl_snk);
  Write(xml_out, t_sink);
  Write(xml_out, sink_mom);
  pop(xml_out);
  push(xml_out,"param7");
  Write(xml_out, Wvf_kind);
  Write(xml_out, wvf_param);
  Write(xml_out, WvfIntPar);
  pop(xml_out);
  push(xml_out,"param8");
  Write(xml_out, AnisoP);
  Write(xml_out, t_dir);
  Write(xml_out, xi_0);
  Write(xml_out, xiF_0);
  pop(xml_out);
  push(xml_out,"param9");
  Write(xml_out, numSeq_src);
  Write(xml_out, Seq_src);
  pop(xml_out);
  push(xml_out,"param10");
  Write(xml_out, seed);
  pop(xml_out);
  push(xml_out,"lattis");
  Write(xml_out, nrow);
  Write(xml_out, boundary);
  Write(xml_out, t_srce);
  pop(xml_out);
  FLUSH_WRITE_NAMELIST(xml_out);


  int ncg_had = 0;			/* Initialise iteration counter */

  /* If we require a shell wave function sink type, determine it now: */
  wvf_type = 0;
  if (Sl_snk == YES)
  {
    switch (Wvf_kind)
    {
    case 3:
      wvf_type = OPTION(GAUGE_INV_GAUSSIAN_WVF);
      break;
    case 4:
      wvf_type = OPTION(WUPPERTAL_WVF);
      break;
    default:
      QDP_error_exit("Unsupported gauge-invariant Wvf_kind[not 3 or 4]", Wvf_kind);
    }
  }

  /*
   *  Allocate the source type
   */
  if (Pt_src == YES)
  {
    src_type = OPTION(POINT_SOURCE);
    if (Sl_src == YES){
      QDPIO::cout << " Warning: shell source ignored; do point source" << endl;
    }
  }
  else if (Sl_src == YES){
    src_type = OPTION(SHELL_SOURCE);
  }
  else{
    QDP_error_exit("Must specify point source or shell source");
  }
 
  /*
   *  Allocate the sink type
   */

  if (Pt_snk == YES){
    snk_type = OPTION(POINT_SINK);
    if (Sl_snk == YES){
      FPRINTF(trm_out," Warning: shell sink ignored; do point sink\n");
    }
  }
  else if (Sl_snk == YES){
    snk_type = OPTION(SHELL_SINK);
  }
  else{
    QDP_error_exit("Must specify point sink or shell sink");
  }
 
  /*
   *  Turn on the boundary conditions through the phase factors.
   */
  setph(boundary);

  multi1d<LatticeColorMatrix> u_tmp = u;
  phfctr(u_tmp);		/* Boundary phases on */

  /*
   *  Now loop over the various kappas
   */
    
  for(loop = 0; loop < numKappa; loop++)
  {
    /*
     * Read the quark propagator
     */
    sprintf(propagator_file,"propagator_%d", loop);
    strcpy(prop_file, propagator_file);
    FPRINTF(stderr,"start read prop_file %s\n",prop_file);
    readqprop (prop_io_location, prop_file, quark_propagator);
    FPRINTF(stderr,"finished read prop_file %s\n",prop_file);

    if (Sl_snk == YES)
    {
      /*
       * Do the sink smearing BEFORE the interpolating operator
       */

      sink_smear2 (u, quark_propagator, wvf_type, wvf_param[loop], WvfIntPar[loop], j_decay);
    }

    for(seq_src_ctr = 0; seq_src_ctr < numSeq_src; seq_src_ctr++)
    {
      FPRINTF(stdout,"Start seqprop calculation for seq_src number %d\n",seq_src_ctr);

      /*
       * Allocate space for the sequential source
       */

      
      /*
       *  Sources 0 -> 9 corresponding to Baryon sequential sources
       *  Sources 10 -> 19 corresponds to a Meson sequential source
       *  Souces  21 -> 29 are additional Baryon ones we thought of
       *
       *  Note that not all the source values are necessarily implemented
       *
       */

      seq_src_value = Seq_src[seq_src_ctr]; /* Assign the particular 
					       source type */


      if(((0 <= seq_src_value) && (seq_src_value <= 9)) ||
	 ((21 <= seq_src_value) && (seq_src_value <= 29))) 
      {
	/*
	 *  Computation of the Baryon sequential source
	 */
      
	barSeqSource (quark_propagator, quark_propagator, quark_prop_src, 
		      t_sink, sink_mom, j_decay, seq_src_value);
      }
      else if ((10 <= seq_src_value) && (seq_src_value <= 20))
      {
	/*
	 * Computation of the Meson sequential source
	 */

	mesonSeqSource (quark_propagator, quark_prop_src, 
			t_sink, sink_mom, j_decay, seq_src_value);
      }
      else{
	QDP_error_exit("Unknown sequential source type", seq_src_value);
      }

      if (Sl_snk == YES)
      {

	/*
	 * Do the sink smearing AFTER the interpolating operator
	 */
	sink_smear2 (u, quark_prop_src, wvf_type, wvf_param[loop], WvfIntPar[loop], j_decay);
      }

      /*
       *  Compute the full propagator.
       */

      seq_quark_prop = 0;
      QuarkProp (u_tmp, quark_prop_src, Kappa[loop], RsdCG[loop], seq_quark_prop, ncg_had);
      
      /*
       *  Write the sequential propagator out to disk
       *
       *  We need to dump some sort of header.  At the very least, we
       *  should dump the kappa values
       */

      
      /*
       *  Now create the name of the propagator value
       */

      sprintf(propagator_file,"seqprop_%d_%d", loop, seq_src_value);
      strcpy(prop_file, propagator_file);
      Kappa_mes = Kappa[loop];
      wrtqprop_Kappa(Kappa_mes);
      wrtqprop (prop_io_location, prop_file, seq_quark_prop);
      
      /*
       *  In the case of the pion, we know that the exponentiated propagator
       *  back to the source should be the pion correlator at time-slice
       *  zero, and so will write this out
       */

      if(seq_src_value == 10)
      {
	seqPionTest(pion_src, seq_quark_prop, t_srce);
	
	push(xml_out,"Seq_propagator_test");
	Write(xml_out, pion_src);
	pop(xml_out);
      }
    } /* end loop over sequential sources */
      
  } /* end loop over the kappa value */

    
  push(xml_out,"Relaxation_Iterations");
  Write(xml_out, ncg_had);
  pop(xml_out);

  pop(xml_out);

  xml_out.close();
  xml_in.close();

  // Time to bolt
  QDP_finalize();

  END_CODE("seqprop");

  exit(0);
}

