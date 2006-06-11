// $Id: cfgtransf.cc,v 3.1 2006-06-11 06:30:33 edwards Exp $
/*! \file
 *  \brief Many-to-many gauge transformation routine
 */

#include "chroma.h"

using namespace Chroma;

//! Many-to-many gauge transformation routine
/*! \defgroup cfgtransf Tranformation routine
 *  \ingroup main
 *
 * Main program for transforming gauge formats
 */

int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);

  START_CODE();
  
  multi1d<int> nrow(Nd);
  QDPIO::cout << "Enter lattice size\n";
  QDPIO::cin >> nrow;
  
  // Setup QDP
  Layout::setLattSize(nrow);
  Layout::create();

//  XMLFileWriter xml_out("cfgtransf.xml");
  Chroma::setXMLOutputFileName("cfgtransf.xml");
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out, "cfgtransf");

  proginfo(xml_out);    // Print out basic program info

  SzinGauge_t  szin_gauge_header;  // In case we want to write out a szin
  bool szin_gauge_header_initP = false;

#if 0
  bool AnisoP;
  QDPIO::cout << "Enter Anisotropy or not [Y/N = 1/0]\n";
  QDPIO::cin >> AnisoP;

  if ( AnisoP )
  {
    Real xi_0;
    int  t_dir;

    QDPIO::cout << "Enter the bare anisitopy factor\n";
    QDPIO::cin >> xi_0;

    QDPIO::cout << "Enter the time direction, t_dir\n";
    QDPIO::cin >> t_dir;

    QDPIO::cerr << "Currently do not support anisotropy" << endl;
    exit(1);
  }
#endif

  int input_type;
  QDPIO::cout << "Enter input Gauge field type\n"
	      << "  (1) Free field\n"
	      << "  (2) Random-transformed free field\n"
	      << "  (3) Hot start (call hotst)\n"
	      << "  (4) SZIN configuration on DV\n"
	      << "  (5) SZIN configuration on FE\n"
	      << "  (6) Illinois staggered configuration on DV\n"
	      << "  (7) MILC configuration on FE\n"
	      << "  (8) Columbia configuration on FE\n"
	      << "  (9) Schroedinger classical background field\n"
	      << " (10) FORTRAN SZIN configuration on FE\n"
	      << " (11) ASCII MILC configuration on FE\n"
	      << " (12) SZIN configuration with fund to adj transf.\n"
	      << " (13) O(3) config to U(1)\n"
	      << " (14) Instanton configuration\n"
	      << " (15) QCD Archive configuration on FE\n"
	      << " (16) MIT gauge configuration on FE\n"
	      << " (17) double prec. MIT gauge configuration on FE\n"
	      << " (18) Kentucky gauge configuration FE\n"
	      << " (19) UKQCD gauge configuration FE\n"
	      << " (20) Single-precision UKQCD gauge configuration FE\n"
	      << " (21) SZIN config in QIO format\n"
	      << " (22) Uniform back ground field\n";
  
  QDPIO::cin >> input_type;
  

  // Schroedinger BC
  Real SchrPhiMult;            /* Multiplier for Schr. BC fields */
  int loop_extent;
  int decay_dir;

  if ( input_type == 9 )
  {
    QDPIO::cout << "This will be a Schroedinger style config used for c_sw measurements\n" << endl;

    QDPIO::cout << "Enter extent of loops in decay direction\n" << endl;
    QDPIO::cin >> loop_extent;;
  
    QDPIO::cout << "Enter multiplier for Schroedinger boundary fields\n";
    QDPIO::cin >> SchrPhiMult;
 
    QDPIO::cout << "Enter the boundary condition direction\n";
    QDPIO::cin >> decay_dir;
  }
  else
  {
    loop_extent = 1;
    decay_dir   = Nd-1;
    SchrPhiMult = 1;
  }


  int output_type;
  QDPIO::cout << "Enter output Gauge field type\n"
	      << "  (1) back-end SZIN (note, on most platforms back and front-end are the same)\n"
	      << "  (2) front-end SZIN\n"
	      << "  (3) MILC config on FE\n"
	      << "  (4) QCD Archive config on FE\n"
	      << "  (5) MIT gauge config on FE\n"
	      << "  (6) Kentucky config on FE\n"
	      << "  (7) SZIN config in QIO SINGLEFILE format\n"
	      << "  (8) SZIN config in QIO MULTIFILE format\n"
	      << "  (9) replicated in time dir SZIN config in szin format\n";
  QDPIO::cin >> output_type;
  
  string cfg_input_file;
  QDPIO::cout << "Enter input file name\n";
  QDPIO::cin >> cfg_input_file;

  string cfg_output_file;
  QDPIO::cout << "Enter output file name\n";
  QDPIO::cin >> cfg_output_file;
  
  bool CGaugeP;			// Flat for Complex Conjugate
  bool HGaugeP;			// Flag for Hermitian Conjugate
  bool RGaugeP;
  multi2d<Real> theta(2,Nd);		// An array of angles for each dir


  if(input_type > 1){
    QDPIO::cout << "Complex conjugate of config?\n";
    QDPIO::cin >> CGaugeP;

    QDPIO::cout << "Hermitian conjugate of config?\n";
    QDPIO::cin >> HGaugeP;

    if ( input_type > 3 )
      {
	QDPIO::cout << "Random gauge transform of config?\n";
	QDPIO::cin >> RGaugeP;
      }
  }


  bool GSmearP;
  QDPIO::cout << "APE gauge smearing?\n";
  QDPIO::cin >> GSmearP;

  int BlkMax = 100;	/* Maximum number of blocking/smearing iterations */
  Real BlkAccu = 1.0-5;	/* Blocking/smearing accuracy */
  int sm_numb;		/* 'Smearing' number */
  Real sm_fact;		/* 'Smearing' factor */
  int j_decay;

  if ( GSmearP )
  {
    BlkAccu = 1.0e-5;
    BlkMax = 100;

    QDPIO::cout << "Enter the direction of decay\n";
    QDPIO::cin >> j_decay;

    QDPIO::cout << "Enter the number of smearing\n";
    QDPIO::cin >> sm_numb;

    QDPIO::cout << "Enter the smearing factor\n";
    QDPIO::cin >> sm_fact;
  }

  bool HypSmearP;
  QDPIO::cout << "HYP gauge smearing?\n";
  QDPIO::cin >> HypSmearP;

  int GFMax;
  Real alpha1;		/* 'Hyp-mearing' parameter */
  Real alpha2;		/* 'Hyp-mearing' parameter */
  Real alpha3;		/* 'Hyp-mearing' parameter */
  if ( HypSmearP )
  {
    BlkAccu = 1.0e-5;
    BlkMax = 100;

    QDPIO::cout << "Enter the number of smearing\n";
    QDPIO::cin >> sm_numb;

    QDPIO::cout << "Enter the smearing parameters alpha1, alpha2 and alpha3\n";
    QDPIO::cin >> alpha1 >> alpha2 >> alpha3;
  }

  bool GFixP;
  QDPIO::cout << "Gauge fixing?\n";
  QDPIO::cin >> GFixP;

  bool OrlxDo;
  Real GFAccu;
  Real OrPara;
  string gauge_rotate_file;

  if ( GFixP )
  {
    QDPIO::cout << "Enter the direction of decay\n";
    QDPIO::cin >> j_decay;

    QDPIO::cout << "Enter the gauge fixing accuracy\n";
    QDPIO::cin >> GFAccu;

    QDPIO::cout << "Enter the maximum number of gauge fixing sweeps\n";
    QDPIO::cin >> GFMax;

    QDPIO::cout << "Want over-relaxation? (yes=YES)\n";
    QDPIO::cin >> OrlxDo;

    QDPIO::cout << "Enter the over-relaxtion parameter\n";
    QDPIO::cin >> OrPara;

    QDPIO::cout << "Enter gauge rotation file name\n";
    QDPIO::cin >> gauge_rotate_file;
  }

  QDPIO::cout << "I am working on it...\n";
  
  push(xml_out,"Lattis");
  write(xml_out, "Nd", Nd);
  write(xml_out, "Nc", Nc);
  write(xml_out, "nrow", nrow);
  pop(xml_out);
  
#if 0
  if ( AnisoP )
  {
    push(xml_out,"Anisotropy_parameters");
    write(xml_out, "AnisoP", AnisoP);
    write(xml_out, "xi_0", xi_0);
    write(xml_out, "t_dir", t_dir);
    pop(xml_out);
  }
#endif

#if 0
  /* Setup Schroedinger boundary fields if desired */
  if (SchrFun > 0)
  {
    SetSFbc(SchrPhiMult, SchrFermP, theta, j_decay);
    push(xml_out,"Schroed_Func_parameters");
    write(xml_out, "SchrFun", SchrFun);
    write(xml_out, "j_decay", j_decay);
    write(xml_out, "SchrPhiMult", SchrPhiMult);
    write(xml_out, "SchrFermP", SchrFermP);
    write(xml_out, "theta", theta);
    pop(xml_out);
  }
#endif


  //
  //  Have params, now find out the source for the gauge field 
  //
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_xml_in, gauge_file_xml_in, gauge_record_xml_in;
    

  switch (input_type)
  {
  case 1:
    push(xml_out,"Free_Field");
    write(xml_out, "input_type", input_type);
    pop(xml_out);
    QDPIO::cout << "Fill u with free field\n";
    u = 1;
    break;

  case 2:
    push(xml_out,"Free_Field_with_random_gauge_transformation");
    write(xml_out, "input_type", input_type);
    pop(xml_out);
    QDPIO::cout << "Fill u with random gauge transformed free field\n";
    u = 1;
    rgauge(u);
    break;

  case 3:
    push(xml_out,"Semi-Haar_measure");
    write(xml_out, "input_type", input_type);
    pop(xml_out);
    QDPIO::cout << "Fill u with semi-Haar\n";
    HotSt(u);
    break;

  case 4:
  case 5:
    push(xml_out,"SZIN_configuration");
    write(xml_out, "input_type", input_type);
    write(xml_out, "cfg_input_file", cfg_input_file);
    pop(xml_out);
    
    QDPIO::cout << "Read SZIN u from FE file " << cfg_input_file << endl;
    readSzin(gauge_xml_in, u, cfg_input_file);
    read(gauge_xml_in, "/szin", szin_gauge_header);
    szin_gauge_header_initP = true;
    break;

// case 6: ancient Illinois staggered

  case 7:
    push(xml_out,"MILC_config");
    write(xml_out, "input_type", input_type);
    write(xml_out, "cfg_input_file", cfg_input_file);
    pop(xml_out);
    QDPIO::cout << "Read MILC u from FE file\n";
    readMILC(gauge_xml_in, u, cfg_input_file);
    break;

//  case 8: Ancient columbia config

#if 0
  case 9:
  {
    push(xml_out,"Schroedinger_BC_config");
    write(xml_out, "input_type", input_type);
    pop(xml_out);

    SchrGaugeBCParams params;
    params.loop_extent = loop_extent;
    params.SchrPhiMult = SchrPhiMult;
    params.decay_dir   = decay_dir;
    SchrNonPertGaugeBC gaugebc(params);
    u = gaugebc.SFBndFld();
  }
  break;
#endif

// case 10: Ancient fortran szin_config

#if 0
    // Not sure I want this ...
  case 11:
    push(xml_out,"Ascii_MILC_config");
    write(xml_out, "input_type", input_type);
    write(xml_out, "cfg_input_file", cfg_input_file);
    pop(xml_out);
    QDPIO::cout << "Read ASCII MILC u from FE file\n";
    readascii (u, cfg_input_file);
    break;
#endif

#if 0
    // Not yet...
  case 12:
    push(xml_out,"Fund_to_adj_SZIN_config");
    write(xml_out, "input_type", input_type);
    write(xml_out, "cfg_input_file", cfg_input_file);
    pop(xml_out);
    QDPIO::cout << "Read fund to adj SZIN from file " << cfg_input_file << endl;
    readFunToAdj (u, cfg_input_file);
    break;
#endif

#if 0
    // Not yet...
  case 13:
    push(xml_out,"O3_to_U1_config");
    write(xml_out, "input_type", input_type);
    write(xml_out, "cfg_input_file", cfg_input_file);
    pop(xml_out);
    QDPIO::cout << "Read O(3) to U(1) SZIN from file " << cfg_input_file << endl;
    ReadO3toU1 (u, cfg_input_file);
    break;
#endif

#if 0
  case 14:
  {
    multi1d<Real> center(Nd);
    Real rho;
    int su2_index;
    int sign;

    if (Nc == 1)
    {
      QDPIO::cerr << "Instanton construction requires Nc>1";
      QDP_abort(1);
    }

    QDPIO::cout << "Enter instanton center coordinates\n";
    QDPIO::cin >> center;
    QDPIO::cout << "Enter instanton size\n";
    QDPIO::cin >> rho;
    QDPIO::cout << "Enter instanton sign (+/-1)\n";
    QDPIO::cin >> sign;
    if (Nc > 2)
    {
      int j = Nc*(Nc-1)/2 - 1;
      QDPIO::cout << "Enter SU(2) subgroup index, 0 .. %d\n",j;
      QDPIO::cin >> su2_index;
    }
    else
      su2_index = 0;

    push(xml_out,"Instanton_config");
    write(xml_out, "input_type", input_type);
    write(xml_out, "center", center);
    write(xml_out, "rho", rho);
    write(xml_out, "su2_index", su2_index);
    write(xml_out, "sign", sign);
    pop(xml_out);
    QDPIO::cout << "Create instanton configuration\n";
    instanton (u, center, rho, su2_index, sign);

  break;
#endif

  case 15:
    push(xml_out,"QCD_Archive_config");
    write(xml_out, "input_type", input_type);
    pop(xml_out);
    QDPIO::cout << "Read QCD Archive u from FE file\n";
    readArchiv(gauge_xml_in, u, cfg_input_file);
    break;

#if 0
    // Someday...
  case 16:
    push(xml_out,"MIT_config");
    write(xml_out, "input_type", input_type);
    pop(xml_out);
    QDPIO::cout << "Read MIT gauge config from FE file\n";
    readmitu (cfg_input_file, u);
    break;
#endif

#if 0
    // Not now...
  case 17:
    push(xml_out,"MIT_double_config");
    write(xml_out, "input_type", input_type);
    pop(xml_out);
    QDPIO::cout << "Read double prec. MIT gauge config from FE file\n";
    readmitdu (cfg_input_file, u);
    savern (seed_old);
    break;
#endif

  case 18:
    push(xml_out,"Kentucky_config");
    write(xml_out, "input_type", input_type);
    pop(xml_out);
    QDPIO::cout << "Read Kentucky gauge config from FE file\n";
    readKYU(u, cfg_input_file);
    break;

#if 0
    // Not now...
  case 19:
    push(xml_out,"UKQCD_config");
    write(xml_out, "input_type", input_type);
    pop(xml_out);
    QDPIO::cout << "Read UKQCD gauge config from FE file\n";
    readukqcd (cfg_input_file, u);
    savern (seed_old);
    break;
#endif

#if 0
    // Not now...
  case 20:
    push(xml_out,"UKQCD_single_config");
    write(xml_out, "input_type", input_type);
    pop(xml_out);
    QDPIO::cout << "Read single-precision UKQCD gauge config from FE file\n";
    readukqcdsingle (cfg_input_file, u);
    savern (seed_old);
    break;
#endif

  case 21:
    push(xml_out,"SZIN_QIO_configuration");
    write(xml_out, "input_type", input_type);
    write(xml_out, "cfg_input_file", cfg_input_file);
    pop(xml_out);
    
    QDPIO::cout << "Read SZIN u from QIO file " << cfg_input_file << endl;
    readGauge(gauge_file_xml_in, gauge_xml_in, u, cfg_input_file,
	      QDPIO_SERIAL);
    read(gauge_xml_in, "/szin", szin_gauge_header);
    szin_gauge_header_initP = true;
    break;

  case 22:
    // Here we have a uniform, diagonal background field
    // with each direction written as 
    // U = diag(exp(i t_1), exp(i t_2), exp(-i(t_1 + t_2))

    for(int mu = 0; mu < Nd; mu++){
      QDPIO::cout << "Enter angles for direction " << mu << endl;
      QDPIO::cin >> theta(0, mu) >> theta(1, mu);
    }

    push(xml_out,"Const_diag_gauge");
    write(xml_out, "input_type", input_type);
    pop(xml_out);
    QDPIO::cout << "Creating constant diagonal gauge\n";
    constgauge(u, theta);
    break;

  default:
    QDP_error_exit("unknown input type", input_type);
  }
    

  // Dump a copy of the input gauge xml
  write(xml_out, "input_gauge_header", gauge_xml_in);
  

  // So what's the plaquette?
  MesPlq(xml_out, "Observables", u);
  xml_out.flush();
  
  if ( RGaugeP )
  {
    rgauge (u);
    MesPlq(xml_out, "Rand_Gtransf_observables", u);
  }

  if(CGaugeP){
    conjgauge(u);
    MesPlq(xml_out, "Conj_GTansf_observables",u);
  }

  if(HGaugeP){
    LatticeColorMatrix u_tmp;
    for(int mu = 0; mu < Nd; mu++){

      multi1d<int> posn(Nd);
      posn = (0,0,0,0);
      ColorMatrix uin;
      ColorMatrix uout;
      u_tmp = adj(u[mu]);

      uin = peekSite(u[mu], posn);
      uout = peekSite(u_tmp, posn);

      QDPIO::cout << "MU IS " << mu << endl << endl;
      Complex element_in, element_out;
      for(int i = 0; i < Nc; i++)
	for(int j = 0; j < Nc; j++){

	  element_in = peekColor(uin, i, j);
	  element_out = peekColor(uout, i ,j);
	  QDPIO::cout << "(i,j)= " << i << j <<
	    ", U is " << element_in << ",U dagger is "
		      << element_out << endl;
	}
      u[mu] = u_tmp;
    
    }

    MesPlq(xml_out, "Herm_Gtransf_observables", u);
  }
    

  if ( GSmearP )
  {
    multi1d<LatticeColorMatrix> u_tmp(Nd);
    if ( j_decay < Nd )
      u_tmp[j_decay] = u[j_decay];

    for(int i = 0; i < sm_numb; ++i)
    {
      for(int mu = 0; mu < Nd; ++mu)
	if ( mu != j_decay )
	  APE_Smear(u, u_tmp[mu], mu, 0, sm_fact, BlkAccu, BlkMax, j_decay);
      
      u = u_tmp;
    }
    
    push(xml_out,"Smearing_parameters");
    write(xml_out, "sm_numb",sm_numb);
    write(xml_out, "sm_fact", sm_fact);
    write(xml_out, "j_decay", j_decay);
    pop(xml_out);

    MesPlq(xml_out, "Smeared_observables", u);
  }

  if ( HypSmearP )
  {
    multi1d<LatticeColorMatrix> u_tmp(Nd);
    for(int i = 0; i < sm_numb; ++i)
    {
      Hyp_Smear(u, u_tmp, alpha1, alpha2, alpha3, BlkAccu, BlkMax);
      u = u_tmp;
    }
    
    push(xml_out,"HypSmearing_parameters");
    write(xml_out, "sm_numb", sm_numb);
    write(xml_out, "alpha1", alpha1);
    write(xml_out, "alpha2", alpha2);
    write(xml_out, "alpha3", alpha3);
    pop(xml_out);

    MesPlq(xml_out, "HypSmeared_observables", u);
  }

  LatticeColorMatrix g;
  int nrl_gf;
  if ( GFixP  )
  {
    coulGauge(u, g, nrl_gf, j_decay, GFAccu, GFMax, OrlxDo, OrPara);

    push(xml_out,"Gauge_fixing_parameters");
    write(xml_out, "j_decay", j_decay);
    write(xml_out, "GFAccu", GFAccu);
    write(xml_out, "GFMax", GFMax);
    write(xml_out, "OrlxDo",OrlxDo);
    write(xml_out, "OrPara", OrPara);
    pop(xml_out);

    MesPlq(xml_out, "Gauge_fixed_observables", u);
    
    push(xml_out,"Relaxation_iterations_in_GFIX");
    write(xml_out, "nrl_gf", nrl_gf);
    pop(xml_out);
  }


  // Compute the plaquette again
  Double w_plaq, s_plaq, t_plaq, link;
  MesPlq (u, w_plaq, s_plaq, t_plaq, link);
  
  xml_out.flush();
  
  // Make a new szin header if desired
  switch (output_type)
  {
  case 1:
  case 2:
  case 7:
  case 8:
  case 9:
  {
    bool new_headerP;
    QDPIO::cout << "Enter new szin header?" << endl;
    QDPIO::cin >> new_headerP;

    if ( new_headerP )
    {
      QDPIO::cout << "Enter TotalTrj\n";
      QDPIO::cin >> szin_gauge_header.TotalTrj;
      QDPIO::cout << "Enter NOver\n";
      QDPIO::cin >> szin_gauge_header.NOver;
      QDPIO::cout << "Enter BetaMC\n";
      QDPIO::cin >> szin_gauge_header.BetaMC;
      QDPIO::cout << "Enter BetaMD\n";
      QDPIO::cin >> szin_gauge_header.BetaMD;
      QDPIO::cout << "Enter KappaMC\n";
      QDPIO::cin >> szin_gauge_header.KappaMC;
      QDPIO::cout << "Enter KappaMD\n";
      QDPIO::cin >> szin_gauge_header.KappaMD;
      QDPIO::cout << "Enter dt\n";
      QDPIO::cin >> szin_gauge_header.dt;
      QDPIO::cout << "Enter MesTrj\n";
      QDPIO::cin >> szin_gauge_header.MesTrj;
      QDPIO::cout << "Enter Nf\n";
      QDPIO::cin >> szin_gauge_header.Nf;
      QDPIO::cout << "Enter Npf\n";
      QDPIO::cin >> szin_gauge_header.Npf;
      QDPIO::cout << "Enter RefMomTrj\n";
      QDPIO::cin >> szin_gauge_header.RefMomTrj;
      QDPIO::cout << "Enter RefFnoiseTrj\n";
      QDPIO::cin >> szin_gauge_header.RefFnoiseTrj;
      QDPIO::cout << "Enter LamPl\n";
      QDPIO::cin >> szin_gauge_header.LamPl;
      QDPIO::cout << "Enter LamMi\n";
      QDPIO::cin >> szin_gauge_header.LamMi;
      QDPIO::cout << "Enter AlpLog\n";
      QDPIO::cin >> szin_gauge_header.AlpLog;
      QDPIO::cout << "Enter AlpExp\n";
      QDPIO::cin >> szin_gauge_header.AlpExp;
      QDPIO::cout << "Enter seed\n";
      QDPIO::cin >> szin_gauge_header.seed;
    }
  }
  break;
  }


  /* Now write parameters to file cfg_output_file */
  switch (output_type)
  {
  case 1:
  case 2:
  {
    writeSzin(szin_gauge_header, u, cfg_output_file);
  }
  break;

  case 3:
  {
    /* Write a MILC format file on FE */
    MILCGauge_t milc_out;
    writeMILC (milc_out, u, cfg_output_file);
  }
  break;

  case 4:
  {
    /* Write a QCD Archive format file on FE */
    ArchivGauge_t arc_out;
    arc_out.w_plaq = w_plaq;
    arc_out.link   = link;
    writeArchiv(arc_out, u, cfg_output_file);
  }
  break;

#if 0
    // Not yet...
  case 5:
    /* Write a MIT gauge format file on FE */
    wrtmitu (cfg_output_file, u);
    break;
#endif

  case 6:
    /* Write a Kentucky gauge format file on FE */
    writeKYU(u, cfg_output_file);
    break;

  case 7:
  case 8:
  {
    QDP_volfmt_t volfmt = (output_type == 7) ? QDPIO_SINGLEFILE : QDPIO_MULTIFILE;
    XMLBufferWriter gauge_file_xml_out, gauge_record_xml_out;

    push(gauge_file_xml_out, "gauge");
    write(gauge_file_xml_out, "id", int(0));   // need something better here
    pop(gauge_file_xml_out);
    write(gauge_record_xml_out, "szin", szin_gauge_header);
    writeGauge(gauge_file_xml_out, gauge_record_xml_out, u, cfg_output_file,
	       volfmt, QDPIO_SERIAL);

    // Save the gauge rotation matrices if gauge fixed
    if ( GFixP  )
    {
      XMLBufferWriter transf_file_xml_out, transf_record_xml_out;

      push(transf_file_xml_out, "gauge_transf");
      write(transf_file_xml_out, "id", int(0));   // need something better here
      pop(transf_file_xml_out);
      
      push(transf_record_xml_out, "szin_transf");
      push(transf_record_xml_out,"Gauge_fixing_parameters");
      write(transf_record_xml_out, "j_decay", j_decay);
      write(transf_record_xml_out, "GFAccu", GFAccu);
      write(transf_record_xml_out, "GFMax", GFMax);
      write(transf_record_xml_out, "OrlxDo",OrlxDo);
      write(transf_record_xml_out, "OrPara", OrPara);
      pop(transf_record_xml_out);
      write(transf_record_xml_out, "szin", szin_gauge_header);
      pop(transf_record_xml_out);

      QDPFileWriter to(transf_file_xml_out,gauge_rotate_file,volfmt,QDPIO_SERIAL,QDPIO_OPEN);

      LatticeColorMatrixF g_f = g;
      write(to,transf_record_xml_out,g_f);         // Always save in single precision!
      close(to);
    }
  }
  break;

  case 9:
  {
    int n_replica;
    QDPIO::cout << "Enter the boundary condition direction\n";
    QDPIO::cin >> j_decay;

    QDPIO::cout << "Number of times to replicat\n";
    QDPIO::cin >> n_replica;

    writeSzinReplica(szin_gauge_header, u, j_decay,
		     n_replica,
		     cfg_output_file);
  }
  break;

  default:
    QDP_error_exit("unknown output type", output_type);
  }

  pop(xml_out);
        
  END_CODE();
  
  // Time to bolt
  Chroma::finalize();

  exit(0);
}
