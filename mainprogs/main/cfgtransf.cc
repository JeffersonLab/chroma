// $Id: cfgtransf.cc,v 1.7 2003-10-16 01:40:11 edwards Exp $
/*! \file
 *  \brief Many-to-many gauge transformation routine
 */

#include "chroma.h"

using namespace QDP;

//! Many-to-many gauge transformation routine
/*! \defgroup cfgtransf Tranformation routine
 *  \ingroup main
 *
 * Main program for transforming gauge formats
 */

int main(int argc, char **argv)
{
  // Put the machine into a known state
  QDP_initialize(&argc, &argv);

  START_CODE("cfgtransf");
  
  multi1d<int> nrow(Nd);
  QDPIO::cout << "Enter lattice size\n";
  QDPIO::cin >> nrow;
  
  // Setup QDP
  Layout::setLattSize(nrow);
  Layout::create();

  XMLFileWriter xml_out("cfgtransf.xml");
  push(xml_out, "cfgtransf");

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
	      << " (20) Single-precision UKQCD gauge configuration FE\n";
  QDPIO::cin >> input_type;
  

#if 0
  Real SchrPhiMult;            /* Multiplier for Schr. BC fields */
  bool SchrFermP;           /* Set Schroedinger fermion phases? */
  multi1d<Real> theta(Nd-1);          /* Angles (units of pi) for Schroedinger phases */

  // Cannot handle Symanzik or Schroedinger BC yet
  if ( input_type == 9 )
  {
    QDPIO::cout << "Enter level of Symanzik improvement of gauge action\n"
		<< "  ( 0) No improvement (1x1 plaquette Wilson action)\n"
		<< "  (-1) Arbitrary coeff for 1x2 rectangle\n"
		<< "  ( 1) Tree level coeff for 1x2 rectangle\n"
		<< "  (-2) Arbitrary coeff. for 1x2 rect. and parallelogram\n"
		<< "  ( 2) 1-loop coeff. for 1x2 rect. and parallelogram\n";
    QDPIO::cin >> GlueImp;
  
    QDPIO::cout << "Enter Schroedinger boundary type\n"
		<< "  (0) Do not use Schroedinger BC\n"
		<< "  (1) Zero field on boundary\n"
		<< "  (2) Canonical alpha collaboration field\n";
    QDPIO::cin >> SchrFun;
  
    QDPIO::cout << "Enter multiplier for Schroedinger boundary fields\n";
    QDPIO::cin >> SchrPhiMult;
 
    QDPIO::cout << "Enter the boundary condition direction\n";
    QDPIO::cin >> j_decay;

    QDPIO::cout << "Use fermion phases in Schroedinger boundary? [0/1]\n";
    QDPIO::cin >> SchrFermP;

    if (SchrFermP)
    {
      QDPIO::cout << "Enter the " << Nc << " thetas\n";
      QDPIO::cin >> theta;
    }
    else
    {
      theta = 0;
    }
  }
  else
  {
    GlueImp = 0;
    SchrFun = 0;
    SchrFermP = NO;
    theta = 0;
    SchrPhiMult = 1;
  }
#endif

  int output_type;
  QDPIO::cout << "Enter output Gauge field type\n"
	      << "  (1) back-end SZIN (note, on most platforms back and front-end are the same)\n"
	      << "  (2) front-end SZIN\n"
	      << "  (3) MILC config on FE\n"
	      << "  (4) QCD Archive config on FE\n"
	      << "  (5) MIT gauge config on FE\n"
	      << "  (6) Kentucky config on FE\n";
  QDPIO::cin >> output_type;
  
  string cfg_input_file;
  QDPIO::cout << "Enter input file name\n";
  QDPIO::cin >> cfg_input_file;

  string cfg_output_file;
  QDPIO::cout << "Enter output file name\n";
  QDPIO::cin >> cfg_output_file;
  
  bool RGaugeP;
  if ( input_type > 3 )
  {
    QDPIO::cout << "Random gauge transform of config?\n";
    QDPIO::cin >> RGaugeP;
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
  }

  QDPIO::cout << "I am working on it...\n";
  
  push(xml_out,"Lattis");
  Write(xml_out, Nd);
  Write(xml_out, Nc);
  Write(xml_out, nrow);
  pop(xml_out);
  
#if 0
  if ( AnisoP )
  {
    push(xml_out,"Anisotropy_parameters");
    Write(xml_out, AnisoP);
    Write(xml_out, xi_0);
    Write(xml_out, t_dir);
    pop(xml_out);
  }
#endif

#if 0
  /* Setup Schroedinger boundary fields if desired */
  if (SchrFun > 0)
  {
    SetSFbc(SchrPhiMult, SchrFermP, theta, j_decay);
    push(xml_out,"Schroed_Func_parameters");
    Write(xml_out, SchrFun);
    Write(xml_out, j_decay);
    Write(xml_out, SchrPhiMult);
    Write(xml_out, SchrFermP);
    Write(xml_out, theta);
    pop(xml_out);
  }
#endif


  //
  //  Have params, now find out the source for the gauge field 
  //
  multi1d<LatticeColorMatrix> u(Nd);
  XMLReader gauge_xml_in;
    

  switch (input_type)
  {
  case 1:
    push(xml_out,"Free_Field");
    Write(xml_out, input_type);
    pop(xml_out);
    QDPIO::cout << "Fill u with free field\n";
    u = 1;
    break;

  case 2:
    push(xml_out,"Free_Field_with_random_gauge_transformation");
    Write(xml_out, input_type);
    pop(xml_out);
    QDPIO::cout << "Fill u with random gauge transformed free field\n";
    u = 1;
    rgauge(u);
    break;

  case 3:
    push(xml_out,"Semi-Haar_measure");
    Write(xml_out, input_type);
    pop(xml_out);
    QDPIO::cout << "Fill u with semi-Haar\n";
    HotSt(u);
    break;

  case 4:
  case 5:
    push(xml_out,"SZIN_configuration");
    Write(xml_out, input_type);
    Write(xml_out, cfg_input_file);
    pop(xml_out);
    
    QDPIO::cout << "Read SZIN u from FE file " << cfg_input_file << endl;
    readSzin(gauge_xml_in, u, cfg_input_file);
    break;

// case 6: ancient Illinois staggered

  case 7:
    push(xml_out,"MILC_config");
    Write(xml_out, input_type);
    Write(xml_out, cfg_input_file);
    pop(xml_out);
    QDPIO::cout << "Read MILC u from FE file\n";
    readMILC(gauge_xml_in, u, cfg_input_file);
    break;

//  case 8: Ancient columbia config

#if 0
  case 9:
    push(xml_out,"Schroedinger_BC_config");
    Write(xml_out, input_type);
    pop(xml_out);
    u = SFBndFld;
    break;
#endif

// case 10: Ancient fortran szin_config

#if 0
    // Not sure I want this ...
  case 11:
    push(xml_out,"Ascii_MILC_config");
    Write(xml_out, input_type);
    Write(xml_out, cfg_input_file);
    pop(xml_out);
    QDPIO::cout << "Read ASCII MILC u from FE file\n";
    readascii (u, cfg_input_file);
    break;
#endif

#if 0
    // Not yet...
  case 12:
    push(xml_out,"Fund_to_adj_SZIN_config");
    Write(xml_out, input_type);
    Write(xml_out, cfg_input_file);
    pop(xml_out);
    QDPIO::cout << "Read fund to adj SZIN from file " << cfg_input_file << endl;
    readFunToAdj (u, cfg_input_file);
    break;
#endif

#if 0
    // Not yet...
  case 13:
    push(xml_out,"O3_to_U1_config");
    Write(xml_out, input_type);
    Write(xml_out, cfg_input_file);
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
    Write(xml_out, input_type);
    Write(xml_out, center);
    Write(xml_out, rho);
    Write(xml_out, su2_index);
    Write(xml_out, sign);
    pop(xml_out);
    QDPIO::cout << "Create instanton configuration\n";
    instanton (u, center, rho, su2_index, sign);

  break;
#endif

  case 15:
    push(xml_out,"QCD_Archive_config");
    Write(xml_out, input_type);
    pop(xml_out);
    QDPIO::cout << "Read QCD Archive u from FE file\n";
    readArchiv(gauge_xml_in, u, cfg_input_file);
    break;

#if 0
    // Someday...
  case 16:
    push(xml_out,"MIT_config");
    Write(xml_out, input_type);
    pop(xml_out);
    QDPIO::cout << "Read MIT gauge config from FE file\n";
    readmitu (cfg_input_file, u);
    break;
#endif

#if 0
    // Not now...
  case 17:
    push(xml_out,"MIT_double_config");
    Write(xml_out, input_type);
    pop(xml_out);
    QDPIO::cout << "Read double prec. MIT gauge config from FE file\n";
    readmitdu (cfg_input_file, u);
    savern (seed_old);
    break;
#endif

#if 0
    // Not now...
  case 18:
    push(xml_out,"Kentucky_config");
    Write(xml_out, input_type);
    pop(xml_out);
    QDPIO::cout << "Read Kentucky gauge config from FE file\n";
    readkyu (cfg_input_file, u);
    savern (seed_old);
    break;
#endif

#if 0
    // Not now...
  case 19:
    push(xml_out,"UKQCD_config");
    Write(xml_out, input_type);
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
    Write(xml_out, input_type);
    pop(xml_out);
    QDPIO::cout << "Read single-precision UKQCD gauge config from FE file\n";
    readukqcdsingle (cfg_input_file, u);
    savern (seed_old);
    break;
#endif

  default:
    QDP_error_exit("unknown input type", input_type);
  }
    

  // Dump a copy of the input gauge xml
  write(xml_out, "input_gauge_header", gauge_xml_in);
  

  // So what's the plaquette?
  Double w_plaq;
  Double s_plaq;
  Double t_plaq;
  Double link;
  multi1d<DComplex> pollp(Nd);		/* Polyakov loop */
  
  MesPlq (u, w_plaq, s_plaq, t_plaq, link);
  for(int mu = 0; mu < Nd; ++mu)
    polylp (u, pollp[mu], mu);
  
  push(xml_out,"Observables");
  Write(xml_out, w_plaq);
  Write(xml_out, s_plaq);
  Write(xml_out, t_plaq);
  Write(xml_out, link);
  Write(xml_out, pollp);
  pop(xml_out);

  xml_out.flush();
  
  if ( RGaugeP )
  {
    rgauge (u);

    MesPlq (u, w_plaq, s_plaq, t_plaq, link);
    for(int mu = 0; mu < Nd; ++mu)
      polylp (u, pollp[mu], mu);
  
    push(xml_out,"Rand_Gtransf_observables");
    Write(xml_out, w_plaq);
    Write(xml_out, s_plaq);
    Write(xml_out, t_plaq);
    Write(xml_out, link);
    Write(xml_out, pollp);
    pop(xml_out);
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
    
    MesPlq(u, w_plaq, s_plaq, t_plaq, link);
    for(int mu = 0; mu < Nd; ++mu)
      polylp(u, pollp[mu], mu);
  
    push(xml_out,"Smearing_parameters");
    Write(xml_out, sm_numb);
    Write(xml_out, sm_fact);
    Write(xml_out, j_decay);
    pop(xml_out);
    push(xml_out,"Smeared_observables");
    Write(xml_out, w_plaq);
    Write(xml_out, s_plaq);
    Write(xml_out, t_plaq);
    Write(xml_out, link);
    Write(xml_out, pollp);
    pop(xml_out);
  }

  if ( HypSmearP )
  {
    multi1d<LatticeColorMatrix> u_tmp(Nd);
    for(int i = 0; i < sm_numb; ++i)
    {
      Hyp_Smear(u, u_tmp, alpha1, alpha2, alpha3, BlkAccu, BlkMax);
      u = u_tmp;
    }
    
    MesPlq(u, w_plaq, s_plaq, t_plaq, link);
    for(int mu = 0; mu < Nd; ++mu)
      polylp(u, pollp[mu], mu);
  
    push(xml_out,"HypSmearing_parameters");
    Write(xml_out, sm_numb);
    Write(xml_out, alpha1);
    Write(xml_out, alpha2);
    Write(xml_out, alpha3);
    pop(xml_out);
    push(xml_out,"HypSmeared_observables");
    Write(xml_out, w_plaq);
    Write(xml_out, s_plaq);
    Write(xml_out, t_plaq);
    Write(xml_out, link);
    Write(xml_out, pollp);
    pop(xml_out);
  }

  int nrl_gf;
  if ( GFixP  )
  {
//    gfix(u, j_decay, GFAccu, GFMax, nrl_gf, OrlxDo, OrPara);
    QDP_error_exit("gauge fixing not implemented");

    MesPlq(u, w_plaq, s_plaq, t_plaq, link);
    for(int mu = 0; mu < Nd; ++mu)
      polylp(u, pollp[mu], mu);
  
    push(xml_out,"Gauge_fixing_parameters");
    Write(xml_out, j_decay);
    Write(xml_out, GFAccu);
    Write(xml_out, GFMax);
    Write(xml_out, OrlxDo);
    Write(xml_out, OrPara);
    pop(xml_out);
    push(xml_out,"Gauge_fixed_observables");
    Write(xml_out, w_plaq);
    Write(xml_out, s_plaq);
    Write(xml_out, t_plaq);
    Write(xml_out, link);
    Write(xml_out, pollp);
    pop(xml_out);
    
    push(xml_out,"Relaxation_iterations_in_GFIX");
    Write(xml_out, nrl_gf);
    pop(xml_out);
  }


  // Compute the plaquette again
  MesPlq (u, w_plaq, s_plaq, t_plaq, link);
  for(int mu = 0; mu < Nd; ++mu)
    polylp (u, pollp[mu], mu);
  
  xml_out.flush();
  
  /* Now write parameters to file cfg_output_file */
  switch (output_type)
  {
  case 1:
  case 2:
  {
    int TotalTrj;
    SzinGauge_t szin_out;
    szinGaugeInit(szin_out);

    QDPIO::cout << "Enter TotalTrj\n";
    QDPIO::cin >> TotalTrj;
    if (TotalTrj > 0)
    {
      QDPIO::cout << "Enter NOver\n";
      QDPIO::cin >> szin_out.NOver;
      QDPIO::cout << "Enter BetaMC\n";
      QDPIO::cin >> szin_out.BetaMC;
      QDPIO::cout << "Enter BetaMD\n";
      QDPIO::cin >> szin_out.BetaMD;
      QDPIO::cout << "Enter KappaMC\n";
      QDPIO::cin >> szin_out.KappaMC;
      QDPIO::cout << "Enter KappaMD\n";
      QDPIO::cin >> szin_out.KappaMD;
      QDPIO::cout << "Enter dt\n";
      QDPIO::cin >> szin_out.dt;
      QDPIO::cout << "Enter MesTrj\n";
      QDPIO::cin >> szin_out.MesTrj;
      QDPIO::cout << "Enter Nf\n";
      QDPIO::cin >> szin_out.Nf;
      QDPIO::cout << "Enter Npf\n";
      QDPIO::cin >> szin_out.Npf;
      QDPIO::cout << "Enter RefMomTrj\n";
      QDPIO::cin >> szin_out.RefMomTrj;
      QDPIO::cout << "Enter RefFnoiseTrj\n";
      QDPIO::cin >> szin_out.RefFnoiseTrj;
      QDPIO::cout << "Enter LamPl\n";
      QDPIO::cin >> szin_out.LamPl;
      QDPIO::cout << "Enter LamMi\n";
      QDPIO::cin >> szin_out.LamMi;
      QDPIO::cout << "Enter AlpLog\n";
      QDPIO::cin >> szin_out.AlpLog;
      QDPIO::cout << "Enter AlpExp\n";
      QDPIO::cin >> szin_out.AlpExp;
      QDPIO::cout << "Enter seed\n";
      QDPIO::cin >> szin_out.seed;
    }

    writeSzin(szin_out, u, cfg_output_file);
  }
  break;

  case 3:
  {
    /* Write a MILC format file on FE */
    MILCGauge_t milc_out;
    MILCGaugeInit(milc_out);
    writeMILC (milc_out, u, cfg_output_file);
  }
  break;

  case 4:
  {
    /* Write a QCD Archive format file on FE */
    ArchivGauge_t arc_out;
    archivGaugeInit(arc_out);
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

#if 0
    // Not yet...
  case 6:
    /* Write a Kentucky gauge format file on FE */
    wrtkyu (cfg_output_file, u);
    break;
#endif

  default:
    QDP_error_exit("unknown output type", output_type);
  }

  pop(xml_out);
        
  // Time to bolt
  QDP_finalize();

  exit(0);
}
