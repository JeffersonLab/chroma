// $Id: readszin.cc,v 1.4 2003-04-02 22:28:21 edwards Exp $

/*! \file
 *  \brief Read in a configuration written by SZIN up to configuration version 7.
 */

#include "chromabase.h"
#include "io/readszin.h"
#include "primitives.h"
#include "qdp_util.h"    // from QDP

using namespace QDP;

//! Read a SZIN configuration file
/*!
 * \ingroup io
 *
 *   Gauge field layout is (fortran ordering)
 *     u(real/imag,color_row,color_col,site,cb,Nd)
 *         = u(2,Nc,Nc,VOL_CB,2,4)
 *
 *
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 * \param seed_old   seed in configuration ( Modify )            
 */    

void readSzin2(multi1d<LatticeColorMatrix>& u, char cfg_file[], Seed& seed_old)
{
  ColorMatrixF u_old;
  multi1d<int> nrow_old(Nd); /* Lattice size (from CFGIN) */
  int Nd_old; /* Number of spacetime dimensions (from CFGIN) */
  int Nc_old; /* Number of colours (from CFGIN) */

  multi2d<Real> wstat(41, 20); /* On-line statistical accumulators */
  multi2d<Real32> wstat_old(41, 20); /* WStat values from CFGIN file */

  int TotalTrj_old; /* Total number of trajectories */
  int TotalCG_old; /* Total number of CG iterations */
  int FermTypeP_old; /* Fermion type (from CFGIN) */
  Real32 MesTrj_old;
  int spec_acc;
  int MesItr_old;
//int TotalItr_old;
  int NOver_old;
  int TotalTry_old;
  int TotalFail_old;
  int Npf_old;
  int RefMomTrj_old;
  int RefFnoiseTrj_old;
  Real32 MesTrj; /* Trajectories per measurement (as a Float) */
  Real32 BetaMC_old;
  Real32 BetaMD_old;
  Real32 bh_old;
  Real32 dt_old;
  Real32 KappaMC_old;
  Real32 KappaMD_old;
  Real32 MassMC_old;
  Real32 MassMD_old;
  Real32 Nf_old;
  Real32 LamPl_old;
  Real32 LamMi_old;
  Real32 AlpLog_old;
  Real32 AlpExp_old;

  int i;
  int j;
  string date;
  string banner;

  int cfg_record_size; /* not used */
  int cfg_version;
  int date_size;
  int banner_size;

  START_CODE("readSzin");

  // Read in the configuration along with relevant information
  BinaryReader cfg_in(cfg_file); // for now, cfg_io_location not used

  read(cfg_in,date_size);
  read(cfg_in,banner_size);
  read(cfg_in,cfg_record_size);

  if( date_size < 1 || date_size > 99)
    QDP_error_exit("Apparently wrong configuration file, date_size=%d",date_size);

  /*
   * This marvelous piece of code is logically equivalent to
   *  read (cfg_in) date(1:date_size), banner(1:banner_size), 
   *      cfg_version;
   */
  for(i=0; i < date_size; ++i)
  {
    read(cfg_in,j);
    date[i] = j;
  }
  date[date_size] = '\0';

  for(i=0; i < banner_size; ++i)
  {
    read(cfg_in,j);
    banner[i] = j;
  }
  banner[banner_size] = '\0';

  read(cfg_in,cfg_version);

  switch(cfg_version) /* just add new cases if the CFG format changes */
  {
  case 1:
    read(cfg_in,Nd_old); 
    read(cfg_in,Nc_old); 
    read(cfg_in,BetaMC_old); 
    read(cfg_in,bh_old); 
    read(cfg_in,dt_old); 
    read(cfg_in,MesTrj);
    read(cfg_in,KappaMC_old);
    TotalTrj_old = 0;
    BetaMD_old = BetaMC_old;
    KappaMD_old = KappaMC_old;
    MassMC_old = 0;
    MassMD_old = 0;
    spec_acc = 1;
    FermTypeP_old = WILSON_FERMIONS;
    NOver_old = 0;
    TotalTry_old = 0;
    TotalFail_old = 0;
    Nf_old = 0;
    Npf_old = 0;
    RefMomTrj_old = 0;
    RefFnoiseTrj_old = 0;
    LamPl_old = 0;
    LamMi_old = 0;
    AlpLog_old = 0;
    AlpExp_old = 0;
    break;
  case 2:
    read(cfg_in,Nd_old); 
    read(cfg_in,Nc_old); 
    read(cfg_in,BetaMC_old); 
    read(cfg_in,bh_old); 
    read(cfg_in,dt_old); 
    read(cfg_in,MesTrj);
    read(cfg_in,KappaMC_old);
    read(cfg_in,TotalCG_old); 
    read(cfg_in,TotalTrj_old); 
    BetaMD_old = BetaMC_old;
    KappaMD_old = KappaMC_old;
    MassMC_old = 0;
    MassMD_old = 0;
    FermTypeP_old = WILSON_FERMIONS;
    NOver_old = 0;
    TotalTry_old = 0;
    TotalFail_old = 0;
    Nf_old = 0;
    Npf_old = 0;
    RefMomTrj_old = 0;
    RefFnoiseTrj_old = 0;
    LamPl_old = 0;
    LamMi_old = 0;
    AlpLog_old = 0;
    AlpExp_old = 0;
    break;
  case 3:
    read(cfg_in,Nd_old); 
    read(cfg_in,Nc_old); 
    read(cfg_in,BetaMC_old); 
    read(cfg_in,bh_old); 
    read(cfg_in,dt_old); 
    read(cfg_in,MesTrj);
    read(cfg_in,KappaMC_old);
    read(cfg_in,TotalCG_old); 
    read(cfg_in,TotalTrj_old); 
    read(cfg_in,spec_acc);
    BetaMD_old = BetaMC_old;
    KappaMD_old = KappaMC_old;
    MassMC_old = 0;
    MassMD_old = 0;
    FermTypeP_old = WILSON_FERMIONS;
    NOver_old = 0;
    TotalTry_old = 0;
    TotalFail_old = 0;
    Nf_old = 0;
    Npf_old = 0;
    RefMomTrj_old = 0;
    RefFnoiseTrj_old = 0;
    LamPl_old = 0;
    LamMi_old = 0;
    AlpLog_old = 0;
    AlpExp_old = 0;
    break;
  case 4:
    read(cfg_in,Nd_old); 
    read(cfg_in,Nc_old); 
    read(cfg_in,BetaMC_old); 
    read(cfg_in,BetaMD_old); 
    read(cfg_in,bh_old); 
    read(cfg_in,dt_old); 
    read(cfg_in,MesTrj);
    read(cfg_in,KappaMC_old); 
    read(cfg_in,KappaMD_old); 
    read(cfg_in,TotalCG_old); 
    read(cfg_in,TotalTrj_old); 
    read(cfg_in,spec_acc);
    MassMC_old = 0;
    MassMD_old = 0;
    FermTypeP_old = WILSON_FERMIONS;
    NOver_old = 0;
    TotalTry_old = 0;
    TotalFail_old = 0;
    Nf_old = 0;
    Npf_old = 0;
    RefMomTrj_old = 0;
    RefFnoiseTrj_old = 0;
    LamPl_old = 0;
    LamMi_old = 0;
    AlpLog_old = 0;
    AlpExp_old = 0;
    break;
  case 5:
    read(cfg_in,FermTypeP_old); 
    read(cfg_in,Nd_old); 
    read(cfg_in,Nc_old);
    read(cfg_in,BetaMC_old); 
    read(cfg_in,BetaMD_old);
    read(cfg_in,KappaMC_old); 
    read(cfg_in,KappaMD_old);
    read(cfg_in,MassMC_old); 
    read(cfg_in,MassMD_old);
    read(cfg_in,dt_old); 
    read(cfg_in,MesTrj_old); 
    read(cfg_in,TotalCG_old); 
    read(cfg_in,TotalTrj_old);
    NOver_old = 0;
    TotalTry_old = 0;
    TotalFail_old = 0;
    Nf_old = 0;
    Npf_old = 0;
    RefMomTrj_old = 0;
    RefFnoiseTrj_old = 0;
    LamPl_old = 0;
    LamMi_old = 0;
    AlpLog_old = 0;
    AlpExp_old = 0;
    break;
  case 6:
    read(cfg_in,FermTypeP_old); 
    read(cfg_in,Nd_old); 
    read(cfg_in,Nc_old); 
    read(cfg_in,BetaMC_old); 
    read(cfg_in,BetaMD_old);
    read(cfg_in,KappaMC_old); 
    read(cfg_in,KappaMD_old);
    read(cfg_in,MassMC_old); 
    read(cfg_in,MassMD_old);
    read(cfg_in,dt_old); 
    read(cfg_in,MesItr_old); 
    read(cfg_in,TotalCG_old); 
    read(cfg_in,TotalTrj_old); 
    read(cfg_in,spec_acc);
    read(cfg_in,NOver_old); 
    read(cfg_in,TotalTry_old); 
    read(cfg_in,TotalFail_old);
    Nf_old = 0;
    Npf_old = 0;
    RefMomTrj_old = 0;
    RefFnoiseTrj_old = 0;
    LamPl_old = 0;
    LamMi_old = 0;
    AlpLog_old = 0;
    AlpExp_old = 0;
    break;

  case 7:
    read(cfg_in,FermTypeP_old);
    read(cfg_in,Nd_old);
    read(cfg_in,Nc_old);
    read(cfg_in,BetaMC_old);
    read(cfg_in,BetaMD_old);

    read(cfg_in,KappaMC_old);
    read(cfg_in,KappaMD_old);
    read(cfg_in,MassMC_old);
    read(cfg_in,MassMD_old);
    read(cfg_in,dt_old);
    read(cfg_in,MesTrj_old);
    read(cfg_in,TotalCG_old);
    read(cfg_in,TotalTrj_old);
    read(cfg_in,spec_acc);

    read(cfg_in,NOver_old);
    read(cfg_in,TotalTry_old);
    read(cfg_in,TotalFail_old);
    read(cfg_in,Nf_old);
    read(cfg_in,Npf_old);
    read(cfg_in,RefMomTrj_old);
    read(cfg_in,RefFnoiseTrj_old);
    read(cfg_in,LamPl_old);
    read(cfg_in,LamMi_old);
    read(cfg_in,AlpLog_old);
    read(cfg_in,AlpExp_old);
    break;
  default:
    QDP_error_exit("configuration file version is invalid: version=%d",cfg_version);
  }

  // Check that old and new parameters are compatible
  if ( Nd_old != Nd )
    QDP_error_exit("number of dimensions specified different from configuration file: Nd_old=%d",
                   Nd_old);

  if ( Nc_old != Nc )
    QDP_error_exit("number of colors specified different from configuration file: Nc_old=%d",
                   Nc_old);

  read(cfg_in,nrow_old);

  for(j = 0; j < Nd; ++j)
    if ( nrow_old[j] != Layout::lattSize()[j] )
      QDP_error_exit("lattice size specified different from configuration file: nrow_old[%d]=%d",
                     j,nrow_old[j]);


  read(cfg_in,seed_old);
  read(cfg_in,wstat_old);


  /*
   *  Szin stores data "checkerboarded".  We must therefore "undo" the checkerboarding
   *  We use as a model the propagator routines
   */

  ColorMatrix u_tmp;
  
  multi1d<int> lattsize_cb = Layout::lattSize();
  lattsize_cb[0] /= 2;		// Evaluate the coords on the checkerboard lattice


  // The slowest moving index is the direction
  for(j = 0; j < Nd; j++)
  {
    for(int cb=0; cb < 2; ++cb)
      for(int sitecb=0; sitecb < Layout::vol()/2; ++sitecb)
      {
	multi1d<int> coord = crtesn(sitecb, lattsize_cb); // The coordinate
      
	// construct the checkerboard offset
	int sum = 0;
	for(int m=1; m<Nd; m++)
	  sum += coord[m];

	// The true lattice x-coord
	coord[0] = 2*coord[0] + ((sum + cb) & 1);

	read(cfg_in,u_old); 	// Read in an SU(3) matrix

	u_tmp = transpose(u_old); // Take the transpost
	pokeSite(u[j], u_tmp, coord); // Put it into the correct place
      }
  }

  cfg_in.close();

  END_CODE("readSzin");
}


QDP_END_NAMESPACE();
