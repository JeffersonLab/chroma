// $Id: readszin.cc,v 1.13 2003-08-27 22:08:41 edwards Exp $

/*! \file
 *  \brief Read in a configuration written by SZIN up to configuration version 7.
 */

#include "chromabase.h"
#include "io/readszin.h"
#include "primitives.h"
#include "qdp_util.h"    // from QDP

#include <string>
using std::string;

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
 * \param xml        xml reader holding config info ( Modify )
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readSzin(XMLReader& xml, multi1d<LatticeColorMatrix>& u, const string& cfg_file)
{
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
  Seed   seed_old;

  int i;
  int j;

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
   * Read in the date & banner. They are written as int's. Use a new
   * instead of declaring a size in the constructor for the strings 
   * to avoid extra space. I found that nulls in the strings made xmlreader
   * blow up down below.
   */
  char *date_tmp = new char[date_size+1];
  for(i=0; i < date_size; ++i)
  {
    read(cfg_in,j);
    date_tmp[i] = j;
  }
  date_tmp[date_size] = '\0';
  string date = date_tmp;
  delete[] date_tmp;

  char *banner_tmp = new char[banner_size+1];
  for(i=0; i < banner_size; ++i)
  {
    read(cfg_in,j);
    banner_tmp[i] = j;
  }
  banner_tmp[banner_size] = '\0';
  string banner = banner_tmp;
  delete[] banner_tmp;

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

  ColorMatrix u_tmp, u_old;
  
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


  // Now, set up the XML header. Do this by first making a buffer
  // writer that is then used to make the reader
  XMLBufferWriter  xml_buf;
  push(xml_buf, "szin");

  write(xml_buf,"cfg_version",cfg_version);

  write(xml_buf,"date",date);
  write(xml_buf,"banner",banner);

  write(xml_buf,"FermTypeP",FermTypeP_old);
  write(xml_buf,"Nd",Nd_old);
  write(xml_buf,"Nc",Nc_old);
  write(xml_buf,"BetaMC",BetaMC_old);
  write(xml_buf,"BetaMD",BetaMD_old);

  write(xml_buf,"KappaMC",KappaMC_old);
  write(xml_buf,"KappaMD",KappaMD_old);
  write(xml_buf,"MassMC",MassMC_old);
  write(xml_buf,"MassMD",MassMD_old);
  write(xml_buf,"dt",dt_old);
  write(xml_buf,"MesTrj",MesTrj_old);
  write(xml_buf,"TotalCG",TotalCG_old);
  write(xml_buf,"TotalTrj",TotalTrj_old);
  write(xml_buf,"spec_acc",spec_acc);

  write(xml_buf,"NOver",NOver_old);
  write(xml_buf,"TotalTry",TotalTry_old);
  write(xml_buf,"TotalFail",TotalFail_old);
  write(xml_buf,"Nf",Nf_old);
  write(xml_buf,"Npf",Npf_old);
  write(xml_buf,"RefMomTrj",RefMomTrj_old);
  write(xml_buf,"RefFnoiseTrj",RefFnoiseTrj_old);
  write(xml_buf,"LamPl",LamPl_old);
  write(xml_buf,"LamMi",LamMi_old);
  write(xml_buf,"AlpLog",AlpLog_old);
  write(xml_buf,"AlpExp",toFloat(AlpExp_old));

  write(xml_buf,"nrow",nrow_old);
  write(xml_buf,"seed",seed_old);

  pop(xml_buf);

  try 
  {
    xml.open(xml_buf);
  }
  catch(const string& e)
  { 
    QDP_error_exit("Error in readszin: %s",e.c_str());
  }

  END_CODE("writeSzin");
}

