// $Id: szin_io.cc,v 1.1 2003-10-08 02:30:25 edwards Exp $

/*! \file
 *  \brief Reader/writers for szin headers
 */

#include "chromabase.h"
#include "io/szin_io.h"
#include "primitives.h"

#include <string>
using std::string;

using namespace QDP;

//! Read a SZIN header from XML
/*!
 * \ingroup io
 *
 *   Gauge field layout is (fortran ordering)
 *     u(real/imag,color_row,color_col,site,cb,Nd)
 *         = u(2,Nc,Nc,VOL_CB,2,4)
 *
 *
 * \param xml        xml reader holding where we will extract header ( Read )
 * \param path       path to structure
 * \param cfg_file   path ( Read )
 */    

//! Source header read
void read(XMLReader& xml, const string& path, SzinGauge_t& header);

//! Source header writer
void write(XMLWriter& xml, const string& path, const SzinGauge_t& header);

void writeSzin(XMLBufferWriter& xml, SzinHeader_t head)
{
  multi1d<int> nrow_old; /* Lattice size (from CFGIN) */
  int Nd_old; /* Number of spacetime dimensions (from CFGIN) */
  int Nc_old; /* Number of colors (from CFGIN) */

  multi2d<Real32> wstat_old(41, 20); /* On-line statistical accumulators */

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

  int cfg_record_size = 0;
  int cfg_version;

  string banner_tmp;
  string date_tmp;

  START_CODE("writeSzin");

  // Extract info from the XML header
  try 
  {
    XMLReader  xml_in_top(xml);
    XMLReader  xml_in(xml_in_top, "/szin");

    read(xml_in,"cfg_version",cfg_version);

    read(xml_in,"date",date_tmp);
    read(xml_in,"banner",banner_tmp);

    read(xml_in,"FermTypeP",FermTypeP_old);
    read(xml_in,"Nd",Nd_old);
    read(xml_in,"Nc",Nc_old);
    read(xml_in,"BetaMC",BetaMC_old);
    read(xml_in,"BetaMD",BetaMD_old);

    read(xml_in,"KappaMC",KappaMC_old);
    read(xml_in,"KappaMD",KappaMD_old);
    read(xml_in,"MassMC",MassMC_old);
    read(xml_in,"MassMD",MassMD_old);
    read(xml_in,"dt",dt_old);
    read(xml_in,"MesTrj",MesTrj_old);
    read(xml_in,"TotalCG",TotalCG_old);
    read(xml_in,"TotalTrj",TotalTrj_old);
    read(xml_in,"spec_acc",spec_acc);

    read(xml_in,"NOver",NOver_old);
    read(xml_in,"TotalTry",TotalTry_old);
    read(xml_in,"TotalFail",TotalFail_old);
    read(xml_in,"Nf",Nf_old);
    read(xml_in,"Npf",Npf_old);
    read(xml_in,"RefMomTrj",RefMomTrj_old);
    read(xml_in,"RefFnoiseTrj",RefFnoiseTrj_old);
    read(xml_in,"LamPl",LamPl_old);
    read(xml_in,"LamMi",LamMi_old);
    read(xml_in,"AlpLog",AlpLog_old);
    read(xml_in,"AlpExp",AlpExp_old);

    read(xml_in,"nrow",nrow_old);
    read(xml_in,"seed",seed_old);
  }
  catch(const string& e)
  { 
    QDP_error_exit("Error in writeszin: %s",e.c_str());
  }


  // Check that old and new parameters are compatible
  if ( Nd_old != Nd )
    QDP_error_exit("number of dimensions specified different from configuration file: Nd_old=%d",
                   Nd_old);

  if ( Nc_old != Nc )
    QDP_error_exit("number of colors specified different from configuration file: Nc_old=%d",
                   Nc_old);

  for(int j = 0; j < Nd; ++j)
    if ( nrow_old[j] != Layout::lattSize()[j] )
      QDP_error_exit("lattice size specified different from configuration file: nrow_old[%d]=%d",
                     j,nrow_old[j]);


  // Write in the configuration along with relevant information
  BinaryWriter cfg_out(cfg_file); // for now, cfg_io_location not used

  int date_size = date_tmp.length() + 1;
  int banner_size = banner_tmp.length() + 1;

  if( date_size < 1 || date_size > 99)
    QDP_error_exit("Apparently wrong configuration file, date_size=%d",date_size);

  write(cfg_out,date_size);
  write(cfg_out,banner_size);
  write(cfg_out,cfg_record_size);

  /*
   * Write out the date & banner. They are written as int's
   */
  for(int i=0; i < date_size; ++i)
  {
    int j = date_tmp.c_str()[i];
    write(cfg_out,j);
  }

  for(int i=0; i < banner_size; ++i)
  {
    int j = banner_tmp.c_str()[i];
    write(cfg_out,j);
  }

  write(cfg_out,cfg_version);

  write(cfg_out,FermTypeP_old);
  write(cfg_out,Nd_old);
  write(cfg_out,Nc_old);
  write(cfg_out,BetaMC_old);
  write(cfg_out,BetaMD_old);
  
  write(cfg_out,KappaMC_old);
  write(cfg_out,KappaMD_old);
  write(cfg_out,MassMC_old);
  write(cfg_out,MassMD_old);
  write(cfg_out,dt_old);
  write(cfg_out,MesTrj_old);
  write(cfg_out,TotalCG_old);
  write(cfg_out,TotalTrj_old);
  write(cfg_out,spec_acc);

  write(cfg_out,NOver_old);
  write(cfg_out,TotalTry_old);
  write(cfg_out,TotalFail_old);
  write(cfg_out,Nf_old);
  write(cfg_out,Npf_old);
  write(cfg_out,RefMomTrj_old);
  write(cfg_out,RefFnoiseTrj_old);
  write(cfg_out,LamPl_old);
  write(cfg_out,LamMi_old);
  write(cfg_out,AlpLog_old);
  write(cfg_out,AlpExp_old);

  write(cfg_out,nrow_old);
  write(cfg_out,seed_old);

  wstat_old = 0;
  write(cfg_out,wstat_old);


  /*
   *  SZIN stores data "checkerboarded".  We must therefore "fake" a checkerboarding
   */

  ColorMatrix u_tmp, u_old;
  
  multi1d<int> lattsize_cb = Layout::lattSize();
  lattsize_cb[0] /= 2;		// Evaluate the coords on the checkerboard lattice

  // The slowest moving index is the direction
  for(int j = 0; j < Nd; j++)
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

	u_old = peekSite(u[j], coord); // Put it into the correct place
	u_tmp = transpose(u_old); // Take the transpose

	write(cfg_out,u_tmp); 	// Write in an SU(3) matrix
      }
  }

  cfg_out.close();
  END_CODE("writeSzin");
}

