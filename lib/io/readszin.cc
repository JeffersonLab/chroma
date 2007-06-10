// $Id: readszin.cc,v 3.1 2007-06-10 14:40:23 edwards Exp $

/*! \file
 *  \brief Read in a configuration written by SZIN up to configuration version 7.
 */

#include "chromabase.h"
#include "io/szin_io.h"
#include "io/readszin.h"
// #include "io/param_io.h"
#include "qdp_util.h"    // from QDP

namespace Chroma {

#define SZIN_WILSON_FERMIONS  1


//! Read a SZIN configuration file
/*!
 * \ingroup io
 *
 *   Gauge field layout is (fortran ordering)
 *     u(real/imag,color_row,color_col,site,cb,Nd)
 *         = u(2,Nc,Nc,VOL_CB,2,4)
 *
 *
 * \param header     structure holding config info ( Modify )
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readSzin(SzinGauge_t& header, multi1d<LatticeColorMatrix>& u, const string& cfg_file)
{
  START_CODE();

  int cfg_record_size; // must read but will ignore - not used
  int date_size;
  int banner_size;
  Real32 bh = 0;       // old beta used for higgs term - not used

  // Read in the configuration along with relevant information
  BinaryFileReader cfg_in(cfg_file); // for now, cfg_io_location not used

  read(cfg_in,date_size);
  read(cfg_in,banner_size);
  read(cfg_in,cfg_record_size);

  if( date_size < 1 || date_size > 99)
  {
    QDPIO::cerr << __func__ 
		<< ": apparently wrong SZIN configuration file, date_size=" 
		<< date_size << endl;
    QDP_abort(1);
  }

  /*
   * Read in the date & banner. They are written as int's. Use a new
   * instead of declaring a size in the constructor for the strings 
   * to avoid extra space. I found that nulls in the strings made xmlreader
   * blow up down below.
   */
  char *date_tmp = new char[date_size+1];
  for(int i=0; i < date_size; ++i)
  {
    int j;
    read(cfg_in,j);
    date_tmp[i] = j;
  }
  date_tmp[date_size] = '\0';
  header.date = date_tmp;
  delete[] date_tmp;

  char *banner_tmp = new char[banner_size+1];
  for(int i=0; i < banner_size; ++i)
  {
    int j;
    read(cfg_in,j);
    banner_tmp[i] = j;
  }
  banner_tmp[banner_size] = '\0';
  header.banner = banner_tmp;
  delete[] banner_tmp;

  read(cfg_in, header.cfg_version);

  switch(header.cfg_version) /* just add new cases if the CFG format changes */
  {
  case 1:
    read(cfg_in, header.Nd); 
    read(cfg_in, header.Nc); 
    read(cfg_in, header.BetaMC); 
    read(cfg_in, bh); 
    read(cfg_in, header.dt); 
    read(cfg_in, header.MesTrj);
    read(cfg_in, header.KappaMC);
    header.BetaMD = header.BetaMC;
    header.KappaMD = header.KappaMC;
    header.FermTypeP = SZIN_WILSON_FERMIONS;
    break;

  case 2:
    read(cfg_in, header.Nd); 
    read(cfg_in, header.Nc); 
    read(cfg_in, header.BetaMC); 
    read(cfg_in, bh); 
    read(cfg_in, header.dt); 
    read(cfg_in, header.MesTrj);
    read(cfg_in, header.KappaMC);
    read(cfg_in, header.TotalCG); 
    read(cfg_in, header.TotalTrj); 
    header.BetaMD = header.BetaMC;
    header.KappaMD = header.KappaMC;
    header.FermTypeP = SZIN_WILSON_FERMIONS;
    break;

  case 3:
    read(cfg_in, header.Nd); 
    read(cfg_in, header.Nc); 
    read(cfg_in, header.BetaMC); 
    read(cfg_in, bh); 
    read(cfg_in, header.dt); 
    read(cfg_in, header.MesTrj);
    read(cfg_in, header.KappaMC);
    read(cfg_in, header.TotalCG); 
    read(cfg_in, header.TotalTrj); 
    read(cfg_in, header.spec_acc);
    header.BetaMD = header.BetaMC;
    header.KappaMD = header.KappaMC;
    header.FermTypeP = SZIN_WILSON_FERMIONS;
    break;

  case 4:
    read(cfg_in, header.Nd); 
    read(cfg_in, header.Nc); 
    read(cfg_in, header.BetaMC); 
    read(cfg_in, header.BetaMD); 
    read(cfg_in, bh); 
    read(cfg_in, header.dt); 
    read(cfg_in, header.MesTrj);
    read(cfg_in, header.KappaMC); 
    read(cfg_in, header.KappaMD); 
    read(cfg_in, header.TotalCG); 
    read(cfg_in, header.TotalTrj); 
    read(cfg_in, header.spec_acc);
    header.FermTypeP = SZIN_WILSON_FERMIONS;
    break;

  case 5:
    read(cfg_in, header.FermTypeP); 
    read(cfg_in, header.Nd); 
    read(cfg_in, header.Nc);
    read(cfg_in, header.BetaMC); 
    read(cfg_in, header.BetaMD);
    read(cfg_in, header.KappaMC); 
    read(cfg_in, header.KappaMD);
    read(cfg_in, header.MassMC); 
    read(cfg_in, header.MassMD);
    read(cfg_in, header.dt); 
    read(cfg_in, header.MesTrj); 
    read(cfg_in, header.TotalCG); 
    read(cfg_in, header.TotalTrj);
    break;

  case 6:
    read(cfg_in, header.FermTypeP); 
    read(cfg_in, header.Nd); 
    read(cfg_in, header.Nc); 
    read(cfg_in, header.BetaMC); 
    read(cfg_in, header.BetaMD);
    read(cfg_in, header.KappaMC); 
    read(cfg_in, header.KappaMD);
    read(cfg_in, header.MassMC); 
    read(cfg_in, header.MassMD);
    read(cfg_in, header.dt); 
    read(cfg_in, header.MesTrj); 
    read(cfg_in, header.TotalCG); 
    read(cfg_in, header.TotalTrj); 
    read(cfg_in, header.spec_acc);
    read(cfg_in, header.NOver); 
    read(cfg_in, header.TotalTry); 
    read(cfg_in, header.TotalFail);
    break;

  case 7:
    read(cfg_in, header.FermTypeP);
    read(cfg_in, header.Nd);
    read(cfg_in, header.Nc);
    read(cfg_in, header.BetaMC);
    read(cfg_in, header.BetaMD);

    read(cfg_in, header.KappaMC);
    read(cfg_in, header.KappaMD);
    read(cfg_in, header.MassMC);
    read(cfg_in, header.MassMD);
    read(cfg_in, header.dt);
    read(cfg_in, header.MesTrj);
    read(cfg_in, header.TotalCG);
    read(cfg_in, header.TotalTrj);
    read(cfg_in, header.spec_acc);

    read(cfg_in, header.NOver);
    read(cfg_in, header.TotalTry);
    read(cfg_in, header.TotalFail);
    read(cfg_in, header.Nf);
    read(cfg_in, header.Npf);
    read(cfg_in, header.RefMomTrj);
    read(cfg_in, header.RefFnoiseTrj);
    read(cfg_in, header.LamPl);
    read(cfg_in, header.LamMi);
    read(cfg_in, header.AlpLog);
    read(cfg_in, header.AlpExp);
    break;

  default:
    QDPIO::cerr << __func__ << ": SZIN configuration file version is invalid: version"
		<< header.cfg_version << endl;
    QDP_abort(1);
  }

  // Check that old and new parameters are compatible
  if ( Nd != header.Nd )
  {
    QDPIO::cerr << __func__ 
		<< ": num dimensions different from SZIN config file: header.Nd=" 
		<< header.Nd << endl;
    QDP_abort(1);
  }

  if ( Nc != header.Nc )
  {
    QDPIO::cerr << __func__ 
		<< ": number of colors specified different from SZIN config file: header.Nc="
		<< header.Nc << endl;
    QDP_abort(1);
  }

  header.nrow.resize(Nd);
  read(cfg_in, header.nrow, Nd);

  for(int j = 0; j < Nd; ++j)
    if ( header.nrow[j] != Layout::lattSize()[j] )
    {
      QDPIO::cerr << __func__ 
		  << ": lattice size specified different from configuration file: nrow[" 
		  << j << "]=" << header.nrow[j] << endl;
      QDP_abort(1);
    }

  read(cfg_in, header.seed);

  multi1d<Real32> wstat(41*20); /* On-line statistical accumulators - throw away */
  read(cfg_in, wstat, wstat.size());    // will not use


  /*
   *  Szin stores data "checkerboarded".  We must therefore "undo" the checkerboarding
   *  We use as a model the propagator routines
   */
  u.resize(Nd);

  multi1d<int> lattsize_cb = Layout::lattSize();
  lattsize_cb[0] /= 2;		// Evaluate the coords on the checkerboard lattice

  // The slowest moving index is the direction
  for(int j = 0; j < Nd; j++)
  {
    LatticeColorMatrixF u_old;
  
    for(int cb=0; cb < 2; ++cb) { 
      for(int sitecb=0; sitecb < Layout::vol()/2; ++sitecb)
      {
	multi1d<int> coord = crtesn(sitecb, lattsize_cb); // The coordinate
      
	// construct the checkerboard offset
	int sum = 0;
	for(int m=1; m<Nd; m++)
	  sum += coord[m];

	// The true lattice x-coord
	coord[0] = 2*coord[0] + ((sum + cb) & 1);

	read(cfg_in, u_old, coord); 	// Read in an SU(3) matrix into coord
      }
    }
    LatticeColorMatrix u_old_prec(u_old);
   
    u[j] = transpose(u_old_prec);            // Take the transpose
  }

  cfg_in.close();

  END_CODE();
}



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
  START_CODE();

  SzinGauge_t header;

  // Read the config and its binary header
  readSzin(header, u, cfg_file);

  // Now, set up the XML header. Do this by first making a buffer
  // writer that is then used to make the reader
  XMLBufferWriter  xml_buf;
  write(xml_buf, "szin", header);

  try 
  {
    xml.open(xml_buf);
  }
  catch(const string& e)
  { 
    QDPIO::cerr << __func__ << ": Error in readszin: " << e.c_str() << endl;
    QDP_abort(1);
  }

  END_CODE();
}

}  // end namespace Chroma
