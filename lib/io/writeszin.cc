// $Id: writeszin.cc,v 1.2 2003-10-08 04:37:50 edwards Exp $

/*! \file
 *  \brief Write out a configuration written by SZIN up to configuration version 7.
 */

#include "chromabase.h"
#include "io/szin_io.h"
#include "io/writeszin.h"
#include "primitives.h"
#include "qdp_util.h"    // from QDP

#include <string>
using std::string;

using namespace QDP;

//! Write a SZIN configuration file
/*!
 * \ingroup io
 *
 *   Gauge field layout is (fortran ordering)
 *     u(real/imag,color_row,color_col,site,cb,Nd)
 *         = u(2,Nc,Nc,VOL_CB,2,4)
 *
 *
 * \param header     structure holding config info ( Modify )
 * \param u          gauge configuration ( Read )
 * \param cfg_file   path ( Read )
 */    

void writeSzin(const SzinGauge_t& header, const multi1d<LatticeColorMatrix>& u, const string& cfg_file)
{
  START_CODE("writeSzin");

  multi2d<Real32> wstat(41, 20); /* On-line statistical accumulators - write junk */
  wstat = 0;

  int cfg_record_size = 0;

  // Write in the configuration along with relevant information
  BinaryWriter cfg_out(cfg_file); // for now, cfg_io_location not used

  int date_size = header.date.length() + 1;
  int banner_size = header.banner.length() + 1;

  if( date_size < 1 || date_size > 99)
    QDP_error_exit("Apparently wrong configuration file, date_size=%d",date_size);

  write(cfg_out, date_size);
  write(cfg_out, banner_size);
  write(cfg_out, cfg_record_size);

  /*
   * Write out the date & banner. They are written as int's
   */
  for(int i=0; i < date_size; ++i)
  {
    int j = header.date.c_str()[i];
    write(cfg_out, j);
  }

  for(int i=0; i < banner_size; ++i)
  {
    int j = header.banner.c_str()[i];
    write(cfg_out, j);
  }

  write(cfg_out, header.cfg_version);

  write(cfg_out, header.FermTypeP);
  write(cfg_out, header.Nd);
  write(cfg_out, header.Nc);
  write(cfg_out, header.BetaMC);
  write(cfg_out, header.BetaMD);
  
  write(cfg_out, header.KappaMC);
  write(cfg_out, header.KappaMD);
  write(cfg_out, header.MassMC);
  write(cfg_out, header.MassMD);
  write(cfg_out, header.dt);
  write(cfg_out, header.MesTrj);
  write(cfg_out, header.TotalCG);
  write(cfg_out, header.TotalTrj);
  write(cfg_out, header.spec_acc);

  write(cfg_out, header.NOver);
  write(cfg_out, header.TotalTry);
  write(cfg_out, header.TotalFail);
  write(cfg_out, header.Nf);
  write(cfg_out, header.Npf);
  write(cfg_out, header.RefMomTrj);
  write(cfg_out, header.RefFnoiseTrj);
  write(cfg_out, header.LamPl);
  write(cfg_out, header.LamMi);
  write(cfg_out, header.AlpLog);
  write(cfg_out, header.AlpExp);

  write(cfg_out, header.nrow);
  write(cfg_out, header.seed);

  write(cfg_out, wstat);


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

	write(cfg_out, u_tmp); 	// Write in an SU(3) matrix
      }
  }

  cfg_out.close();
  END_CODE("writeSzin");
}


//! Write a SZIN configuration file
/*!
 * \ingroup io
 *
 *   Gauge field layout is (fortran ordering)
 *     u(real/imag,color_row,color_col,site,cb,Nd)
 *         = u(2,Nc,Nc,VOL_CB,2,4)
 *
 *
 * \param xml        xml writer holding config info ( Read )
 * \param u          gauge configuration ( Read )
 * \param cfg_file   path ( Read )
 */    

void writeSzin(XMLBufferWriter& xml, multi1d<LatticeColorMatrix>& u, const string& cfg_file)
{
  START_CODE("writeSzin");

  SzinGauge_t header;
  XMLReader  xml_in(xml);   // use the buffer writer to instantiate a reader
  read(xml_in, "/szin", header);

  writeSzin(header, u, cfg_file);

  END_CODE("writeSzin");
}

