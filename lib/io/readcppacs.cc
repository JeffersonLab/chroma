// $Id: readcppacs.cc,v 3.1 2007-06-10 14:40:23 edwards Exp $

/*! \file
 *  \brief Read a CPPACS gauge configuration 
 */

#include "chromabase.h"
#include "io/cppacs_io.h"
#include "io/readcppacs.h"
#include "qdp_util.h"    // from QDP

namespace Chroma {

//! Read a CPPACCPPACS configuration file
/*!
 * \ingroup io
 * 
 * Based on the C-routine ToSTDConf.c available at
 * http://www.lqa.rccp.tsukuba.ac.jp/tools.html
 *
 * Author: Andreas Juettner
 * 
 * \param header     structure holding config info ( Modify )
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readCPPACS(CPPACSGauge_t& header, multi1d<LatticeColorMatrix>& u, const string& cfg_file)
{
  START_CODE();
#define NBC 144  /* = 8*2*3*3 = number of bytes to copy */
#define NBK   0  /* = 8*2*3*0 = number of bytes to skip */
#define NBS 576  /* = 8*2*3*3*4 = number of bytes for each site */
typedef unsigned short int  n_uint16_t;
typedef unsigned int        n_uint32_t;
typedef unsigned long int   n_uint64_t;
     char *file ;
     int NX,NY,NZ,NT ;
     int ix0,ix1,iy0,iy1,iz0,iz1,it0,it1 ;
     double *config ;
  int one=1, four=4, eight=8 ;
  int data_endian = 1, machine_endian ;
  int fd ; 
  char header_ildg[1020] ;
  int c,ix,iy,iz,it,mu,n,ndata ;
  int Nsp, Ntm, Ntx, Nty, Ntz, Ntt ;
  int xseek0,xseek1,yseek0,yseek1,zseek0,zseek1,tseek0,tseek1;
  int pos;
  long lSize ;
  char *conf ;
  int site ;
  long size;

#ifdef PLAQ
  double wss,wst ;
#endif
  NT = Layout::lattSize()[3];
  NZ = Layout::lattSize()[2];
  NY = Layout::lattSize()[1];
  NX = Layout::lattSize()[0];
  BinaryFileReader cfg_in(cfg_file)  ;
  int magic_number ;
  read(cfg_in, magic_number);

  bool byterev = false;
  if( magic_number != 19920410 ) {
    QDPIO::cout << "readCPPACS: magic number invalid" << endl;
    QDPIO::cout << "Trying byte reversal" << endl;
    QDPUtil::byte_swap((void *)&magic_number, sizeof(int), 1 );
    byterev = true;
  }
    QDPIO::cout << "Magic number: " << magic_number << endl;


  if( magic_number != 19920410){
    QDP_error_exit("readCPPACS: unexpected magic number");
  }

  if(byterev)
    QDPIO::cout<<"Doing bytereversal on the links...\n" ;

  // read in the header 
  cfg_in.readArray(header_ildg,1,1020);
  /* read header and check file name */

  /* check lattice size */
  header_ildg[19] = '\0' ;  Ntm = atoi(&header_ildg[17]) ;
  header_ildg[16] = '\0' ;  Nsp = atoi(&header_ildg[14]) ;  
  if( Ntm != NT ) {
    fprintf(stderr,"readCPPACS: NT mismatch  %d %d\n",NT,Ntm) ;
    exit(1) ; }
  if( Nsp != NZ ) {
    fprintf(stderr,"readCPPACS: NX mismatch  %d %d\n",NX,Nsp) ;
    exit(1) ; }
  if( Nsp != NY ) {
    fprintf(stderr,"readCPPACS: NY mismatch  %d %d\n",NY,Nsp) ;
    exit(1) ; }
  if( Nsp != NX ) {
    fprintf(stderr,"readCPPACS: NZ mismatch  %d %d\n",NZ,Nsp) ;
    exit(1) ; }

// Read in SU(3) matrices.
// the running of the lexicographic index and the su3_matrix layout in CPPACS 
// gauge configs is  the same as in QDP 
// (the fastest running direction is the 0th direction, corresponding to the 
// x-direction)

  u = zero ; 

  LatticeColorMatrixD  uu ; 

  ColorMatrixD  uuuD ; 
  ColorMatrix  uuu ; 

  for(int site=0; site < Layout::vol(); ++site)
  {
    multi1d<int> coord = crtesn(site, Layout::lattSize()); // The coordinate
    // read in a single site
    for(int mu=0; mu < Nd; ++mu) 
      {  
	read(cfg_in, uuuD );  
	if(byterev){
	  QDPUtil::byte_swap((void *)&uuuD.elem(),sizeof(double),2*Nc*Nc);
	}
	uuu = uuuD ; 
	pokeSite(u[mu],uuu,coord); 
      }    

  }

  cfg_in.close();

  END_CODE();
}



//! Read a CPPACS configuration file
/*!
 * \ingroup io
 *
 * \param xml        xml reader holding config info ( Modify )
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readCPPACS(XMLReader& xml, multi1d<LatticeColorMatrix>& u, const string& cfg_file)
{
  START_CODE();

  CPPACSGauge_t header;

  // Read the config and its binary header
  readCPPACS(header, u, cfg_file);

#if 0
  // Now, set up the XML header. Do this by first making a buffer
  // writer that is then used to make the reader
  XMLBufferWriter  xml_buf;
//  write(xml_buf, "CPPACS", header);

  try 
 {
    xml.open(xml_buf);
  }
  catch(const string& e)
  { 
    QDP_error_exit("Error in readCPPACS: %s",e.c_str());
  }
#endif

  END_CODE();
}

}  // end namespace Chroma

