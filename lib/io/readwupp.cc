/*! \file
 *  \brief Read in a configuration written by BMW up to configuration 
           version 7 (not the compressed format).
 */

#include "chromabase.h"

#include "io/readwupp.h"
// #include "io/param_io.h"
#include "qdp_util.h"    // from QDP
#include "util/gauge/unit_check.h"

namespace Chroma {


//! Read a WUPP configuration file
/*!
 * \ingroup io
 *  Code written by Sara collins
 *  Need to first run  a conversion code supplied by Kalman Szabo
 *
 * \param header     structure holding config info ( Modify )
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    

  void readWupp(multi1d<LatticeColorMatrix>& u, const string& cfg_file) 
{
  START_CODE();

  int NT = Layout::lattSize()[3];
  int NZ = Layout::lattSize()[2];
  int NY = Layout::lattSize()[1];
  int NX = Layout::lattSize()[0];

  BinaryFileReader cfg_in(cfg_file)  ;

  int nxw,nyw,nzw,ntw;

  read(cfg_in, nxw);
  read(cfg_in, nyw);
  read(cfg_in, nzw);
  read(cfg_in, ntw);

  bool byterev = false;
  if( nxw != NX ) {
    QDPIO::cout << "nxw " << nxw <<endl;
    QDPIO::cout << "readWupp: lattice dimension invalid" << endl;
    QDPIO::cout << "Trying byte reversal" << endl;
    QDPUtil::byte_swap((void *)&nxw, sizeof(int), 1 );
    byterev = true;
  }
    QDPIO::cout << "NX: " << nxw << endl;


  if( nxw != NX){
    QDP_error_exit("readWupp: unexpected lattice dimension");
  }

  if(byterev)
    {
      QDPIO::cout<<"Doing bytereversal on the links...\n" ;
      QDPUtil::byte_swap((void *)&nyw, sizeof(int), 1 );
      QDPUtil::byte_swap((void *)&nzw, sizeof(int), 1 );
      QDPUtil::byte_swap((void *)&ntw, sizeof(int), 1 );
      QDPIO::cout << " Ny Nz Nt " << nyw << " " << nzw << " " <<ntw <<endl;
    }

  if( ntw != NT ) {
    fprintf(stderr,"readWupp: NT mismatch  %d %d\n",NT,ntw) ;
    exit(1) ; }
  if( nyw != NY ) {
    fprintf(stderr,"readWupp: NY mismatch  %d %d\n",NY,nyw) ;
    exit(1) ; }
  if( nzw != NX ) {
    fprintf(stderr,"readWupp: NZ mismatch  %d %d\n",NZ,nzw) ;
    exit(1) ; }

// Read in SU(3) matrices.
// the running of the lexicographic index and the su3_matrix layout in WUPP
// gauge configs is  the same as in QDP 
// (the fastest running direction is the 0th direction, corresponding to the 
// x-direction)

  u = zero ; 

  ColorMatrix  uu ; 
  ColorMatrixF  uuF ; 

  QDPIO::cout<<"Reading the gauge fields...\n" ;
  for(int site=0; site < Layout::vol(); ++site)
  {
     multi1d<int> coord = crtesn(site, Layout::lattSize()); // The coordinate
    // read in a single site
    for(int mu=0; mu < Nd; ++mu) 
      {  
	read(cfg_in, uuF );  
	if(byterev){
	  QDPUtil::byte_swap((void *)&uuF.elem(),sizeof(float),2*Nc*Nc);
	}
	uu = uuF;
	/*
	if(site < 10 && mu == 0){
	  QDPIO::cout << "site "<<site<<endl;
	  for(int ic = 0; ic < Nc; ic++){
	    for(int jc = 0; jc < Nc; jc++){
	      QDPIO::cout << "ic "<<ic<<" "<< jc<<" "<< peekColor(uuF,ic,jc)<<endl;
	      QDPIO::cout << "ic "<<ic<<" "<< jc<<" "<< peekColor(uu,ic,jc)<<endl;
	    }
	  }
	}
	*/
	pokeSite(u[mu],uu,coord); 
      }    

  }

  cfg_in.close();

  END_CODE();
}



//! Read a Wupp configuration file
/*!
 * \ingroup io
 *
 *
 * \param xml        xml reader holding config info ( Modify )
 * \param u          gauge configuration ( Modify )
 * \param cfg_file   path ( Read )
 */    

void readWupp(XMLReader& xml, multi1d<LatticeColorMatrix>& u, const string& cfg_file)
{
  START_CODE();

  // Read the config and its binary header
  readWupp(u, cfg_file);

#if 0
  // Now, set up the XML header. Do this by first making a buffer
  // writer that is then used to make the reader
  XMLBufferWriter  xml_buf;
  const string ss("wupper_config_read") ;
  write(xml, ss );


  write(xml_buf, ss );

  try 
  {
    xml.open(xml_buf);
  }
  catch(const string& e)
  { 
    QDPIO::cerr << __func__ << ": Error in readwupp: " << e.c_str() << endl;
    QDP_abort(1);
  }

#endif

  END_CODE();
}

}  // end namespace Chroma
