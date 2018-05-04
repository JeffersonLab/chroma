
/*! \file
 *  \brief Read a CERN gauge configuration 
 */

#include "chromabase.h"
#include "meas/glue/mesplq.h"

namespace Chroma
{
  //  SU3 matrix is row major (column fastest) in CERN format

  void pokeCernLink(const double *linkbuffer, size_t bytes, multi1d<int> sitecoord,
		    LatticeColorMatrix& Umufield)
  { 
    ColorMatrix U;
#if  defined(ARCH_PARSCALAR) || defined(ARCH_PARSCALARVEC)
    int recvnode=Layout::nodeNumber(sitecoord);
    if (recvnode!=0){
      if (Layout::primaryNode())
	QDPInternal::sendToWait((void *)linkbuffer,recvnode,bytes);
      if (Layout::nodeNumber() == recvnode)
	QDPInternal::recvFromWait((void *)linkbuffer,0,bytes);}
    if (Layout::nodeNumber() == recvnode){
#endif
      int count=0;
      for (int row=0;row<Nc;row++)
	for (int col=0;col<Nc;col++){
	  Complex sitecomp=cmplx(Real64(linkbuffer[count]),Real64(linkbuffer[count+1]));
	  pokeColor(U,sitecomp,row,col);
	  count+=2;}
#if  defined(ARCH_PARSCALAR) || defined(ARCH_PARSCALARVEC)
    }
#endif
    pokeSite(Umufield,U,sitecoord);
  }



  //  Read gauge field in CERN format into "u".
  //  Details about CERN format:
  //     - little endian, ints are 4 bytes, doubles are 8 bytes
  //     - 4 ints NT,NX,NY,NZ
  //     - 1 double for average plaquette
  //     - links as SU3 matrices (row major, double precision complex)
  //         in the following order:
  //            - the 8 links in directions +0,-0,...,+3,-3 at the first 
  //              odd point, the second odd point, and so on. 
  //            - The order of the point (x0,x1,x2,x3) with
  //               Cartesian coordinates in the range 
  //                 0<=x0<N0,...,0<=x3<N3 is determined by the index
  //                   ix=x3+N3*x2+N2*N3*x1+N1*N2*N3*x0,
  //               where N0,N1,N2,N3 are the global lattice sizes


  void readCERN(multi1d<LatticeColorMatrix>& u, const std::string& cfg_file)
  {
    if ((sizeof(int)!=4)||(sizeof(double)!=8)){
      QDPIO::cout << "CERN files contain 4-byte ints, 8-byte doubles"<<std::endl;
      QDP_abort(1);}
    if (QDP::Nd!=4){
      QDPIO::cout << "readCERN only supported for 4 space-time dimensions"<<std::endl;
      QDP_abort(1);}

    int NT = QDP::Layout::lattSize()[3];
    int NZ = QDP::Layout::lattSize()[2];
    int NY = QDP::Layout::lattSize()[1];
    int NX = QDP::Layout::lattSize()[0];
    QDPIO::cout <<std::endl<< "Beginning read of CERN gauge field ("<<NX<<" x "<<NY
		<<" x "<<NZ<<" ) x "<<NT<<std::endl;
    StopWatch rtimer; rtimer.start(); 
    StopWatch iotimer;
 
    //  CERN field always stored as little endian
    iotimer.start();
    BinaryFileReader fin(cfg_file);
    uint Ncern[4]; 
    fin.readArrayLittleEndian((char*)&Ncern[0],sizeof(int),4);  // assumes little endian, no checksums
    iotimer.stop();
  
    if ((Ncern[0]!=NT)||(Ncern[1]!=NX)||(Ncern[2]!=NY)||(Ncern[3]!=NZ)){
      QDPIO::cout << "readCERN: Lattice size mismatch " << std::endl;
      QDPIO::cout << "read "<<Ncern[1]<<" "<<Ncern[2]<<" "<<Ncern[3]
		  <<" "<<Ncern[0]<<std::endl;
      QDPIO::cout << "Chroma wants "<<NX<<" "<<NY<<" "<<NZ<<" "<<NT<<std::endl;
      QDP_abort(1); }
    double plaq;
    iotimer.start();
    fin.readArrayLittleEndian((char*)&plaq,sizeof(double),1);
    iotimer.stop();
    plaq/=3.0;

    size_t linkdbles=2*QDP::Nc*QDP::Nc;
    size_t ndir=2*QDP::Nd;
    size_t nelem=ndir*linkdbles;
    size_t dbsize=sizeof(double);
    size_t linkbytes=dbsize*linkdbles;
    double sitebuffer[nelem];
    for (int it=0;it<NT;it++)
      for (int ix=0;ix<NX;ix++)
	for (int iy=0;iy<NY;iy++)
	  for (int iz=0;iz<NZ;iz++)
	    if ((ix+iy+iz+it)%2){

	      multi1d<int> coord(QDP::Nd);
	      coord[0]=ix; coord[1]=iy; coord[2]=iz; coord[3]=it;
	      iotimer.start();
	      // reads on primary
	      fin.readArrayPrimaryNodeLittleEndian((char*)&sitebuffer[0],dbsize,nelem);
	      iotimer.stop();

	      pokeCernLink(&sitebuffer[0],linkbytes,coord,u[3]);
	      pokeCernLink(&sitebuffer[2*linkdbles],linkbytes,coord,u[0]);
	      pokeCernLink(&sitebuffer[4*linkdbles],linkbytes,coord,u[1]);
	      pokeCernLink(&sitebuffer[6*linkdbles],linkbytes,coord,u[2]);

	      int keep=coord[3]; 
	      if (coord[3]==0) coord[3]=NT-1; else coord[3]--;
	      pokeCernLink(&sitebuffer[linkdbles],linkbytes,coord,u[3]);
	      coord[3]=keep;
	      keep=coord[0];
	      if (coord[0]==0) coord[0]=NX-1; else coord[0]--;
	      pokeCernLink(&sitebuffer[3*linkdbles],linkbytes,coord,u[0]);
	      coord[0]=keep;
	      keep=coord[1];
	      if (coord[1]==0) coord[1]=NY-1; else coord[1]--;
	      pokeCernLink(&sitebuffer[5*linkdbles],linkbytes,coord,u[1]);
	      coord[1]=keep;
	      if (coord[2]==0) coord[2]=NZ-1; else coord[2]--;
	      pokeCernLink(&sitebuffer[7*linkdbles],linkbytes,coord,u[2]);
	    }

    QDPIO::cout << "readCERN: plaq read: " << plaq << std::endl;
    Double w_plaq,s_plaq,t_plaq,link;
    MesPlq(u, w_plaq, s_plaq, t_plaq,link);  
    QDPIO::cout << "readCERN: plaq recomputed " << w_plaq << std::endl;
    rtimer.stop();
    QDPIO::cout << "Read of CERN gauge field done in "<<rtimer.getTimeInSeconds()<<" seconds"<<std::endl;
    QDPIO::cout << "Time of IO operations = "<<iotimer.getTimeInSeconds()<<" seconds"<<std::endl<<std::endl;
  }



}  // end namespace Chroma

