
/*! \file
 *  \brief Read a CERN gauge configuration 
 */

#include "chromabase.h"
#include "meas/glue/mesplq.h"
#include <vector>
using namespace std;

namespace Chroma
{


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


  
  
void localPokeCernLink(const multi1d<int>& coord, const double *buf,
                       LatticeColorMatrix& Umufield)
{
 int localsiteindex=Layout::linearSiteIndex(coord);
 ColorMatrix locColMat;
 int count=0;
 for (int row=0;row<QDP::Nc;row++)
 for (int col=0;col<QDP::Nc;col++){
    Complex sitecomp=cmplx(Real64(buf[count]),Real64(buf[count+1]));
    pokeColor(locColMat,sitecomp,row,col);
    count+=2;}
 Umufield.elem(localsiteindex) = locColMat.elem();
}



void pokeCernLinks(const double *rdbuf, int lcbz, int linkdbles, multi1d<int>& coord,
                   int NX, int NY, int NZ, int NT, multi1d<LatticeColorMatrix>& u, 
                   StopWatch& cmtimer)
{
    // reorder the sites in the buffer
 vector<double> cmbuf(8*lcbz*linkdbles);
 int count=0;
#if  defined(ARCH_PARSCALAR) || defined(ARCH_PARSCALARVEC)
 if (Layout::primaryNode()){
#endif
   for (int dir=0;dir<8;dir+=2)
   for (int lz=0;lz<lcbz;lz++){
      int indstart=(8*lz+dir)*linkdbles;
      for (int k=0;k<linkdbles;k++)
         cmbuf[count++]=rdbuf[indstart+k];}
   for (int dir=7;dir>0;dir-=2)
   for (int lz=lcbz-1;lz>=0;lz--){
      int indstart=(8*lz+dir)*linkdbles;
      for (int k=0;k<linkdbles;k++)
         cmbuf[count++]=rdbuf[indstart+k];}
#if  defined(ARCH_PARSCALAR) || defined(ARCH_PARSCALARVEC)
 }
#endif
   
  int zkeep,ykeep,xkeep,tkeep,xshift,yshift,zshift,shift;
           
#if  defined(ARCH_PARSCALAR) || defined(ARCH_PARSCALARVEC)

    //  get nodes involved
 int cnode=Layout::nodeNumber(coord);
 tkeep=coord[3]; 
 if (coord[3]==0) coord[3]=NT-1; else coord[3]--;
 int tmnode=Layout::nodeNumber(coord);
 coord[3]=tkeep;
 xkeep=coord[0];
 if (coord[0]==0) coord[0]=NX-1; else coord[0]--;
 int xmnode=Layout::nodeNumber(coord);
 coord[0]=xkeep;
 ykeep=coord[1];
 if (coord[1]==0) coord[1]=NY-1; else coord[1]--;
 int ymnode=Layout::nodeNumber(coord);
 coord[1]=ykeep;
 zkeep=coord[2];
 if (coord[2]==0) coord[2]=NZ-1; else coord[2]--;
 int zmnode=Layout::nodeNumber(coord);
 coord[2]=zkeep;
 
 cmtimer.start();
 size_t dbsize=sizeof(double);
 if ((cnode==tmnode)&&(cnode==xmnode)&&(cnode==ymnode)&&(cnode==zmnode)){
    size_t bytes=8*lcbz*linkdbles*dbsize;
    if (cnode!=0){
      if (Layout::primaryNode())
	QDPInternal::sendToWait((void *)&cmbuf[0],cnode,bytes);
      if (Layout::nodeNumber() == cnode)
	QDPInternal::recvFromWait((void *)&cmbuf[0],0,bytes);}}
 else{
    if (cnode!=0){
      size_t bytes=(5*lcbz-((cnode!=zmnode)?1:0))*linkdbles*dbsize;
      if (Layout::primaryNode())
	QDPInternal::sendToWait((void *)&cmbuf[0],cnode,bytes);
      if (Layout::nodeNumber() == cnode)
	QDPInternal::recvFromWait((void *)&cmbuf[0],0,bytes);}
    count=(5*lcbz-1)*linkdbles;
    if ((zmnode!=cnode)&&(zmnode!=0)){
      size_t bytes=linkdbles*dbsize;
      if (Layout::primaryNode())
	QDPInternal::sendToWait((void *)&cmbuf[count],zmnode,bytes);
      if (Layout::nodeNumber() == zmnode)
	QDPInternal::recvFromWait((void *)&cmbuf[count],0,bytes);}
    count=5*lcbz*linkdbles;
    size_t bytes=lcbz*linkdbles*dbsize;
    if (ymnode!=0){
      if (Layout::primaryNode())
	QDPInternal::sendToWait((void *)&cmbuf[count],ymnode,bytes);
      if (Layout::nodeNumber() == ymnode)
	QDPInternal::recvFromWait((void *)&cmbuf[count],0,bytes);}
    count=6*lcbz*linkdbles;
    if (xmnode!=0){
      if (Layout::primaryNode())
	QDPInternal::sendToWait((void *)&cmbuf[count],xmnode,bytes);
      if (Layout::nodeNumber() == xmnode)
	QDPInternal::recvFromWait((void *)&cmbuf[count],0,bytes);}
    count=7*lcbz*linkdbles;
    if (tmnode!=0){
      if (Layout::primaryNode())
	QDPInternal::sendToWait((void *)&cmbuf[count],tmnode,bytes);
      if (Layout::nodeNumber() == tmnode)
	QDPInternal::recvFromWait((void *)&cmbuf[count],0,bytes);}} 
 cmtimer.stop();
 
 if (Layout::nodeNumber()==cnode){
#endif
    zkeep=coord[2];
    xshift=lcbz*linkdbles;
    yshift=2*lcbz*linkdbles;
    zshift=3*lcbz*linkdbles;
    shift=0;
    for (int lz=0;lz<lcbz;lz++,shift+=linkdbles){
       localPokeCernLink(coord,&cmbuf[xshift+shift],u[0]);
       localPokeCernLink(coord,&cmbuf[yshift+shift],u[1]);
       localPokeCernLink(coord,&cmbuf[zshift+shift],u[2]);
       localPokeCernLink(coord,&cmbuf[shift],u[3]);
       coord[2]+=2;}
    coord[2]=zkeep+2*(lcbz-1)-1;
    shift=4*lcbz*linkdbles;
    for (int lz=1;lz<lcbz;lz++,shift+=linkdbles){
       localPokeCernLink(coord,&cmbuf[shift],u[2]);
       coord[2]-=2;} 
    coord[2]=zkeep;
#if  defined(ARCH_PARSCALAR) || defined(ARCH_PARSCALARVEC)
    }  

 if (Layout::nodeNumber()==zmnode){
#endif
    zkeep=coord[2];
    if (coord[2]==0) coord[2]=NZ-1; else coord[2]--;
    shift=(5*lcbz-1)*linkdbles;
    localPokeCernLink(coord,&cmbuf[shift],u[2]);
    coord[2]=zkeep;
#if  defined(ARCH_PARSCALAR) || defined(ARCH_PARSCALARVEC)
    }

 if (Layout::nodeNumber()==ymnode){
#endif
    zkeep=coord[2];
    ykeep=coord[1];
    coord[2]=zkeep+2*(lcbz-1);
    if (coord[1]==0) coord[1]=NY-1; else coord[1]--;
    shift=5*lcbz*linkdbles;
    for (int lz=0;lz<lcbz;lz++,shift+=linkdbles){
       localPokeCernLink(coord,&cmbuf[shift],u[1]);
       coord[2]-=2;}
    coord[1]=ykeep;
    coord[2]=zkeep;
#if  defined(ARCH_PARSCALAR) || defined(ARCH_PARSCALARVEC)
    }

 if (Layout::nodeNumber()==xmnode){
#endif
    zkeep=coord[2];
    xkeep=coord[0];
    coord[2]=zkeep+2*(lcbz-1);
    if (coord[0]==0) coord[0]=NX-1; else coord[0]--;
    shift=6*lcbz*linkdbles;
    for (int lz=0;lz<lcbz;lz++,shift+=linkdbles){
       localPokeCernLink(coord,&cmbuf[shift],u[0]);
       coord[2]-=2;}
    coord[0]=xkeep;
    coord[2]=zkeep;
#if  defined(ARCH_PARSCALAR) || defined(ARCH_PARSCALARVEC)
    }

 if (Layout::nodeNumber()==tmnode){
#endif
    zkeep=coord[2];
    tkeep=coord[3];
    coord[2]=zkeep+2*(lcbz-1);
    if (coord[3]==0) coord[3]=NT-1; else coord[3]--;
    shift=7*lcbz*linkdbles;
    for (int lz=0;lz<lcbz;lz++,shift+=linkdbles){
       localPokeCernLink(coord,&cmbuf[shift],u[3]);
       coord[2]-=2;}
    coord[3]=tkeep;
    coord[2]=zkeep; 
#if  defined(ARCH_PARSCALAR) || defined(ARCH_PARSCALARVEC)
    }
 
#endif
}




void readCERN(multi1d<LatticeColorMatrix>& u, const std::string& cfg_file)
{
 if ((sizeof(int)!=4)||(sizeof(double)!=8)){
   QDPIO::cout << "CERN files contain 4-byte ints, 8-byte doubles"<<std::endl;
   QDP_abort(1);}
 if (QDP::Nc!=3){
   QDPIO::cout << "readCERN only supports Nc=3"<<std::endl;
   QDP_abort(1);}
 if (QDP::Nd!=4){
   QDPIO::cout << "readCERN only supported for 4 space-time dimensions"<<std::endl;
   QDP_abort(1);}

 int NT = QDP::Layout::lattSize()[3];
 int NZ = QDP::Layout::lattSize()[2];
 int NY = QDP::Layout::lattSize()[1];
 int NX = QDP::Layout::lattSize()[0];
 if (NZ%2){
    QDPIO::cout << "NZ must be even"<<endl;
    QDP_abort(1); }
   
 QDPIO::cout <<std::endl<< "Beginning read of CERN gauge field ("<<NX<<" x "<<NY
             <<" x "<<NZ<<" ) x "<<NT<<std::endl;
 StopWatch rtimer; rtimer.start(); 
 StopWatch iotimer;
 StopWatch cmtimer;
  
 //  CERN field always stored as little endian
 iotimer.start();
 BinaryFileReader fin(cfg_file);
 int Ncern[4]; 
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

 multi1d<int> gridsize=Layout::subgridLattSize();
 int LZ=gridsize[2];
 size_t dbsize=sizeof(double);
 size_t linkdbles=2*QDP::Nc*QDP::Nc;
 size_t rdlinks=4*NZ;   // This is 8 directions *  (NZ/2) odd sites in the z-direction; NZ must be even
 size_t rdsize=rdlinks*linkdbles;
 vector<double> rdbuf(rdsize);
 
 multi1d<int> coord(QDP::Nd);
 for (int it=0;it<NT;it++){
   coord[3]=it;
   for (int ix=0;ix<NX;ix++){
     coord[0]=ix;
     for (int iy=0;iy<NY;iy++){
       coord[1]=iy;
	    // read entire z on primary
       iotimer.start();
       fin.readArrayPrimaryNodeLittleEndian((char*)&rdbuf[0],dbsize,rdsize);
       iotimer.stop();

       int izstart=((it+ix+iy)%2)?0:1;
       int zshift=0;
       int zz=0;
       if (LZ%2){   // careful if LZ is odd
          zz=(izstart)?1:-1;
          }
       for (int iz=izstart;iz<NZ;iz+=LZ+zz){
          int lcbz=(LZ-zz)/2;
          coord[2]=iz;
          pokeCernLinks(&rdbuf[zshift],lcbz,linkdbles,coord,NX,NY,NZ,NT,u,cmtimer);
          zshift+=8*lcbz*linkdbles;
          if (zz!=0) zz=-zz;
          }}}}
 
 rtimer.stop();
 QDPIO::cout << "readCERN:       plaq read: " << plaq << std::endl;
 Double w_plaq,s_plaq,t_plaq,link;
 MesPlq(u, w_plaq, s_plaq, t_plaq,link);  
 QDPIO::cout << "readCERN: plaq recomputed: " << w_plaq << std::endl;
 QDPIO::cout << "Read of CERN gauge field done in "<<rtimer.getTimeInSeconds()<<" seconds"<<std::endl;
 QDPIO::cout << "Time of IO operations = "<<iotimer.getTimeInSeconds()<<" seconds"<<std::endl;
 QDPIO::cout << "Time of communications = "<<cmtimer.getTimeInSeconds()<<" seconds"<<std::endl<<std::endl;
#if  defined(ARCH_PARSCALAR) || defined(ARCH_PARSCALARVEC)
 QMP_barrier();
#endif
}


}  // end namespace Chroma

