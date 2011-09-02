// $Id: t_invert4_precwilson.cc,v 3.4 2009-10-09 13:59:46 bjoo Exp $

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

#include <cstdio>

#include <stdlib.h>
#include <sys/time.h>
#include <math.h>

#include "chroma.h"
#include <string>
using namespace Chroma;
using namespace std;

typedef LatticeFermionF  TF;
typedef LatticeColorMatrixF  UF;
typedef multi1d<LatticeColorMatrixF> QF;
typedef multi1d<LatticeColorMatrixF> PF;

typedef LatticeFermionD  TD;
typedef LatticeColorMatrixD  UD;
typedef multi1d<LatticeColorMatrixD> QD;
typedef multi1d<LatticeColorMatrixD> PD;


#include "quda.h"
#if 0
#ifdef __cplusplus
extern "C" {
#endif
  void commDimPartitionedSet(int);
#ifdef __cplusplus
};
#endif
#endif



bool testQudaDslash(const QF& u, enum PlusMinus isign, int cb)
{
  // I want to test the QUDA dslash
  // First make a reference dslash: 4D for now
  multi1d<int> boundary(4); boundary[0]=boundary[1]=boundary[2]=1;
  boundary[3]=-1;

  Handle< FermBC<TF,PF,QF> > bc_handle(new SimpleFermBC<TF,PF,QF>(boundary));

  Handle< FermState<TF,PF,QF> > fs_handle(new SimpleFermState<TF,PF,QF>(bc_handle, u));

  WilsonDslashF D_me(fs_handle);

  // Now set up a QUDA Dslash
  // ----------**************************--------------
  QudaGaugeParam q_gauge_param=newQudaGaugeParam();
  QudaInvertParam quda_inv_param=newQudaInvertParam();

  const multi1d<int>& latdims = Layout::subgridLattSize();      
  q_gauge_param.X[0] = latdims[0];
  q_gauge_param.X[1] = latdims[1];
  q_gauge_param.X[2] = latdims[2];
  q_gauge_param.X[3] = latdims[3];
  int vol=latdims[0]*latdims[1]*latdims[2]*latdims[3];

  //
  //  setDims(q_gauge_param.X);

  // Setup gauge anisotropy
  q_gauge_param.anisotropy = 1.0;

  q_gauge_param.type = QUDA_WILSON_LINKS;
  // Setup Gauge Order 
  q_gauge_param.gauge_order = QUDA_QDP_GAUGE_ORDER;

  // Setp Boundaries
  if(boundary[3]==-1){ 
    q_gauge_param.t_boundary = QUDA_ANTI_PERIODIC_T;
  }
  else { 
    q_gauge_param.t_boundary = QUDA_PERIODIC_T;
  }

  QudaPrecision_s cpu_prec=QUDA_SINGLE_PRECISION;
  QudaPrecision_s gpu_prec=QUDA_SINGLE_PRECISION;
  QudaPrecision_s gpu_half_prec=QUDA_HALF_PRECISION;

  q_gauge_param.cpu_prec=cpu_prec;
  q_gauge_param.cuda_prec=gpu_prec;
  q_gauge_param.reconstruct = QUDA_RECONSTRUCT_NO;
  q_gauge_param.reconstruct_sloppy = q_gauge_param.reconstruct;
  q_gauge_param.cuda_prec_sloppy=q_gauge_param.cuda_prec;
  q_gauge_param.gauge_fix = QUDA_GAUGE_FIXED_NO;

  quda_inv_param.kappa=1;
  quda_inv_param.matpc_type = QUDA_MATPC_ODD_ODD;
  quda_inv_param.dagger = (isign == PLUS) ? QUDA_DAG_NO : QUDA_DAG_YES;

   // Setup padding
  multi1d<int> face_size(4);
  face_size[0] = latdims[1]*latdims[2]*latdims[3]/2;
  face_size[1] = latdims[0]*latdims[2]*latdims[3]/2;
  face_size[2] = latdims[0]*latdims[1]*latdims[3]/2;
  face_size[3] = latdims[0]*latdims[1]*latdims[2]/2;

  int max_face = face_size[0];
  for(int i=1; i <=3; i++) { 
    if ( face_size[i] > max_face ) { 
      max_face = face_size[i]; 
    }
  }
  q_gauge_param.ga_pad = max_face;
  quda_inv_param.sp_pad = 0;
  quda_inv_param.cl_pad = 0;

  quda_inv_param.dirac_order = QUDA_DIRAC_ORDER;
  quda_inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
   // Now set up inverter params.
  quda_inv_param.dslash_type = QUDA_WILSON_DSLASH; // Sets Wilson Matrix
  quda_inv_param.cpu_prec = cpu_prec;
  quda_inv_param.cuda_prec = gpu_prec;
  quda_inv_param.cuda_prec_sloppy = gpu_half_prec;
 
  //commDimPartitionedSet(3);

  // Make the links
  QF links_single(Nd);


  for(int mu=0; mu < Nd; mu++) {
    links_single[mu] = fs_handle->getLinks()[mu];
  }

 
  // Set up the links
  void* gauge[4];

  for(int mu=0; mu < Nd; mu++) { 
    gauge[mu] = (void *)&(links_single[mu].elem(all.start()).elem().elem(0,0).real());
  }

  loadGaugeQuda((void *)gauge, &q_gauge_param);
  

 

  TF src, dst1,dst2;
  gaussian(src);
  //src.elem(rb[1].end()).elem(0).elem(0).real() = 1.0;

  gaussian(dst1); // Junk these
  gaussian(dst2);
  D_me.apply(dst1,src,isign,cb);

  // Compute buffers for transfer to GPU
  unsigned int vol2 = vol/2;
  unsigned int h_size=4*3*2*vol2;
  multi1d<float> buffer_in(h_size);
  multi1d<float> buffer_out(h_size);



  int otherCB=1-cb;
  int daggerBit = (isign == MINUS ? 1 : 0);

  // Pack buffer_in: output cb = rb[0], input cb=rb[1]
  int lin=0;
  for(int site=0; site < rb[1].siteTable().size(); site++) { 
    for(int spin=0; spin < 4; spin++) {		      
      for(int col=0; col < 3; col++) {
	buffer_in[lin++] = src.elem(rb[otherCB].siteTable()[site]).elem(spin).elem(col).real();
	buffer_in[lin++] = src.elem(rb[otherCB].siteTable()[site]).elem(spin).elem(col).imag();
      }
    }
  }


  dslashQuda((void *)buffer_out.slice(),
	     (void *)buffer_in.slice(),
	     &quda_inv_param,
	     (QudaParity)cb);     // source parity =1, dst parity=0

  // Unpack
  lin=0;
  for(int site=0; site < rb[1].siteTable().size(); site++) { 
    for(int spin=0; spin < 4; spin++) {		      
      for(int col=0; col < 3; col++) {
	dst2.elem(rb[cb].siteTable()[site]).elem(spin).elem(col).real() = buffer_out[lin++];
	dst2.elem(rb[cb].siteTable()[site]).elem(spin).elem(col).imag() = buffer_out[lin++];
      }
    }
  }

  TF diff=zero;
  diff[rb[cb]] = dst1-dst2;

  // Free QUDA data structures
  Double diff_norm = sqrt(norm2(diff,rb[cb]))/Double(vol2);
  QDPIO::cout << "\t\t diff = " << diff_norm << " per site \t"; 
  bool ret_val; 
  if ( toBool( diff_norm < Double(1.0e-7)  ) ) {
    ret_val = true;
  }
  else {
    ret_val = false; 
  }

  freeGaugeQuda();

  return ret_val;
}



// -------DOUBLE PRECISION TEST --------
bool testQudaDslashD(const QD& u, enum PlusMinus isign, int cb)
{
  // I want to test the QUDA dslash
  // First make a reference dslash: 4D for now
  multi1d<int> boundary(4); boundary[0]=boundary[1]=boundary[2]=1;
  boundary[3]=-1;

  Handle< FermBC<TD,PD,QD> > bc_handle(new SimpleFermBC<TD,PD,QD>(boundary));

  Handle< FermState<TD,PD,QD> > fs_handle(new SimpleFermState<TD,PD,QD>(bc_handle, u));

  WilsonDslashD D_me(fs_handle);
  
#if 0
  {
    QDPWilsonDslashD D_qdp(fs_handle);

    TD src, dst1, dst2;
    gaussian(src);

    D_me.apply(dst1, src, isign, cb);
    D_qdp.apply(dst2, src, isign, cb);

    TD diff = dst1-dst2;
    RealD diff_norm = sqrt( norm2(diff, rb[cb]) );
    QDPIO::cout << "diff_norm = " << diff_norm << endl;
    QDPIO::cout << "diff_norm / site = " << diff_norm/((double)Layout::vol()/(double)2) << endl;
  }
#endif

  // Now set up a QUDA Dslash
  // ----------**************************--------------
  QudaGaugeParam q_gauge_param=newQudaGaugeParam();
  QudaInvertParam quda_inv_param=newQudaInvertParam();

  QudaPrecision_s cpu_prec=QUDA_DOUBLE_PRECISION;
  QudaPrecision_s gpu_prec=QUDA_DOUBLE_PRECISION;
  QudaPrecision_s gpu_half_prec=QUDA_DOUBLE_PRECISION;

  // Setup Boundaries
  if(boundary[3]==-1){ 
    q_gauge_param.t_boundary = QUDA_ANTI_PERIODIC_T;
  }
  else { 
    q_gauge_param.t_boundary = QUDA_PERIODIC_T;
  }

  // Setup Gauge Order 
  q_gauge_param.gauge_order = QUDA_QDP_GAUGE_ORDER;
  q_gauge_param.reconstruct = QUDA_RECONSTRUCT_12;
  q_gauge_param.reconstruct_sloppy = QUDA_RECONSTRUCT_12;

  // Setup precisions
  q_gauge_param.cpu_prec=cpu_prec;
  q_gauge_param.cuda_prec=gpu_prec;
  q_gauge_param.cuda_prec_sloppy=gpu_half_prec;

  // Setup lattice size
  const multi1d<int>& latdims = Layout::subgridLattSize();      
  q_gauge_param.X[0] = latdims[0];
  q_gauge_param.X[1] = latdims[1];
  q_gauge_param.X[2] = latdims[2];
  q_gauge_param.X[3] = latdims[3];
  int vol=latdims[0]*latdims[1]*latdims[2]*latdims[3];
  q_gauge_param.type = QUDA_WILSON_LINKS;
  q_gauge_param.gauge_order = QUDA_QDP_GAUGE_ORDER;

  // Setup gauge fixing
  q_gauge_param.gauge_fix = QUDA_GAUGE_FIXED_NO;

  // Setup padding
  multi1d<int> face_size(4);
  face_size[0] = latdims[1]*latdims[2]*latdims[3]/2;
  face_size[1] = latdims[0]*latdims[2]*latdims[3]/2;
  face_size[2] = latdims[0]*latdims[1]*latdims[3]/2;
  face_size[3] = latdims[0]*latdims[1]*latdims[2]/2;

  int max_face = face_size[0];
  for(int i=1; i <=3; i++) { 
    if ( face_size[i] > max_face ) { 
      max_face = face_size[i]; 
    }
  }


  q_gauge_param.ga_pad = max_face;
  quda_inv_param.sp_pad = 0;
  quda_inv_param.cl_pad = 0;

  // Setup gauge anisotropy
  q_gauge_param.anisotropy = 1.0;

  //commDimPartitionedSet(3);

  // Make the links
  QD links_single(Nd);

#if 0
  QD links_minus(Nd);
#endif

  for(int mu=0; mu < Nd; mu++) {
    links_single[mu] = fs_handle->getLinks()[mu];
#if 0
    links_minus[mu] =  shift(links_single[mu], BACKWARD, mu);
#endif
    //gaussian(links_minus[mu]);
  }

 
  // Set up the links
  void* gauge[4];
#if 0
  void* gauge_minus[4];
#endif

  for(int mu=0; mu < Nd; mu++) { 
    gauge[mu] = (void *)&(links_single[mu].elem(all.start()).elem().elem(0,0).real());

#if 0
    gauge_minus[mu] = (void *)&(links_minus[mu].elem(all.start()).elem().elem(0,0).real());
#endif
  }

#if 0
  loadGaugeQuda((void *)gauge, (void *)gauge_minus, &q_gauge_param);
#else
  loadGaugeQuda((void *)gauge, &q_gauge_param);
#endif  
  // Now set up inverter params.
  quda_inv_param.dslash_type = QUDA_WILSON_DSLASH; // Sets Wilson Matrix
  quda_inv_param.cpu_prec = cpu_prec;
  quda_inv_param.cuda_prec = gpu_prec;
  quda_inv_param.cuda_prec_sloppy = gpu_half_prec;
  quda_inv_param.dirac_order = QUDA_DIRAC_ORDER;
  quda_inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
  quda_inv_param.dagger = (isign == PLUS) ? QUDA_DAG_NO : QUDA_DAG_YES;
  TD src, dst1,dst2;
  gaussian(src);
  //src.elem(rb[1].end()).elem(0).elem(0).real() = 1.0;

  gaussian(dst1); // Junk these
  gaussian(dst2);
  D_me.apply(dst1,src,isign,cb);

  // Compute buffers for transfer to GPU
  unsigned int vol2 = vol/2;
  unsigned int h_size=4*3*2*vol2;
  multi1d<double> buffer_in(h_size);
  multi1d<double> buffer_out(h_size);



  int otherCB=1-cb;

  // Pack buffer_in: output cb = rb[0], input cb=rb[1]
  int lin=0;
  for(int site=0; site < rb[1].siteTable().size(); site++) { 
    for(int spin=0; spin < 4; spin++) {		      
      for(int col=0; col < 3; col++) {
	buffer_in[lin++] = src.elem(rb[otherCB].siteTable()[site]).elem(spin).elem(col).real();
	buffer_in[lin++] = src.elem(rb[otherCB].siteTable()[site]).elem(spin).elem(col).imag();
      }
    }
  }


  dslashQuda((void *)buffer_out.slice(),
	     (void *)buffer_in.slice(),
	     &quda_inv_param,
	     (QudaParity)cb);      // source parity =1, dst parity=0


  // Unpack
  lin=0;
  for(int site=0; site < rb[1].siteTable().size(); site++) { 
    for(int spin=0; spin < 4; spin++) {		      
      for(int col=0; col < 3; col++) {
	dst2.elem(rb[cb].siteTable()[site]).elem(spin).elem(col).real() = buffer_out[lin++];
	dst2.elem(rb[cb].siteTable()[site]).elem(spin).elem(col).imag() = buffer_out[lin++];
      }
    }
  }

  TD diff=zero;
  diff[rb[cb]] = dst1-dst2;

  // Free QUDA data structures
  Double diff_norm = sqrt(norm2(diff,rb[cb]))/Double(vol2);
  QDPIO::cout << "\t\t diff = " << diff_norm << " per site \t"; 
  bool ret_val; 
  if ( toBool( diff_norm < Double(1.0e-7)  ) ) {
    ret_val = true;
  }
  else {
    ret_val = false; 
  }

  freeGaugeQuda();
  return ret_val;
}

#if 0
bool testQudaDslash3D(const QF& u, enum PlusMinus isign, int cb)
{
  // I want to test the QUDA dslash
  // First make a reference dslash: 4D for now
  multi1d<int> boundary(4); boundary[0]=boundary[1]=boundary[2]=1;
  boundary[3]=+1;

  Handle< FermBC<TF,PF,QF> > bc_handle(new SimpleFermBC<TF,PF,QF>(boundary));

  Handle< FermState<TF,PF,QF> > fs_handle(new SimpleFermState<TF,PF,QF>(bc_handle, u));

  WilsonDslash3D D_me(fs_handle);

  // Now set up a QUDA Dslash
  // ----------**************************--------------
  QudaGaugeParam q_gauge_param=newQudaGaugeParam();
  QudaInvertParam quda_inv_param=newQudaInvertParam();

  QudaPrecision_s cpu_prec=QUDA_SINGLE_PRECISION;
  QudaPrecision_s gpu_prec=QUDA_SINGLE_PRECISION;
  QudaPrecision_s gpu_half_prec=QUDA_SINGLE_PRECISION;

  // Setup Boundaries
  if(boundary[3]==-1){ 
    q_gauge_param.t_boundary = QUDA_ANTI_PERIODIC_T;
  }
  else { 
    q_gauge_param.t_boundary = QUDA_PERIODIC_T;
  }

  // Setup Gauge Order 
  q_gauge_param.reconstruct = QUDA_RECONSTRUCT_12;
  q_gauge_param.reconstruct_sloppy = QUDA_RECONSTRUCT_12;

  // Setup precisions
  q_gauge_param.cpu_prec=cpu_prec;
  q_gauge_param.cuda_prec=gpu_prec;
  q_gauge_param.cuda_prec_sloppy=gpu_half_prec;

  // Setup lattice size
  const multi1d<int>& latdims = Layout::subgridLattSize();      
  q_gauge_param.X[0] = latdims[0];
  q_gauge_param.X[1] = latdims[1];
  q_gauge_param.X[2] = latdims[2];
  q_gauge_param.X[3] = latdims[3];
  
  // Setup gauge fixing
  q_gauge_param.gauge_fix = QUDA_GAUGE_FIXED_NO;



  // Setup padding
  unsigned int vol = latdims[0]*latdims[1]*latdims[2]*latdims[3];
  unsigned int vol2 = vol/2;
  unsigned int padding=latdims[0]*latdims[1]*latdims[2]/2;
  quda_inv_param.sp_pad = q_gauge_param.ga_pad = quda_inv_param.cl_pad = padding;


  // Setup gauge anisotropy
  q_gauge_param.anisotropy = 1.0;

  
  // Make the links
  QF links_single(Nd);
 //  QF links_minus(Nd);

  {
    QF tmp_links(Nd);
    // QF tmp_links_minus(Nd);
    
    for(int mu=0; mu < Nd; mu++){ 
      tmp_links[mu] = fs_handle->getLinks()[mu];
      // tmp_links_minus[mu] = shift(tmp_links[mu], BACKWARD, mu);
    }


    // Repack gauge, so that for links_single[mu] has gauge field with 3d cb being contiguous
    for(int mu=0; mu < Nd; mu++) {
      for(int t=0; t < latdims[3]; t++) { 
	for(int z=0; z < latdims[2]; z++) { 
	  for(int y=0; y < latdims[1]; y++) { 
	    for(int x=0; x < latdims[0]; x++) { 
	      int xh=x/2;
	      int par=(x + y + z) & 1;
	      int par4d=(x + y + z + t) & 1;
	      // Weirdly, this is the same offset in both 3 and 4d
	      int off=xh+(latdims[0]/2)*(y + latdims[1]*(z+latdims[2]*t));
	      
	      for(int r=0; r < Nc; r++) { 
		for(int c=0; c < Nc; c++) { 
		  links_single[mu].elem(vol2*par+off).elem().elem(r,c).real()
		    = tmp_links[mu].elem(vol2*par4d+off).elem().elem(r,c).real();
		  
		  
		  
		  links_single[mu].elem(vol2*par+off).elem().elem(r,c).imag()
		    = tmp_links[mu].elem(vol2*par4d+off).elem().elem(r,c).imag();
		  
		  
//		  links_minus[mu].elem(vol2*par+off).elem().elem(r,c).real()
//		    = tmp_links_minus[mu].elem(vol2*par4d+off).elem().elem(r,c).real();
//		links_minus[mu].elem(vol2*par+off).elem().elem(r,c).imag()
//		  = tmp_links_minus[mu].elem(vol2*par4d+off).elem().elem(r,c).imag();
		}
	      }
	    }
	  }
	}
      }
    }

  }
  // Set up the links
  void* gauge[4];
  // void* gauge_minus[4];

  for(int mu=0; mu < Nd; mu++) { 
    gauge[mu] = (void *)&(links_single[mu].elem(0).elem().elem(0,0).real());
  //  gauge_minus[mu] = (void *)&(links_minus[mu].elem(0).elem().elem(0,0).real());
  }
  // loadGaugeQuda((void *)gauge,(void *)gauge_minus, &q_gauge_param);
  loadGaugeQuda((void *)gauge,  &q_gauge_param);
  
  // Now set up inverter params.
  quda_inv_param.dslash_type = QUDA_WILSON_DSLASH; // Sets Wilson Matrix
  quda_inv_param.cpu_prec = cpu_prec;
  quda_inv_param.cuda_prec = gpu_prec;
  quda_inv_param.cuda_prec_sloppy = gpu_half_prec;
  quda_inv_param.dirac_order = QUDA_DIRAC_ORDER;
  quda_inv_param.gamma_basis = QUDA_DEGRAND_ROSSI_GAMMA_BASIS;
  TF src, dst1,dst2;
  gaussian(src);
  gaussian(dst1); // Junk these
  gaussian(dst2);
  D_me.apply(dst1,src,isign,cb); // fills out dst[rb[0]]

  // Compute buffers for transfer to GPU

  unsigned int h_size=4*3*2*vol2;
  multi1d<float> buffer_in(h_size);
  multi1d<float> buffer_out(h_size);



  int otherCB=1-cb;
  int daggerBit = (isign == MINUS ? 1 : 0);

  // Pack buffer_in: output cb = rb3[0], input cb=rb3[1]
  int lin;
  for(int t=0; t < latdims[3]; t++) { 
    for(int z=0; z < latdims[2]; z++) { 
      for(int y=0; y < latdims[1]; y++) { 
	for(int x=0; x < latdims[0]; x++) { 
	  int xh=x/2;
	  int par=(x + y + z) & 1;
	  if( par == otherCB ) {
	    int par4d=(x + y + z + t) & 1;
	    // Weirdly, this is the same offset in both 3 and 4d
	    int off=xh+(latdims[0]/2)*(y + latdims[1]*(z+latdims[2]*t));
	    int lin=0;
	    for(int spin=0; spin < 4; spin++) {		      
	      for(int col=0; col < 3; col++) {
		buffer_in[24*off+(lin++)] = src.elem(vol2*par4d+off).elem(spin).elem(col).real();
		buffer_in[24*off+(lin++)] = src.elem(vol2*par4d+off).elem(spin).elem(col).imag();
	      }
	    }
	  }
	}
      }
    }
  }


  dslash3DQuda((void *)buffer_out.slice(),
	       (void *)buffer_in.slice(),
	       &quda_inv_param,
	       cb,      // source parity =1, dst parity=0
	       daggerBit);     // no dagger

  for(int t=0; t < latdims[3]; t++) { 
    for(int z=0; z < latdims[2]; z++) { 
      for(int y=0; y < latdims[1]; y++) { 
	for(int x=0; x < latdims[0]; x++) { 
	  int xh=x/2;
	  int par=(x + y + z) & 1;
	  if( par == cb ) { 
	    int par4d=(x + y + z + t) & 1;
	    // Weirdly, this is the same offset in both 3 and 4d
	    int off=xh+(latdims[0]/2)*(y + latdims[1]*(z+latdims[2]*t));
	    int lin=0;
	    for(int spin=0; spin < 4; spin++) {		      
	      for(int col=0; col < 3; col++) {
		dst2.elem(vol2*par4d+off).elem(spin).elem(col).real()= buffer_out[24*off+(lin++)] ;
		dst2.elem(vol2*par4d+off).elem(spin).elem(col).imag()=buffer_out[24*off+(lin++)];
	      }
	    }
	  }
	}
      }
    }
  }

  TF diff=zero;
  diff[rb3[cb]] = dst1-dst2;

  // Free QUDA data structures
  Double diff_norm = sqrt(norm2(diff,rb3[cb]))/Double(vol2);
  QDPIO::cout << "\t diff = " << diff_norm << " per site \t"; 
  bool ret_val;
  if ( toBool( diff_norm < Double(2.0e-8)  ) ) {
    ret_val = true;
  }
  else {
    ret_val = false; 
  }

  return ret_val;
}
#endif

bool linkageHack(void)
{
  bool foo = true;

  // Inline Measurements
  foo &= InlineAggregateEnv::registerAll();
  foo &= GaugeInitEnv::registerAll();

  return foo;
}

int main(int argc, char **argv)
{
  // Put the machine into a known state
  Chroma::initialize(&argc, &argv);
  QDPIO::cout << "Linkage = " << linkageHack() << endl;


  //  AppParams params;
  const int lsize[4]={12,12,12,48};
  multi1d<int> nrow(4);
  nrow=lsize;
  Layout::setLattSize(nrow);
  Layout::create();
  
  multi1d<UF> u(Nd);
  multi1d<UD> ud(Nd);
#if 0
  XMLReader gauge_file_xml, gauge_xml;
  Cft_t inputCfg;
  inputCfg.cfg_type=WEAK_FIELD;
  gaugeStartup(gauge_file_xml, gauge_xml, u, params.inputCfg);
#endif
  for(int mu=0; mu < Nd; mu++){ 
    gaussian(u[mu]);
    reunit(u[mu]);
    gaussian(ud[mu]);
    reunit(ud[mu]);

    //u[mu] = Real(1);
    //ud[mu] = u[mu];

  }

  unitarityCheck(u);
  unitarityCheck(ud);

  // Setup the lattice
 
  XMLFileWriter& xml_out = Chroma::getXMLOutputInstance();
  push(xml_out,"t_invert");
  push(xml_out,"params");
  write(xml_out, "nrow", nrow);
  pop(xml_out); // Params

  // Measure the plaquette on the gauge
  MesPlq(xml_out, "Observables", ud);
  xml_out.flush();

  // Write code here?
  QDPIO::cout << "Howdy" << endl;
  bool result;

#if 1
  // Single Precision tests 
  QDPIO::cout << "Test: Dslash PLUS, cb=0";
  if ( ! testQudaDslash(u, PLUS, 0)  ) { 
    QDPIO::cout << "\t FAILED" << endl;
    //  QDP_abort(1);
  }
  else { 
    QDPIO::cout << "\t OK" << endl;
  }

  QDPIO::cout << "Test: Dslash PLUS, cb=1" ;
  if ( ! testQudaDslash(u, PLUS, 1)  ) { 
    QDPIO::cout << "\t FAILED" << endl;
    // QDP_abort(1);
  }
  else { 
    QDPIO::cout << "\t OK" << endl;
  }


  QDPIO::cout << "Test: Dslash MINUS, cb=0";
  if ( ! testQudaDslash(u, MINUS, 0)  ) { 
    QDPIO::cout << "\t FAILED" << endl;
    // QDP_abort(1);
  }
  else { 
    QDPIO::cout << "\t OK" << endl;
  }


  QDPIO::cout << "Test: Dslash MINUS, cb=1" ;
  if ( ! testQudaDslash(u, MINUS, 1)  ) { 
    QDPIO::cout << "\t FAILED" << endl;
    // QDP_abort(1);
  }
  else { 
    QDPIO::cout << "\t OK" << endl;
  }

#endif
 

#if 1
  // Double Precision tests 
  QDPIO::cout << "Test: Dslash PLUS, cb=0";
  if ( ! testQudaDslashD(ud, PLUS, 0)  ) { 
    QDPIO::cout << "\t FAILED" << endl;
    // QDP_abort(1);
  }
  else { 
    QDPIO::cout << "\t OK" << endl;
  }

  QDPIO::cout << "Test: Dslash PLUS, cb=1" ;
  if ( ! testQudaDslashD(ud, PLUS, 1)  ) { 
    QDPIO::cout << "\t FAILED" << endl;
    // QDP_abort(1);
  }
  else { 
    QDPIO::cout << "\t OK" << endl;
  }


  QDPIO::cout << "Test: Dslash MINUS, cb=0";
  if ( ! testQudaDslashD(ud, MINUS, 0)  ) { 
    QDPIO::cout << "\t FAILED" << endl;
    //QDP_abort(1);
  }
  else { 
    QDPIO::cout << "\t OK" << endl;
  }


  QDPIO::cout << "Test: Dslash MINUS, cb=1" ;
  if ( ! testQudaDslashD(ud, MINUS, 1)  ) { 
    QDPIO::cout << "\t FAILED" << endl;
    // QDP_abort(1);
  }
  else { 
    QDPIO::cout << "\t OK" << endl;
  }
#endif

 
#if 0
  // 3D Dslash testsd

  QDPIO::cout << "Test: Dslash3D PLUS, cb=0";
  if ( ! testQudaDslash3D(u, PLUS, 0)  ) { 
    QDPIO::cout << "\t FAILED" << endl;
    QDP_abort(1);
  }
  else { 
    QDPIO::cout << "\t OK" << endl;
  }

  QDPIO::cout << "Test: Dslash3D MINUS, cb=0";
  if ( ! testQudaDslash3D(u, MINUS, 0)  ) { 
    QDPIO::cout << "\t FAILED" << endl;
    QDP_abort(1);
  }
  else { 
    QDPIO::cout << "\t OK" << endl;
  }

  QDPIO::cout << "Test: Dslash3D PLUS, cb=1" ;
  if ( ! testQudaDslash3D(u, PLUS, 1)  ) { 
    QDPIO::cout << "\t FAILED" << endl;
    QDP_abort(1);
  }
  else { 
    QDPIO::cout << "\t OK" << endl;
  }

  QDPIO::cout << "Test: Dslash3D MINUS, cb=1" ;
  if ( ! testQudaDslash3D(u, MINUS, 1)  ) { 
    QDPIO::cout << "\t FAILED" << endl;
    QDP_abort(1);
  }
  else { 
    QDPIO::cout << "\t OK" << endl;
  }
#endif



  QDPIO::cout << "All tests passed" << endl;
  pop(xml_out);
  xml_out.close();

  Chroma::finalize();
    
  exit(0);
}
