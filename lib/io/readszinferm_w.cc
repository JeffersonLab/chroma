// $Id: readszinferm_w.cc,v 3.1 2007-06-10 14:40:23 edwards Exp $: readszinqprop_w.cc,v 1.6 2003/04/30 21:19:33 edwards Exp $
/*!
 * @file
 * @brief  Read an old SZIN-style (checkerboarded) lattice Dirac fermion
 */

#include "chromabase.h"
#include "io/readszinferm_w.h"

#include "qdp_util.h"   // from QDP++

namespace Chroma {

//! Read a SZIN fermion. This is a simple memory dump reader.
/*!
 * \ingroup io
 *
 * \param q          lattice fermion ( Modify )
 * \param file       path ( Read )
 */    

void readSzinFerm(LatticeFermion& q, const string& file)
{
  BinaryFileReader cfg_in(file);

#if 1
  multi1d<DFermion> sitebuf_prec(Layout::vol());  // Buffer for a fermion
  multi1d<Fermion>  sitebuf(Layout::vol());

  read(cfg_in, sitebuf_prec, Layout::vol());

  /* Downcast here unfortunately Fermion s(DFermion x) doesn't work
     only DFermion s(Fermion x) */
  /* However, since I am not doing any site peeking, perhaps this will
     be OK */

  for(int site=0; site < Layout::vol(); site++) { 
	DFermion site_tmp_d = sitebuf_prec[site];
	Fermion site_tmp_r;

	for(int spin=0; spin < Ns; spin++) { 
	  DColorVector d_col_vec;
	  ColorVector col_vec;
	  
	  d_col_vec = peekSpin(site_tmp_d, spin);
          for(int color=0; color < Nc; color++) {
	     DComplex elem = peekColor(d_col_vec, color); 
	     Double elem_re = real(elem);
	     Double elem_im = imag(elem);

	     Complex elem_r = cmplx(Real(elem_re), Real(elem_im));
	     pokeColor(col_vec, elem_r, color);
          }
          pokeSpin(site_tmp_r, col_vec, spin);
        }	  
        sitebuf[site] = site_tmp_r;
  }
  // sitebuf_prec no longer needed. Free up memory
  sitebuf_prec.resize(0); 

  /* Now get out of lexicographic order */   
  multi1d<int> lattsize_cb = Layout::lattSize();
  lattsize_cb[0] /= 2;  // checkerboard in the x-direction in szin

  int index = 0; 
  for(int cb=0; cb < 2; ++cb)
  {
    for(int sitecb=0; sitecb < Layout::vol()/2; ++sitecb)
    {
      multi1d<int> coord = crtesn(sitecb, lattsize_cb);

      // construct the checkerboard offset
      int sum = 0;
      for(int m=1; m<Nd; m++) {
	sum += coord[m];
      }
      // The true lattice x-coord
      coord[0] = 2*coord[0] + ((sum + cb) & 1);

      pokeSite(q, sitebuf[index++], coord);

    }
  }

#endif

  cfg_in.close();
}

}  // end namespace Chroma
