#include "chroma.h"

#include <iostream>
#include "actions/gauge/gaugestates/simple_gaugestate.h"
#include "actions/gauge/gaugeacts/aniso_sym_spatial_gaugeact.h"
#include "actions/gauge/gaugeacts/aniso_sym_temporal_gaugeact.h"
#include "actions/gauge/gaugeacts/rect_gaugeact.h"
#include "actions/gauge/gaugeacts/plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/aniso_sym_shared_functions.h"

using namespace Chroma;


int main(int argc, char *argv[])
{
  // Initialise QDP
  Chroma::initialize(&argc, &argv);

  // Setup a small lattice
  const int nrow_arr[] = {2, 2, 2, 2};
  multi1d<int> nrow(Nd);
  nrow=nrow_arr;
  Layout::setLattSize(nrow);
  Layout::create();


  multi1d<LatticeColorMatrix> u(Nd);
  {
    XMLReader file_xml;
    XMLReader config_xml;
    Cfg_t foo; foo.cfg_type=CFG_TYPE_DISORDERED;
    // Cfg_t foo; foo.cfg_type=CFG_TYPE_SZIN; foo.cfg_file="./CFGIN";
    gaugeStartup(file_xml, config_xml, u, foo);
  }

  XMLFileWriter xml_out("./XMLDAT");
  push(xml_out, "t_forces");

  // Get Periodic Gauge Boundaries
  typedef multi1d<LatticeColorMatrix> P;
  typedef P Q;

  Handle<GaugeBC<P, Q> > gbc(new PeriodicGaugeBC);
  Real beta = Real(2.0);
  Real u_s = Real(0.78);
  Real u_t = Real(0.96);

  Real u_s_2 = u_s*u_s;
  Real u_s_4 = u_s_2*u_s_2;

  Real u_t_2 = u_t*u_t;
  Real u_t_4 = u_t_2*u_t_2;
  Real u_s_6 = u_s_2*u_s_4;

  // Coefficients for the plaquette term (eq 4 in hep-lat/9911003)
  Real plaq_c_s = beta * Real(5)  / ( Real(3) * u_s_4 );
  
  // I want to set up a plaquette and rectangle gauge action...
  Real plaq_c_t = beta * Real(4) / ( Real(3) * u_s_2 * u_t_2 );
  

  Real rect_c_s = - beta / ( Real(12)*u_s_6 );
  Real rect_c_t_2 = - beta / ( Real(12)*u_s_4*u_t_2);
  
  Real rect_c_t_1 = - beta / ( Real(12)*u_s_2*u_t_4);
  
  bool no_temporal_2link = true; // Do temporal 2 links for testing.

  AnisoParam_t aniso;
  aniso.anisoP = false;
  aniso.t_dir = 2;
  aniso.xi_0 = 2.6;
  aniso.nu = 0.8;


  Handle< CreateGaugeState<P,Q> > cgs = new CreateSimpleGaugeState<P,Q>( gbc );
  Handle<PlaqGaugeAct> plaq = new PlaqGaugeAct(cgs, plaq_c_s, plaq_c_t, aniso);
  Handle<RectGaugeAct>  rect = new RectGaugeAct(cgs, rect_c_s, rect_c_t_1, rect_c_t_2, no_temporal_2link, aniso);

  // Create a gauge state
  Handle< GaugeState<P,Q> > gs = plaq->createState(u);

  // Get the force from the plaquette
  P ds_plaq; 
  P ds_rect;

  plaq->deriv(ds_plaq, gs);
  rect->deriv(ds_rect, gs);

  for(int mu=0; mu < Nd; mu++) { 
    ds_plaq[mu] += ds_rect[mu];
  }

  P ds_u;


  ds_u.resize(Nd);
  ds_u = zero;

  
  Real c_munu_plaq;
  Real c_munu_rect;

  const Q& u_bc = gs->getLinks();

  Real xi_0;

  if( aniso.anisoP == true ) { 
    xi_0 = aniso.xi_0;
  }
  else { 
    xi_0 = Real(1);
  }

  multi1d<LatticeColorMatrix> ds_tmp(Nd);
  ds_tmp = zero;
  

  for(int mu = 0; mu < Nd; mu++) { 
    for(int nu = 0 ; nu < Nd; nu++) { 
      if( mu == nu ) continue;
      


      if( mu != aniso.t_dir && nu != aniso.t_dir ) { 
	// both mu and nu are spatial
	c_munu_plaq = plaq_c_s/xi_0;
	c_munu_rect = rect_c_s/xi_0;
      }
      else if( mu == aniso.t_dir ) { 
	// nu is spatial mu is temporal (2 link in time dir)
	c_munu_plaq = plaq_c_t * xi_0;
	c_munu_rect = rect_c_t_1 * xi_0;
      }
      else {
	// mu is spatial nu is temporal (1 link in time dir)
	c_munu_plaq = plaq_c_t * xi_0;
	c_munu_rect = rect_c_t_2 * xi_0;
      }


      AnisoSym::deriv_part(mu, 
			   nu, 
			   aniso.t_dir,
			   c_munu_plaq, 
			   c_munu_rect, 
			   no_temporal_2link,
			   ds_tmp,
			   u_bc);

      
    }
  }

  for(int mu=0; mu < Nd; mu++) { 
      ds_u[mu] += u_bc[mu]*ds_tmp[mu];
  }


  P diff_plaq;
  P diff_rect;

  diff_plaq.resize(Nd);
  diff_rect.resize(Nd);

  for(int mu=0; mu < Nd; mu++) {
    diff_plaq[mu] = ds_plaq[mu] - ds_u[mu];
  }
  
  write(xml_out, "DIFF_PLAQ", diff_plaq);

  for(int mu=0; mu < Nd; mu++) {
    QDPIO::cout << "Diff Plaq["<<mu<<"] = " << norm2(diff_plaq[mu]) << endl;
  }


  LatticeReal lgimp=zero;
  
  // Check on actions:
  for(int mu = 0; mu < Nd; mu++) { 
    for(int nu = 0 ; nu < Nd; nu++) { 
      if( mu == nu ) continue;
      


      if( mu != aniso.t_dir && nu != aniso.t_dir ) { 
	// both mu and nu are spatial
	c_munu_plaq = plaq_c_s/xi_0;
	c_munu_rect = rect_c_s/xi_0;
      }
      else if( mu == aniso.t_dir ) { 
	// nu is spatial mu is temporal (2 link in time dir)
	c_munu_plaq = plaq_c_t * xi_0;
	c_munu_rect = rect_c_t_1 * xi_0;
      }
      else {
	// mu is spatial nu is temporal (1 link in time dir)
	c_munu_plaq = plaq_c_t *  xi_0;
	c_munu_rect = rect_c_t_2 *  xi_0;
      }
      

      AnisoSym::S_part(mu, 
		       nu, 
		       aniso.t_dir,
		       c_munu_plaq, 
		       c_munu_rect, 
		       no_temporal_2link,
		       lgimp,
		       u_bc);

    }
  }
  
  Double s_manual = sum(lgimp);
  s_manual *= -Double(1)/Double(Nc);

  Double s_old = plaq->S(gs) + rect->S(gs);
  

  QDPIO::cout << "Manual S=  " << s_manual << endl;
  QDPIO::cout << "New    S=  " << s_old << endl;
  QDPIO::cout << "Diff    =  " << s_old - s_manual << endl;
    // Finish
  Chroma::finalize();

  exit(0);
}

