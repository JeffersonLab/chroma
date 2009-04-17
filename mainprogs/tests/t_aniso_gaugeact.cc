#include "chroma.h"
#include "actions/gauge/gaugebcs/periodic_gaugebc.h"
#include "actions/gauge/gaugestates/simple_gaugestate.h"
#include "actions/gauge/gaugeacts/plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/spatial_two_plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/plaq_plus_spatial_two_plaq_gaugeact.h"
#include "actions/gauge/gaugeacts/aniso_spectrum_gaugeact.h"

using namespace QDP;

typedef multi1d<LatticeColorMatrix> U;

int main(int argc, char **argv)
{ 
// Initialise QDP
  Chroma::initialize(&argc, &argv);

  multi1d<int> nrow(Nd);

  nrow[0] = nrow[1] = nrow[2] = 4;
  nrow[3] = 8;

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

  Handle< GaugeBC< U,U> > gbc = new PeriodicGaugeBC<U,U>();
  Handle< CreateGaugeState<U, U> > cgs = new CreateSimpleGaugeState<U,U>(gbc);
  

  AnisoParam_t aniso; 
  aniso.anisoP = true;
  aniso.t_dir = 2;
  aniso.xi_0 = 3;

  Real coeff_p_s = -0.5;
  Real coeff_p_t = -0.4;
  Real coeff_2plaq = -0.6;

  PlaqGaugeAct  S_plaq(cgs, coeff_p_s, coeff_p_t, aniso);
  SpatialTwoPlaqGaugeAct S_pp(cgs, coeff_2plaq, aniso);

  PlaqPlusSpatialTwoPlaqGaugeAct S_all(cgs, coeff_p_s, coeff_p_t, coeff_2plaq, aniso);

  Handle< GaugeState<U,U> > gstate = S_plaq.createState(u);

  Double S_all_S = S_all.S(gstate);
  Double S_old_S = S_plaq.S(gstate)+S_pp.S(gstate);

  Double norm_diff = ( S_all_S - S_old_S ) / S_old_S;

  QDPIO::cout <<"( S_new - S_old ) / S_old  = " <<  norm_diff << endl;

  U F_old(Nd);
  U F_new(Nd);
  U F_tmp(Nd);

  S_all.deriv(F_new, gstate);
  S_plaq.deriv(F_old, gstate); 
  S_pp.deriv(F_tmp, gstate);
  for(int mu=0; mu < Nd; mu++) { 
    F_old[mu] += F_tmp[mu];
  }

  for(int mu=0; mu < Nd; mu++) { 
    LatticeColorMatrix f_d = F_new[mu]- F_old[mu];

    QDPIO::cout <<" F_old = " << sqrt(norm2(F_old[mu])) << " F_new=" << sqrt(norm2(F_new[mu])) << "  (F_new-F_old) = " << sqrt(norm2(f_d)) << " (F_new - F_old) / F_old = " << sqrt(norm2(f_d)/norm2(F_old[mu]))<< endl;
  }

  
    Real u_s = Real(0.75793);
    Real u_t = Real(1);
    Real beta = Real(1.5);
    Real omega = Real(3);

    Real u_s_2 = u_s * u_s;
    Real u_s_4 = u_s_2 * u_s_2;
    Real u_s_6 = u_s_4 * u_s_2;
    Real u_s_8 = u_s_4 * u_s_4;

    // temporal powers
    Real u_t_2 = u_t * u_t;
    Real u_t_4 = u_t_2 * u_t_2;

    // Coefficients for the plaquette term (eq 4 in hep-lat/9911003)
    Real plaq_c_s = beta * Real(5) * ( Real(1) + omega ) / ( Real(3) * u_s_4 );
    Real plaq_c_t = beta * Real(4) / ( Real(3) * u_s_2 * u_t_2 );
    coeff_2plaq = Real(-5)*beta*omega/(Real(3)*u_s_8);
    Real rect_c_s = - beta / ( Real(12)*u_s_6 );
    Real rect_c_t_2 = - beta / ( Real(12)*u_s_4*u_t_2);

    // Loops that are long int the time direction ought to be ommitted
    bool no_temporal_2link = true;
    Real rect_c_t_1 = 0; // Specify a zero coefficient (skipped anyway)

    PlaqGaugeAct  S_p_1(cgs, plaq_c_s, plaq_c_t, aniso);
    SpatialTwoPlaqGaugeAct S_pp_1(cgs, coeff_2plaq, aniso);
    RectGaugeAct           S_r(cgs, rect_c_s, rect_c_t_1, rect_c_t_2, no_temporal_2link, aniso);

    AnisoSpectrumGaugeActParams g_p;
    g_p.beta = beta;
    g_p.u_s = u_s;
    g_p.u_t = u_t;
    g_p.omega = omega;
    g_p.aniso = aniso;

    AnisoSpectrumGaugeAct S_aniso(cgs, g_p);

    Double S_candidate = S_p_1.S(gstate)+S_pp_1.S(gstate)+S_r.S(gstate);
    Double S_ref = S_aniso.S(gstate);

    QDPIO::cout << " S_candidate ="<< S_candidate << " S_aniso=" << S_ref <<  " diff =" << S_candidate - S_ref << "  Rel_diff=" << (S_candidate - S_ref)/S_ref << endl;

  S_aniso.deriv(F_new, gstate);
  S_p_1.deriv(F_old, gstate); 
  S_pp_1.deriv(F_tmp, gstate);
  for(int mu=0; mu < Nd; mu++) { 
    F_old[mu] += F_tmp[mu];
  }
  S_r.deriv(F_tmp, gstate);
  for(int mu=0; mu < Nd; mu++) { 
    F_old[mu] += F_tmp[mu];
  }

  for(int mu=0; mu < Nd; mu++) { 
    LatticeColorMatrix f_d = F_new[mu]- F_old[mu];

    QDPIO::cout <<" F_old = " << sqrt(norm2(F_old[mu])) << " F_new=" << sqrt(norm2(F_new[mu])) << "  (F_new-F_old) = " << sqrt(norm2(f_d)) << " (F_new - F_old) / F_old = " << sqrt(norm2(f_d)/norm2(F_old[mu]))<< endl;
  }
    
  Chroma::finalize();

  exit(0);
}
